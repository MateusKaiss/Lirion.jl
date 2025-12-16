function init_dc_bus_partition!(
    part::PartitionDC,
    sys_zone,
    optimizer_factory,
    rho_theta::Float64,
)
    prepare_hotpaths!(part)
    m = build_model(sys_zone, DC; optimizer = optimizer_factory)
    part.model = m
    part.obj_base = m[:generation_cost]
    part.theta_var = m[:theta]
    part.objective_params = sys_zone.objective
    part.has_solution = false

    _set_model_silent!(m)
    for I in keys(part.theta_var)
        set_start_value(part.theta_var[I], 0.0)
    end

    _build_quad_bus_expr!(part, rho_theta)
    @objective(m, Min, part.obj_base + part.quad_bus_expr)

    return nothing
end

# BRANCH/DC init (p_flow on cut edges only)
function init_dc_branch_partition!(
    part::PartitionDC,
    sys_zone,
    optimizer_factory,
    rhoe::Float64,
    cut_edges,
)
    prepare_hotpaths!(part)
    part.zone_edges = _zone_cut_edges(part, cut_edges)

    m = build_model(sys_zone, DC; optimizer = optimizer_factory)
    part.model = m
    part.obj_base = m[:generation_cost]
    part.theta_var = m[:theta]
    part.objective_params = sys_zone.objective
    part.has_solution = false

    _set_model_silent!(m)
    for I in keys(part.theta_var)
        set_start_value(part.theta_var[I], 0.0)
    end

    _build_quad_branch_expr!(part, rhoe, part.zone_edges)
    @objective(m, Min, part.obj_base + part.quad_branch_expr)

    return nothing
end

function init_dc_hybrid_partition!(
    part::PartitionDC,
    sys_zone,
    optimizer_factory,
    rho_theta::Float64,
    rhoe::Float64,
    cut_edges,
)
    prepare_hotpaths!(part)
    part.zone_edges = _zone_cut_edges(part, cut_edges)

    m = build_model(sys_zone, DC; optimizer = optimizer_factory)
    part.model = m
    part.obj_base = m[:generation_cost]
    part.theta_var = m[:theta]
    part.objective_params = sys_zone.objective
    part.has_solution = false

    _set_model_silent!(m)
    for I in keys(part.theta_var)
        set_start_value(part.theta_var[I], 0.0)
    end

    _build_quad_bus_expr!(part, rho_theta)
    _build_quad_branch_expr!(part, rhoe, part.zone_edges)

    @objective(m, Min, part.obj_base + part.quad_bus_expr + part.quad_branch_expr)

    return nothing
end

function run_admm(config::DistributedOPFModel, sys::ADMMSystem)
    t0_ns = time_ns()
    iter_done = 0

    @assert config.base_model == DC "ADMM bus/DC requer base_model == DC."

    partitions = build_dc_partitions(sys)
    consensus_buses = sys.shared_buses
    consensus_buses_sorted, consensus_buses_idx = prepare_dense_maps(consensus_buses)

    rho = Ref(config.rho_theta)
    eps_abs = config.eps_primal
    eps_rel = config.eps_dual
    adaptive = _adaptive_on(config)

    optimizer_factory = resolve_optimizer_dc(config)

    for (zone_id, part) in partitions
        init_dc_bus_partition!(
            part,
            sys.parts[zone_id],
            optimizer_factory,
            rho[],
        )
    end

    zone_ids = collect(keys(partitions))

    duals_dense = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(consensus_buses_sorted)) for z in zone_ids
    )
    lambda_zone_buf = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(consensus_buses_sorted)) for z in zone_ids
    )

    theta_consensus_vec = zeros(Float64, length(consensus_buses_sorted))
    theta_consensus_prev_vec = similar(theta_consensus_vec)

    res_primal_vec = Float64[]
    res_dual_vec = Float64[]
    eps_primal_vec = Float64[]
    eps_dual_vec = Float64[]
    rho_vec = Float64[]
    theta_avg_history = Dict{Int, Vector{Float64}}()
    theta_local_history = Dict{Int, Dict{Int, Vector{Float64}}}()

    prev_rho = rho[]

    for k in 1:config.max_iter
        _FAST_USE_LOGS && println("\nIteração $k (bus/DC)")

        if rho[] != prev_rho
            for (_, part) in partitions
                _build_quad_bus_expr!(part, rho[])
            end
            prev_rho = rho[]
        end

        nz = length(zone_ids)
        results = Vector{Any}(undef, nz)

        Threads.@threads for i in 1:nz
            z = zone_ids[i]
            part = partitions[z]

            lambdaz = gather_zone_lambda!(lambda_zone_buf[z], duals_dense[z])

            lin = _build_linear_bus_expr(
                part,
                consensus_buses_sorted,
                consensus_buses_idx,
                theta_consensus_vec,
                lambdaz,
                rho[],
            )

            m = part.model::JuMP.Model
            @objective(m, Min, part.obj_base + part.quad_bus_expr + lin)
            optimize!(m)

            gen_cost_expr = value(m[:generation_cost])
            real_cost = compute_real_generation_cost(m, part.objective_params)

            thetas = snapshot_theta_shared!(part)
            thetasnap = Dict{Int, Float64}()
            lbs = part.shared_buses_local
            @inbounds for i2 in eachindex(lbs)
                thetasnap[lbs[i2]] = thetas[i2]
            end

            results[i] = (gen_cost_expr = gen_cost_expr, real_cost = real_cost, thetasnap = thetasnap)
        end

        for (i, z) in enumerate(zone_ids)
            part = partitions[z]
            out = results[i]
            part.has_solution = true
            part.generation_cost_model_expr = out.gen_cost_expr
            part.real_generation_cost = out.real_cost

            @inbounds for (b, thetab) in out.thetasnap
                push!(
                    get!(
                        get!(theta_local_history, z, Dict{Int, Vector{Float64}}()),
                        b,
                        Float64[],
                    ),
                    thetab,
                )
            end
        end

        theta_consensus_prev_vec .= theta_consensus_vec
        compute_theta_consensus_vec!(theta_consensus_vec, partitions, consensus_buses_sorted, consensus_buses_idx)
        @inbounds for b in consensus_buses_sorted
            j = consensus_buses_idx[b]
            push!(get!(theta_avg_history, b, Float64[]), theta_consensus_vec[j])
        end

        for z in zone_ids
            update_duals_theta_vec!(
                duals_dense[z],
                partitions[z],
                theta_consensus_vec,
                rho[],
                consensus_buses_sorted,
                consensus_buses_idx,
            )
        end

        res_primal = primal_residual_bus_dc(partitions, theta_consensus_vec, consensus_buses_sorted, consensus_buses_idx)
        res_dual = dual_residual_bus_dc(theta_consensus_vec, theta_consensus_prev_vec)
        push!(res_primal_vec, res_primal)
        push!(res_dual_vec, res_dual)

        eps_pri = eps_primal_bus_dc(partitions, theta_consensus_vec, consensus_buses_sorted, consensus_buses_idx, eps_abs, eps_rel)
        eps_du = eps_dual_bus_dc(theta_consensus_vec, theta_consensus_prev_vec, eps_abs, eps_rel)
        push!(eps_primal_vec, eps_pri)
        push!(eps_dual_vec, eps_du)

        _FAST_USE_LOGS && println("‖r‖ = $res_primal | eps_pri = $eps_pri")
        _FAST_USE_LOGS && println("‖s‖ = $res_dual   | eps_du  = $eps_du")

        push!(rho_vec, rho[])

        iter_done = k
        if res_primal < eps_pri && res_dual < eps_du
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        if adaptive
            adaptive_rho_dc!(
                rho,
                res_primal,
                res_dual;
                μ_incr = config.mu_incr,
                μ_decr = config.mu_decr,
                τ_incr = config.tau_incr,
                τ_decr = config.tau_decr,
            )
        end
    end

    save_dc_bus_plots(
        res_primal_vec,
        res_dual_vec,
        eps_primal_vec,
        eps_dual_vec,
        rho_vec,
        theta_avg_history;
        stem = "admm",
    )

    total_real_cost = sum(part.real_generation_cost for part in values(partitions))
    total_expr_cost = sum(part.generation_cost_model_expr for part in values(partitions))

    println("\nResumo (bus/DC):")
    println("→ Total (value(:generation_cost)) = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("rho_theta inicial: $(config.rho_theta) | rho_theta final: $(rho[])")
    elapsed_s = (time_ns() - t0_ns) / 1e9
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (bus/DC): $(_pretty_time(elapsed_s))")

    return total_real_cost, partitions, iter_done, elapsed_s
end


#################################
# ADMM: BRANCH/DC               #
#################################

function run_admm_branch(
    config::DistributedOPFModel,
    sys::ADMMSystem;
    use_legacy_rho::Bool = false,
    μ_incr::Float64 = 2.0,
    μ_decr::Float64 = 2.0,
    τ_incr::Float64 = 1.05,
    τ_decr::Float64 = 0.98,
)
    t0_ns = time_ns()
    iter_done = 0

    @assert config.base_model == DC "ADMM-branch/DC requer base_model == DC."

    partitions = build_dc_partitions(sys)
    cut_edges = sys.cut_edges

    rho = Ref(config.rho_flow)
    eps_abs = config.eps_primal
    eps_rel = config.eps_dual
    adaptive = _adaptive_on(config)

    # usa o que está no config (mantém seu comportamento original)
    μ_incr = config.mu_incr
    μ_decr = config.mu_decr
    τ_incr = config.tau_incr
    τ_decr = config.tau_decr

    z_flow = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    z_flow_prev = Dict{Int, Float64}()
    duals = Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))

    res_primal_vec = Float64[]
    res_dual_vec = Float64[]
    eps_primal_vec = Float64[]
    eps_dual_vec = Float64[]
    rho_vec = Float64[]

    optimizer_factory = resolve_optimizer_dc(config)

    # init por zona (hotpaths + zone_edges + model + quad)
    for (zone_id, part) in partitions
        init_dc_branch_partition!(
            part,
            sys.parts[zone_id],
            optimizer_factory,
            rho[],
            cut_edges,
        )
    end

    zone_ids = collect(keys(partitions))
    prev_rho = rho[]

    for k in 1:config.max_iter
        _FAST_USE_LOGS && println("\n[ADMM-branch/DC] Iteração $k")

        if rho[] != prev_rho
            for (_, part) in partitions
                _build_quad_branch_expr!(part, rho[], part.zone_edges)
            end
            prev_rho = rho[]
        end

        nz = length(zone_ids)
        results = Vector{Any}(undef, nz)

        Threads.@threads for i in 1:nz
            z = zone_ids[i]
            part = partitions[z]

            lin = _build_linear_branch_expr(part, part.zone_edges, z_flow, duals[z], rho[])

            m = part.model::JuMP.Model
            @objective(m, Min, part.obj_base + part.quad_branch_expr + lin)
            optimize!(m)

            gen_cost_expr = value(m[:generation_cost])
            real_cost = compute_real_generation_cost(m, part.objective_params)

            results[i] = (gen_cost_expr = gen_cost_expr, real_cost = real_cost)
        end

        for (i, z) in enumerate(zone_ids)
            part = partitions[z]
            out = results[i]
            part.has_solution = true
            part.generation_cost_model_expr = out.gen_cost_expr
            part.real_generation_cost = out.real_cost
        end

        z_flow_prev = copy(z_flow)
        compute_z_flow_consensus!(z_flow, partitions, cut_edges)
        update_duals_flow!(duals, partitions, z_flow, cut_edges, rho[])

        res_primal = primal_residual_branch_dc(partitions, z_flow, cut_edges)
        res_dual = dual_residual_branch_dc(z_flow, z_flow_prev)
        push!(res_primal_vec, res_primal)
        push!(res_dual_vec, res_dual)

        eps_pri = eps_primal_branch_dc(partitions, z_flow, cut_edges, eps_abs, eps_rel)
        eps_du = eps_dual_branch_dc(z_flow, z_flow_prev, eps_abs, eps_rel)
        push!(eps_primal_vec, eps_pri)
        push!(eps_dual_vec, eps_du)

        _FAST_USE_LOGS && println("‖r‖ = $res_primal | eps_pri = $eps_pri")
        _FAST_USE_LOGS && println("‖s‖ = $res_dual   | eps_dual = $eps_du")

        push!(rho_vec, rho[])

        iter_done = k
        if res_primal < eps_pri && res_dual < eps_du
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        if adaptive
            if use_legacy_rho
                if res_primal > μ_incr * res_dual
                    rho[] *= τ_incr
                elseif res_dual > μ_decr * res_primal
                    rho[] *= τ_decr
                end
            else
                adaptive_rho_dc!(
                    rho,
                    res_primal,
                    res_dual;
                    μ_incr = μ_incr,
                    μ_decr = μ_decr,
                    τ_incr = τ_incr,
                    τ_decr = τ_decr,
                )
            end
        end
    end

    save_dc_branch_plots(
        res_primal_vec,
        res_dual_vec,
        eps_primal_vec,
        eps_dual_vec,
        rho_vec;
        stem = "admm_branch",
    )

    total_real_cost = sum(part.real_generation_cost for part in values(partitions))
    total_expr_cost = sum(part.generation_cost_model_expr for part in values(partitions))

    println("\n[ADMM-branch/DC] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("rhoe_final = $(rho[])")
    elapsed_s = (time_ns() - t0_ns) / 1e9
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (branch/DC): $(_pretty_time(elapsed_s))")

    return total_real_cost, partitions, iter_done, elapsed_s
end


#################################
# ADMM: HÍBRIDO (BUS + BRANCH)  #
#################################

function run_admm_hybrid(config::DistributedOPFModel, sys::ADMMSystem)
    t0_ns = time_ns()
    iter_done = 0

    @assert config.base_model == DC "Hybrid ADMM/DC requer DC."

    partitions = build_dc_partitions(sys)
    consensus_buses = sys.shared_buses
    cut_edges = sys.cut_edges

    consensus_buses_sorted, consensus_buses_idx = prepare_dense_maps(consensus_buses)

    rho_theta = Ref(config.rho_theta)
    rhoe = Ref(config.rho_flow)
    eps_abs = config.eps_primal
    eps_rel = config.eps_dual
    adaptive = _adaptive_on(config)

    theta_consensus_vec = zeros(Float64, length(consensus_buses_sorted))
    theta_consensus_prev_vec = similar(theta_consensus_vec)

    z_flow = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    z_flow_prev = Dict{Int, Float64}()

    duals_theta = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(consensus_buses_sorted)) for z in keys(partitions)
    )
    duals_flow = Dict{Int, Dict{Int, Float64}}(
        z => Dict{Int, Float64}() for z in keys(partitions)
    )

    res_primal_vec = Float64[]
    res_dual_vec = Float64[]
    eps_primal_vec = Float64[]
    eps_dual_vec = Float64[]
    rho_theta_vec = Float64[]
    rhoe_vec = Float64[]
    theta_avg_history = Dict{Int, Vector{Float64}}()

    optimizer_factory = resolve_optimizer_dc(config)

    # init por zona (hotpaths + zone_edges + model + quad)
    for (zone_id, part) in partitions
        init_dc_hybrid_partition!(
            part,
            sys.parts[zone_id],
            optimizer_factory,
            rho_theta[],
            rhoe[],
            cut_edges,
        )
    end

    zone_ids = collect(keys(partitions))
    lambdatheta_zone_buf = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(consensus_buses_sorted)) for z in zone_ids
    )

    prev_rho_theta = rho_theta[]
    prev_rhoe = rhoe[]

    for k in 1:config.max_iter
        _FAST_USE_LOGS && println("\n[ADMM-hybrid/DC] Iteração $k")

        if rho_theta[] != prev_rho_theta
            for (_, part) in partitions
                _build_quad_bus_expr!(part, rho_theta[])
            end
            prev_rho_theta = rho_theta[]
        end
        if rhoe[] != prev_rhoe
            for (_, part) in partitions
                _build_quad_branch_expr!(part, rhoe[], part.zone_edges)
            end
            prev_rhoe = rhoe[]
        end

        nz = length(zone_ids)
        results = Vector{Any}(undef, nz)

        Threads.@threads for i in 1:nz
            z = zone_ids[i]
            part = partitions[z]

            lambdatheta = gather_zone_lambda!(lambdatheta_zone_buf[z], duals_theta[z])

            lin_bus = _build_linear_bus_expr(
                part,
                consensus_buses_sorted,
                consensus_buses_idx,
                theta_consensus_vec,
                lambdatheta,
                rho_theta[],
            )
            lin_branch = _build_linear_branch_expr(
                part,
                part.zone_edges,
                z_flow,
                duals_flow[z],
                rhoe[],
            )

            m = part.model::JuMP.Model
            @objective(
                m,
                Min,
                part.obj_base + part.quad_bus_expr + part.quad_branch_expr +
                (lin_bus + lin_branch),
            )
            optimize!(m)

            gen_cost_expr = value(m[:generation_cost])
            real_cost = compute_real_generation_cost(m, part.objective_params)
            results[i] = (gen_cost_expr = gen_cost_expr, real_cost = real_cost)
        end

        for (i, z) in enumerate(zone_ids)
            part = partitions[z]
            out = results[i]
            part.has_solution = true
            part.generation_cost_model_expr = out.gen_cost_expr
            part.real_generation_cost = out.real_cost
        end

        theta_consensus_prev_vec .= theta_consensus_vec
        z_flow_prev = copy(z_flow)

        compute_theta_consensus_vec!(theta_consensus_vec, partitions, consensus_buses_sorted, consensus_buses_idx)
        @inbounds for b in consensus_buses_sorted
            j = consensus_buses_idx[b]
            push!(get!(theta_avg_history, b, Float64[]), theta_consensus_vec[j])
        end
        compute_z_flow_consensus!(z_flow, partitions, cut_edges)

        for z in zone_ids
            update_duals_theta_vec!(
                duals_theta[z],
                partitions[z],
                theta_consensus_vec,
                rho_theta[],
                consensus_buses_sorted,
                consensus_buses_idx,
            )
        end
        update_duals_flow!(duals_flow, partitions, z_flow, cut_edges, rhoe[])

        res_pri_bus = primal_residual_bus_dc(partitions, theta_consensus_vec, consensus_buses_sorted, consensus_buses_idx)
        res_pri_br = primal_residual_branch_dc(partitions, z_flow, cut_edges)
        res_primal = hypot(res_pri_bus, res_pri_br)

        res_du_bus = dual_residual_bus_dc(theta_consensus_vec, theta_consensus_prev_vec)
        res_du_br = dual_residual_branch_dc(z_flow, z_flow_prev)
        res_dual = hypot(res_du_bus, res_du_br)

        push!(res_primal_vec, res_primal)
        push!(res_dual_vec, res_dual)

        epsb_pri = eps_primal_bus_dc(partitions, theta_consensus_vec, consensus_buses_sorted, consensus_buses_idx, eps_abs, eps_rel)
        epse_pri = eps_primal_branch_dc(partitions, z_flow, cut_edges, eps_abs, eps_rel)
        push!(eps_primal_vec, hypot(epsb_pri, epse_pri))

        epsb_du = eps_dual_bus_dc(theta_consensus_vec, theta_consensus_prev_vec, eps_abs, eps_rel)
        epse_du = eps_dual_branch_dc(z_flow, z_flow_prev, eps_abs, eps_rel)
        push!(eps_dual_vec, hypot(epsb_du, epse_du))

        _FAST_USE_LOGS && println("‖r‖ = $(last(res_primal_vec)) | eps_pri = $(last(eps_primal_vec))")
        _FAST_USE_LOGS && println("‖s‖ = $(last(res_dual_vec))   | eps_du  = $(last(eps_dual_vec))")

        push!(rho_theta_vec, rho_theta[])
        push!(rhoe_vec, rhoe[])

        iter_done = k
        if res_primal < last(eps_primal_vec) && res_dual < last(eps_dual_vec)
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        if adaptive
            adaptive_rho_dc!(
                rho_theta,
                res_pri_bus,
                res_du_bus;
                μ_incr = config.mu_incr,
                μ_decr = config.mu_decr,
                τ_incr = config.tau_incr,
                τ_decr = config.tau_decr,
            )
            adaptive_rho_dc!(
                rhoe,
                res_pri_br,
                res_du_br;
                μ_incr = config.mu_incr,
                μ_decr = config.mu_decr,
                τ_incr = config.tau_incr,
                τ_decr = config.tau_decr,
            )
        end
    end

    # Plots / resumo
    plot_admm_convergence(
        res_primal_vec,
        res_dual_vec,
        eps_primal_vec,
        eps_dual_vec,
        "plot_convergencia_admm_hybrid.html",
    )

    PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(; x = 1:length(rhoe_vec), y = rhoe_vec, type = :scatter, mode = "lines+markers", name = "rhoe");
                p.layout.title = "Evolução de rhoe";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "rhoe";
                p.layout.yaxis.type = "log";
                p
            )
        )(),
        "plot_rho_flow_hybrid.html",
    )

    plot_theta_avg_heatmap(theta_avg_history, "plot_theta_avg_heatmap_hybrid.html")

    total_real_cost = sum(part.real_generation_cost for part in values(partitions))
    total_expr_cost = sum(part.generation_cost_model_expr for part in values(partitions))

    println("\n[ADMM-hybrid/DC] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("rho_theta_final = $(rho_theta[]) | rhoe_final = $(rhoe[])")
    elapsed_s = (time_ns() - t0_ns) / 1e9
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (hybrid/DC): $(_pretty_time(elapsed_s))")

    return total_real_cost, partitions, iter_done, elapsed_s
end