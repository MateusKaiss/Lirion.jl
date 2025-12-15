function init_dc_bus_partition!(
    part::PartitionDC,
    sys_zone,
    optimizer_factory,
    ρθ::Float64,
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

    _build_quad_bus_expr!(part, ρθ)
    @objective(m, Min, part.obj_base + part.quad_bus_expr)

    return nothing
end

# BRANCH/DC init (p_flow on cut edges only)
function init_dc_branch_partition!(
    part::PartitionDC,
    sys_zone,
    optimizer_factory,
    ρe::Float64,
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

    _build_quad_branch_expr!(part, ρe, part.zone_edges)
    @objective(m, Min, part.obj_base + part.quad_branch_expr)

    return nothing
end

# HYBRID/DC init (θ shared + p_flow on cut edges)
function init_dc_hybrid_partition!(
    part::PartitionDC,
    sys_zone,
    optimizer_factory,
    ρθ::Float64,
    ρe::Float64,
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

    _build_quad_bus_expr!(part, ρθ)
    _build_quad_branch_expr!(part, ρe, part.zone_edges)

    @objective(m, Min, part.obj_base + part.quad_bus_expr + part.quad_branch_expr)

    return nothing
end


############################
# ADMM: BUS/DC             #
############################

function run_admm(config::DistributedOPFModel, sys::ADMMSystem)
    t0_ns = time_ns()
    iter_done = 0

    @assert config.base_model == DC "ADMM bus/DC requer base_model == DC."

    partitions = build_dc_partitions(sys)
    shared_buses = sys.shared_buses

    # mapas densos (BUS)
    shared_buses_sorted, b2pos = prepare_dense_maps(shared_buses)

    ρ = Ref(config.rho_theta)
    ϵ_abs = config.ε_primal
    ϵ_rel = config.ε_dual
    adaptive = _adaptive_on(config)

    optimizer_factory = resolve_optimizer_dc(config)

    # init por zona (hotpaths + model + quad)
    for (zone_id, part) in partitions
        init_dc_bus_partition!(
            part,
            sys.parts[zone_id],
            optimizer_factory,
            ρ[],
        )
    end

    zone_ids = collect(keys(partitions))

    duals_dense = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in zone_ids
    )
    λ_zone_buf = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in zone_ids
    )

    θ̄_vec = zeros(Float64, length(shared_buses_sorted))
    θ̄_prev_vec = similar(θ̄_vec)

    res_primal_vec = Float64[]
    res_dual_vec = Float64[]
    eps_primal_vec = Float64[]
    eps_dual_vec = Float64[]
    rho_vec = Float64[]
    θ_avg_history = Dict{Int, Vector{Float64}}()
    θ_local_history = Dict{Int, Dict{Int, Vector{Float64}}}()

    prev_rho = ρ[]

    for k in 1:config.max_iter
        _FAST_USE_LOGS && println("\nIteração $k (bus/DC)")

        if ρ[] != prev_rho
            for (_, part) in partitions
                _build_quad_bus_expr!(part, ρ[])
            end
            prev_rho = ρ[]
        end

        nz = length(zone_ids)
        results = Vector{Any}(undef, nz)

        Threads.@threads for i in 1:nz
            z = zone_ids[i]
            part = partitions[z]

            λz = gather_zone_lambda!(λ_zone_buf[z], duals_dense[z])

            lin = _build_linear_bus_expr(
                part,
                shared_buses_sorted,
                b2pos,
                θ̄_vec,
                λz,
                ρ[],
            )

            m = part.model::JuMP.Model
            @objective(m, Min, part.obj_base + part.quad_bus_expr + lin)
            optimize!(m)

            gen_cost_expr = value(m[:generation_cost])
            real_cost = compute_real_generation_cost(m, part.objective_params)

            θs = snapshot_theta_shared!(part)
            θsnap = Dict{Int, Float64}()
            lbs = part.shared_buses_local
            @inbounds for i2 in eachindex(lbs)
                θsnap[lbs[i2]] = θs[i2]
            end

            results[i] = (gen_cost_expr = gen_cost_expr, real_cost = real_cost, θsnap = θsnap)
        end

        for (i, z) in enumerate(zone_ids)
            part = partitions[z]
            out = results[i]
            part.has_solution = true
            part.generation_cost_model_expr = out.gen_cost_expr
            part.real_generation_cost = out.real_cost

            @inbounds for (b, θb) in out.θsnap
                push!(
                    get!(
                        get!(θ_local_history, z, Dict{Int, Vector{Float64}}()),
                        b,
                        Float64[],
                    ),
                    θb,
                )
            end
        end

        θ̄_prev_vec .= θ̄_vec
        compute_theta_consensus_vec!(θ̄_vec, partitions, shared_buses_sorted, b2pos)
        @inbounds for b in shared_buses_sorted
            j = b2pos[b]
            push!(get!(θ_avg_history, b, Float64[]), θ̄_vec[j])
        end

        for z in zone_ids
            update_duals_theta_vec!(
                duals_dense[z],
                partitions[z],
                θ̄_vec,
                ρ[],
                shared_buses_sorted,
                b2pos,
            )
        end

        res_primal = primal_residual_bus_dc(partitions, θ̄_vec, shared_buses_sorted, b2pos)
        res_dual = dual_residual_bus_dc(θ̄_vec, θ̄_prev_vec)
        push!(res_primal_vec, res_primal)
        push!(res_dual_vec, res_dual)

        ϵ_pri = eps_primal_bus_dc(partitions, θ̄_vec, shared_buses_sorted, b2pos, ϵ_abs, ϵ_rel)
        ϵ_du = eps_dual_bus_dc(θ̄_vec, θ̄_prev_vec, ϵ_abs, ϵ_rel)
        push!(eps_primal_vec, ϵ_pri)
        push!(eps_dual_vec, ϵ_du)

        _FAST_USE_LOGS && println("‖r‖ = $res_primal | ε_pri = $ϵ_pri")
        _FAST_USE_LOGS && println("‖s‖ = $res_dual   | ε_du  = $ϵ_du")

        push!(rho_vec, ρ[])

        iter_done = k
        if res_primal < ϵ_pri && res_dual < ϵ_du
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        if adaptive
            adaptive_rho_dc!(
                ρ,
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
        θ_avg_history;
        stem = "admm",
    )

    total_real_cost = sum(part.real_generation_cost for part in values(partitions))
    total_expr_cost = sum(part.generation_cost_model_expr for part in values(partitions))

    println("\nResumo (bus/DC):")
    println("→ Total (value(:generation_cost)) = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("ρθ inicial: $(config.rho_theta) | ρθ final: $(ρ[])")
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

    ρ = Ref(config.rho_flow)
    ϵ_abs = config.ε_primal
    ϵ_rel = config.ε_dual
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
            ρ[],
            cut_edges,
        )
    end

    zone_ids = collect(keys(partitions))
    prev_rho = ρ[]

    for k in 1:config.max_iter
        _FAST_USE_LOGS && println("\n[ADMM-branch/DC] Iteração $k")

        if ρ[] != prev_rho
            for (_, part) in partitions
                _build_quad_branch_expr!(part, ρ[], part.zone_edges)
            end
            prev_rho = ρ[]
        end

        nz = length(zone_ids)
        results = Vector{Any}(undef, nz)

        Threads.@threads for i in 1:nz
            z = zone_ids[i]
            part = partitions[z]

            lin = _build_linear_branch_expr(part, part.zone_edges, z_flow, duals[z], ρ[])

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
        update_duals_flow!(duals, partitions, z_flow, cut_edges, ρ[])

        res_primal = primal_residual_branch_dc(partitions, z_flow, cut_edges)
        res_dual = dual_residual_branch_dc(z_flow, z_flow_prev)
        push!(res_primal_vec, res_primal)
        push!(res_dual_vec, res_dual)

        ϵ_pri = eps_primal_branch_dc(partitions, z_flow, cut_edges, ϵ_abs, ϵ_rel)
        ϵ_du = eps_dual_branch_dc(z_flow, z_flow_prev, ϵ_abs, ϵ_rel)
        push!(eps_primal_vec, ϵ_pri)
        push!(eps_dual_vec, ϵ_du)

        _FAST_USE_LOGS && println("‖r‖ = $res_primal | ε_pri = $ϵ_pri")
        _FAST_USE_LOGS && println("‖s‖ = $res_dual   | ε_dual = $ϵ_du")

        push!(rho_vec, ρ[])

        iter_done = k
        if res_primal < ϵ_pri && res_dual < ϵ_du
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        if adaptive
            if use_legacy_rho
                if res_primal > μ_incr * res_dual
                    ρ[] *= τ_incr
                elseif res_dual > μ_decr * res_primal
                    ρ[] *= τ_decr
                end
            else
                adaptive_rho_dc!(
                    ρ,
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
    println("ρe_final = $(ρ[])")
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
    shared_buses = sys.shared_buses
    cut_edges = sys.cut_edges

    shared_buses_sorted, b2pos = prepare_dense_maps(shared_buses)

    ρθ = Ref(config.rho_theta)
    ρe = Ref(config.rho_flow)
    ϵ_abs = config.ε_primal
    ϵ_rel = config.ε_dual
    adaptive = _adaptive_on(config)

    θ̄_vec = zeros(Float64, length(shared_buses_sorted))
    θ̄_prev_vec = similar(θ̄_vec)

    z_flow = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    z_flow_prev = Dict{Int, Float64}()

    duals_theta = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )
    duals_flow = Dict{Int, Dict{Int, Float64}}(
        z => Dict{Int, Float64}() for z in keys(partitions)
    )

    res_primal_vec = Float64[]
    res_dual_vec = Float64[]
    eps_primal_vec = Float64[]
    eps_dual_vec = Float64[]
    rhoθ_vec = Float64[]
    rhoe_vec = Float64[]
    θ_avg_history = Dict{Int, Vector{Float64}}()

    optimizer_factory = resolve_optimizer_dc(config)

    # init por zona (hotpaths + zone_edges + model + quad)
    for (zone_id, part) in partitions
        init_dc_hybrid_partition!(
            part,
            sys.parts[zone_id],
            optimizer_factory,
            ρθ[],
            ρe[],
            cut_edges,
        )
    end

    zone_ids = collect(keys(partitions))
    λθ_zone_buf = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in zone_ids
    )

    prev_ρθ = ρθ[]
    prev_ρe = ρe[]

    for k in 1:config.max_iter
        _FAST_USE_LOGS && println("\n[ADMM-hybrid/DC] Iteração $k")

        if ρθ[] != prev_ρθ
            for (_, part) in partitions
                _build_quad_bus_expr!(part, ρθ[])
            end
            prev_ρθ = ρθ[]
        end
        if ρe[] != prev_ρe
            for (_, part) in partitions
                _build_quad_branch_expr!(part, ρe[], part.zone_edges)
            end
            prev_ρe = ρe[]
        end

        nz = length(zone_ids)
        results = Vector{Any}(undef, nz)

        Threads.@threads for i in 1:nz
            z = zone_ids[i]
            part = partitions[z]

            λθ = gather_zone_lambda!(λθ_zone_buf[z], duals_theta[z])

            lin_bus = _build_linear_bus_expr(
                part,
                shared_buses_sorted,
                b2pos,
                θ̄_vec,
                λθ,
                ρθ[],
            )
            lin_branch = _build_linear_branch_expr(
                part,
                part.zone_edges,
                z_flow,
                duals_flow[z],
                ρe[],
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

        θ̄_prev_vec .= θ̄_vec
        z_flow_prev = copy(z_flow)

        compute_theta_consensus_vec!(θ̄_vec, partitions, shared_buses_sorted, b2pos)
        @inbounds for b in shared_buses_sorted
            j = b2pos[b]
            push!(get!(θ_avg_history, b, Float64[]), θ̄_vec[j])
        end
        compute_z_flow_consensus!(z_flow, partitions, cut_edges)

        for z in zone_ids
            update_duals_theta_vec!(
                duals_theta[z],
                partitions[z],
                θ̄_vec,
                ρθ[],
                shared_buses_sorted,
                b2pos,
            )
        end
        update_duals_flow!(duals_flow, partitions, z_flow, cut_edges, ρe[])

        res_pri_bus = primal_residual_bus_dc(partitions, θ̄_vec, shared_buses_sorted, b2pos)
        res_pri_br = primal_residual_branch_dc(partitions, z_flow, cut_edges)
        res_primal = hypot(res_pri_bus, res_pri_br)

        res_du_bus = dual_residual_bus_dc(θ̄_vec, θ̄_prev_vec)
        res_du_br = dual_residual_branch_dc(z_flow, z_flow_prev)
        res_dual = hypot(res_du_bus, res_du_br)

        push!(res_primal_vec, res_primal)
        push!(res_dual_vec, res_dual)

        ϵb_pri = eps_primal_bus_dc(partitions, θ̄_vec, shared_buses_sorted, b2pos, ϵ_abs, ϵ_rel)
        ϵe_pri = eps_primal_branch_dc(partitions, z_flow, cut_edges, ϵ_abs, ϵ_rel)
        push!(eps_primal_vec, hypot(ϵb_pri, ϵe_pri))

        ϵb_du = eps_dual_bus_dc(θ̄_vec, θ̄_prev_vec, ϵ_abs, ϵ_rel)
        ϵe_du = eps_dual_branch_dc(z_flow, z_flow_prev, ϵ_abs, ϵ_rel)
        push!(eps_dual_vec, hypot(ϵb_du, ϵe_du))

        _FAST_USE_LOGS && println("‖r‖ = $(last(res_primal_vec)) | ε_pri = $(last(eps_primal_vec))")
        _FAST_USE_LOGS && println("‖s‖ = $(last(res_dual_vec))   | ε_du  = $(last(eps_dual_vec))")

        push!(rhoθ_vec, ρθ[])
        push!(rhoe_vec, ρe[])

        iter_done = k
        if res_primal < last(eps_primal_vec) && res_dual < last(eps_dual_vec)
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        if adaptive
            adaptive_rho_dc!(
                ρθ,
                res_pri_bus,
                res_du_bus;
                μ_incr = config.mu_incr,
                μ_decr = config.mu_decr,
                τ_incr = config.tau_incr,
                τ_decr = config.tau_decr,
            )
            adaptive_rho_dc!(
                ρe,
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
                p(; x = 1:length(rhoe_vec), y = rhoe_vec, type = :scatter, mode = "lines+markers", name = "ρe");
                p.layout.title = "Evolução de ρe";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "ρe";
                p.layout.yaxis.type = "log";
                p
            )
        )(),
        "plot_rho_flow_hybrid.html",
    )

    plot_theta_avg_heatmap(θ_avg_history, "plot_theta_avg_heatmap_hybrid.html")

    total_real_cost = sum(part.real_generation_cost for part in values(partitions))
    total_expr_cost = sum(part.generation_cost_model_expr for part in values(partitions))

    println("\n[ADMM-hybrid/DC] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("ρθ_final = $(ρθ[]) | ρe_final = $(ρe[])")
    elapsed_s = (time_ns() - t0_ns) / 1e9
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (hybrid/DC): $(_pretty_time(elapsed_s))")

    return total_real_cost, partitions, iter_done, elapsed_s
end