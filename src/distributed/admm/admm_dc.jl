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

function init_dc_bus_partition!(
    partitions::Dict{Int, PartitionDC},
    sys::ADMMSystem,
    zone_ids::Vector{Int},
    optimizer_factory,
    rho_theta::Float64,
)
    for z in zone_ids
        init_dc_bus_partition!(partitions[z], sys.parts[z], optimizer_factory, rho_theta)
    end
    return nothing
end

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

function init_dc_branch_partition!(
    partitions::Dict{Int, PartitionDC},
    sys::ADMMSystem,
    zone_ids::Vector{Int},
    optimizer_factory,
    rhoe::Float64,
    cut_edges,
)
    for z in zone_ids
        init_dc_branch_partition!(
            partitions[z],
            sys.parts[z],
            optimizer_factory,
            rhoe,
            cut_edges,
        )
    end
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

function init_dc_hybrid_partition!(
    partitions::Dict{Int, PartitionDC},
    sys::ADMMSystem,
    zone_ids::Vector{Int},
    optimizer_factory,
    rho_theta::Float64,
    rhoe::Float64,
    cut_edges,
)
    for z in zone_ids
        init_dc_hybrid_partition!(
            partitions[z],
            sys.parts[z],
            optimizer_factory,
            rho_theta,
            rhoe,
            cut_edges,
        )
    end
    return nothing
end

function run_admm(config::DistributedOPFModel, sys::ADMMSystem)
    t0_ns = time_ns()

    @assert config.base_model == DC "ADMM bus/DC requer base_model == DC."

    shared_buses = sys.shared_buses
    max_iter = config.max_iter

    rho = Ref(config.rho_theta)
    eps_abs = config.eps_primal
    eps_rel = config.eps_dual
    adaptive = _adaptive_on(config)

    optimizer_factory = resolve_optimizer_dc(config)

    partitions = build_dc_partitions(sys)
    zone_ids = collect(keys(partitions))

    shared_buses_sorted, shared_buses_idx = sorted_ids_and_pos(shared_buses)
    n_shared_buses = length(shared_buses_sorted)
    n_zones = length(zone_ids)

    dual_theta = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, n_shared_buses) for z in zone_ids
    )
    dual_theta_buf = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, n_shared_buses) for z in zone_ids
    )

    theta_consensus = zeros(Float64, n_shared_buses)
    theta_consensus_prev = similar(theta_consensus)

    res_primal_hist = Float64[]
    res_dual_hist = Float64[]
    eps_primal_hist = Float64[]
    eps_dual_hist = Float64[]
    rho_hist = Float64[]

    theta_avg_history = Dict{Int, Vector{Float64}}()
    theta_local_history = Dict{Int, Dict{Int, Vector{Float64}}}()

    prev_rho = rho[]
    results = Vector{
        NamedTuple{(:gen_cost_expr, :real_cost, :thetasnap),
            Tuple{Float64, Float64, Dict{Int, Float64}}},
    }(
        undef,
        n_zones,
    )

    init_dc_bus_partition!(partitions, sys, zone_ids, optimizer_factory, rho[])

    iter_done = 0
    for iter in 1:max_iter
        iter_done = iter
        _FAST_USE_LOGS && println("\nIteração $iter (bus/DC)")

        if rho[] != prev_rho
            for part in values(partitions)
                _build_quad_bus_expr!(part, rho[])
            end
            prev_rho = rho[]
        end

        Threads.@threads for i in 1:n_zones
            z = zone_ids[i]
            part = partitions[z]
            m = part.model::JuMP.Model

            lambdaz = gather_zone_lambda!(dual_theta_buf[z], dual_theta[z])

            lin = _build_linear_bus_expr(
                part,
                shared_buses_sorted,
                shared_buses_idx,
                theta_consensus,
                lambdaz,
                rho[],
            )

            @objective(m, Min, part.obj_base + part.quad_bus_expr + lin)
            optimize!(m)

            gen_cost_expr = value(m[:generation_cost])
            real_cost = compute_real_generation_cost(m, part.objective_params)

            thetas = snapshot_theta_shared!(part)
            lbs = part.shared_buses_local
            thetasnap = Dict{Int, Float64}(b => thetas[j] for (j, b) in pairs(lbs))

            results[i] = (
                gen_cost_expr = gen_cost_expr,
                real_cost = real_cost,
                thetasnap = thetasnap,
            )
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

        theta_consensus_prev .= theta_consensus
        compute_theta_consensus_vec!(
            theta_consensus,
            partitions,
            shared_buses_sorted,
            shared_buses_idx,
        )

        @inbounds for b in shared_buses_sorted
            j = shared_buses_idx[b]
            push!(get!(theta_avg_history, b, Float64[]), theta_consensus[j])
        end

        for z in zone_ids
            update_duals_theta_vec!(
                dual_theta[z],
                partitions[z],
                theta_consensus,
                rho[],
                shared_buses_sorted,
                shared_buses_idx,
            )
        end

        r = primal_residual_bus_dc(
            partitions,
            theta_consensus,
            shared_buses_sorted,
            shared_buses_idx,
        )
        s = dual_residual_bus_dc(theta_consensus, theta_consensus_prev)

        eps_r = eps_primal_bus_dc(
            partitions,
            theta_consensus,
            shared_buses_sorted,
            shared_buses_idx,
            eps_abs,
            eps_rel,
        )
        eps_s = eps_dual_bus_dc(theta_consensus, theta_consensus_prev, eps_abs, eps_rel)

        push!(res_primal_hist, r)
        push!(res_dual_hist, s)
        push!(eps_primal_hist, eps_r)
        push!(eps_dual_hist, eps_s)
        push!(rho_hist, rho[])

        _FAST_USE_LOGS && println("‖r‖ = $r | eps_pri = $eps_r")
        _FAST_USE_LOGS && println("‖s‖ = $s | eps_du  = $eps_s")

        if r < eps_r && s < eps_s
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        if adaptive
            adaptive_rho_dc!(
                rho, r, s;
                mu_incr = config.mu_incr,
                mu_decr = config.mu_decr,
                tau_incr = config.tau_incr,
                tau_decr = config.tau_decr,
            )
        end
    end

    save_dc_bus_plots(
        res_primal_hist, res_dual_hist, eps_primal_hist, eps_dual_hist,
        rho_hist, theta_avg_history;
        stem = "admm",
    )

    total_real_cost = sum(p.real_generation_cost for p in values(partitions))
    total_expr_cost = sum(p.generation_cost_model_expr for p in values(partitions))

    elapsed_s = (time_ns() - t0_ns) / 1e9
    println("\nResumo (bus/DC):")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("rho_theta inicial: $(config.rho_theta) | rho_theta final: $(rho[])")
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (bus/DC): $(_pretty_time(elapsed_s))")

    return total_real_cost, partitions, iter_done, elapsed_s
end

function run_admm_branch(
    config::DistributedOPFModel,
    sys::ADMMSystem;
    use_legacy_rho::Bool = false,
    mu_incr::Float64 = 2.0,
    mu_decr::Float64 = 2.0,
    tau_incr::Float64 = 1.05,
    tau_decr::Float64 = 0.98,
)
    t0_ns = time_ns()

    @assert config.base_model == DC "ADMM-branch/DC requer base_model == DC."

    cut_edges = sys.cut_edges
    max_iter = config.max_iter

    rho = Ref(config.rho_flow)
    eps_abs = config.eps_primal
    eps_rel = config.eps_dual
    adaptive = _adaptive_on(config)

    # se você quer SEMPRE usar o config, sobrescreve os kwargs:
    mu_incr = config.mu_incr
    mu_decr = config.mu_decr
    tau_incr = config.tau_incr
    tau_decr = config.tau_decr

    optimizer_factory = resolve_optimizer_dc(config)

    partitions = build_dc_partitions(sys)
    zone_ids = collect(keys(partitions))
    n_zones = length(zone_ids)

    z_flow = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    z_flow_prev = Dict{Int, Float64}()

    duals = Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in zone_ids)

    res_primal_hist = Float64[]
    res_dual_hist = Float64[]
    eps_primal_hist = Float64[]
    eps_dual_hist = Float64[]
    rho_hist = Float64[]

    prev_rho = rho[]

    results = Vector{NamedTuple{(:gen_cost_expr, :real_cost),
        Tuple{Float64, Float64}}}(
        undef,
        n_zones,
    )

    init_dc_branch_partition!(
        partitions,
        sys,
        zone_ids,
        optimizer_factory,
        rho[],
        cut_edges,
    )

    iter_done = 0
    for iter in 1:max_iter
        iter_done = iter
        _FAST_USE_LOGS && println("\n[ADMM-branch/DC] Iteração $iter")

        if rho[] != prev_rho
            for part in values(partitions)
                _build_quad_branch_expr!(part, rho[], part.zone_edges)
            end
            prev_rho = rho[]
        end

        Threads.@threads for i in 1:n_zones
            z = zone_ids[i]
            part = partitions[z]
            m = part.model::JuMP.Model

            lin = _build_linear_branch_expr(part, part.zone_edges, z_flow, duals[z], rho[])

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

        r = primal_residual_branch_dc(partitions, z_flow, cut_edges)
        s = dual_residual_branch_dc(z_flow, z_flow_prev)

        eps_r = eps_primal_branch_dc(partitions, z_flow, cut_edges, eps_abs, eps_rel)
        eps_s = eps_dual_branch_dc(z_flow, z_flow_prev, eps_abs, eps_rel)

        push!(res_primal_hist, r)
        push!(res_dual_hist, s)
        push!(eps_primal_hist, eps_r)
        push!(eps_dual_hist, eps_s)
        push!(rho_hist, rho[])

        _FAST_USE_LOGS && println("‖r‖ = $r | eps_pri = $eps_r")
        _FAST_USE_LOGS && println("‖s‖ = $s | eps_du  = $eps_s")

        if r < eps_r && s < eps_s
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        if adaptive
            if use_legacy_rho
                if r > mu_incr * s
                    rho[] *= tau_incr
                elseif s > mu_decr * r
                    rho[] *= tau_decr
                end
            else
                adaptive_rho_dc!(
                    rho, r, s;
                    mu_incr = mu_incr,
                    mu_decr = mu_decr,
                    tau_incr = tau_incr,
                    tau_decr = tau_decr,
                )
            end
        end
    end

    save_dc_branch_plots(
        res_primal_hist, res_dual_hist, eps_primal_hist, eps_dual_hist,
        rho_hist; stem = "admm_branch",
    )

    total_real_cost = sum(p.real_generation_cost for p in values(partitions))
    total_expr_cost = sum(p.generation_cost_model_expr for p in values(partitions))

    elapsed_s = (time_ns() - t0_ns) / 1e9
    println("\n[ADMM-branch/DC] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("rho_flow inicial: $(config.rho_flow) | rho_flow final: $(rho[])")
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (branch/DC): $(_pretty_time(elapsed_s))")

    return total_real_cost, partitions, iter_done, elapsed_s
end

function run_admm_hybrid(config::DistributedOPFModel, sys::ADMMSystem)
    t0_ns = time_ns()

    @assert config.base_model == DC "Hybrid ADMM/DC requer DC."

    shared_buses = sys.shared_buses
    cut_edges = sys.cut_edges
    max_iter = config.max_iter

    rho_theta = Ref(config.rho_theta)
    rho_flow = Ref(config.rho_flow)
    eps_abs = config.eps_primal
    eps_rel = config.eps_dual
    adaptive = _adaptive_on(config)

    optimizer_factory = resolve_optimizer_dc(config)

    partitions = build_dc_partitions(sys)
    zone_ids = collect(keys(partitions))
    n_zones = length(zone_ids)

    shared_buses_sorted, shared_buses_idx = sorted_ids_and_pos(shared_buses)
    n_shared_buses = length(shared_buses_sorted)

    theta_consensus = zeros(Float64, n_shared_buses)
    theta_consensus_prev = similar(theta_consensus)

    z_flow = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    z_flow_prev = Dict{Int, Float64}()

    duals_theta =
        Dict{Int, Vector{Float64}}(z => zeros(Float64, n_shared_buses) for z in zone_ids)
    dual_theta_buf =
        Dict{Int, Vector{Float64}}(z => zeros(Float64, n_shared_buses) for z in zone_ids)
    duals_flow = Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in zone_ids)

    res_primal_hist = Float64[]
    res_dual_hist = Float64[]
    eps_primal_hist = Float64[]
    eps_dual_hist = Float64[]
    rho_theta_hist = Float64[]
    rho_flow_hist = Float64[]

    theta_avg_history = Dict{Int, Vector{Float64}}()

    prev_rho_theta = rho_theta[]
    prev_rho_flow = rho_flow[]

    results = Vector{NamedTuple{(:gen_cost_expr, :real_cost),
        Tuple{Float64, Float64}}}(
        undef,
        n_zones,
    )

    init_dc_hybrid_partition!(
        partitions,
        sys,
        zone_ids,
        optimizer_factory,
        rho_theta[],
        rho_flow[],
        cut_edges,
    )

    iter_done = 0
    for iter in 1:max_iter
        iter_done = iter
        _FAST_USE_LOGS && println("\n[ADMM-hybrid/DC] Iteração $iter")

        if rho_theta[] != prev_rho_theta
            for part in values(partitions)
                _build_quad_bus_expr!(part, rho_theta[])
            end
            prev_rho_theta = rho_theta[]
        end
        if rho_flow[] != prev_rho_flow
            for part in values(partitions)
                _build_quad_branch_expr!(part, rho_flow[], part.zone_edges)
            end
            prev_rho_flow = rho_flow[]
        end

        Threads.@threads for i in 1:n_zones
            z = zone_ids[i]
            part = partitions[z]
            m = part.model::JuMP.Model

            lambdatheta = gather_zone_lambda!(dual_theta_buf[z], duals_theta[z])

            lin_bus = _build_linear_bus_expr(
                part,
                shared_buses_sorted,
                shared_buses_idx,
                theta_consensus,
                lambdatheta,
                rho_theta[],
            )
            lin_branch = _build_linear_branch_expr(
                part,
                part.zone_edges,
                z_flow,
                duals_flow[z],
                rho_flow[],
            )

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

        theta_consensus_prev .= theta_consensus
        z_flow_prev = copy(z_flow)

        compute_theta_consensus_vec!(
            theta_consensus,
            partitions,
            shared_buses_sorted,
            shared_buses_idx,
        )
        @inbounds for b in shared_buses_sorted
            j = shared_buses_idx[b]
            push!(get!(theta_avg_history, b, Float64[]), theta_consensus[j])
        end
        compute_z_flow_consensus!(z_flow, partitions, cut_edges)

        for z in zone_ids
            update_duals_theta_vec!(
                duals_theta[z],
                partitions[z],
                theta_consensus,
                rho_theta[],
                shared_buses_sorted,
                shared_buses_idx,
            )
        end
        update_duals_flow!(duals_flow, partitions, z_flow, cut_edges, rho_flow[])

        r_bus = primal_residual_bus_dc(
            partitions,
            theta_consensus,
            shared_buses_sorted,
            shared_buses_idx,
        )
        r_br = primal_residual_branch_dc(partitions, z_flow, cut_edges)
        r_tot = hypot(r_bus, r_br)

        s_bus = dual_residual_bus_dc(theta_consensus, theta_consensus_prev)
        s_br = dual_residual_branch_dc(z_flow, z_flow_prev)
        s_tot = hypot(s_bus, s_br)

        eps_bus_r = eps_primal_bus_dc(
            partitions,
            theta_consensus,
            shared_buses_sorted,
            shared_buses_idx,
            eps_abs,
            eps_rel,
        )
        eps_br_r = eps_primal_branch_dc(partitions, z_flow, cut_edges, eps_abs, eps_rel)
        eps_r_tot = hypot(eps_bus_r, eps_br_r)

        eps_bus_s = eps_dual_bus_dc(theta_consensus, theta_consensus_prev, eps_abs, eps_rel)
        eps_br_s = eps_dual_branch_dc(z_flow, z_flow_prev, eps_abs, eps_rel)
        eps_s_tot = hypot(eps_bus_s, eps_br_s)

        push!(res_primal_hist, r_tot)
        push!(res_dual_hist, s_tot)
        push!(eps_primal_hist, eps_r_tot)
        push!(eps_dual_hist, eps_s_tot)
        push!(rho_theta_hist, rho_theta[])
        push!(rho_flow_hist, rho_flow[])

        _FAST_USE_LOGS && println("‖r‖ = $r_tot | eps_pri = $eps_r_tot")
        _FAST_USE_LOGS && println("‖s‖ = $s_tot | eps_du  = $eps_s_tot")

        if r_tot < eps_r_tot && s_tot < eps_s_tot
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        if adaptive
            adaptive_rho_dc!(
                rho_theta, r_bus, s_bus;
                mu_incr = config.mu_incr,
                mu_decr = config.mu_decr,
                tau_incr = config.tau_incr,
                tau_decr = config.tau_decr,
            )
            adaptive_rho_dc!(
                rho_flow, r_br, s_br;
                mu_incr = config.mu_incr,
                mu_decr = config.mu_decr,
                tau_incr = config.tau_incr,
                tau_decr = config.tau_decr,
            )
        end
    end

    plot_admm_convergence(
        res_primal_hist, res_dual_hist, eps_primal_hist, eps_dual_hist,
        "plot_convergencia_admm_hybrid.html",
    )

    PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(; x = 1:length(rho_flow_hist), y = rho_flow_hist,
                    type = :scatter, mode = "lines+markers", name = "rho_flow");
                p.layout.title = "Evolução de rho_flow";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "rho_flow";
                p.layout.yaxis.type = "log";
                p
            )
        )(),
        "plot_rho_flow_hybrid.html",
    )

    plot_theta_avg_heatmap(theta_avg_history, "plot_theta_avg_heatmap_hybrid.html")

    total_real_cost = sum(p.real_generation_cost for p in values(partitions))
    total_expr_cost = sum(p.generation_cost_model_expr for p in values(partitions))

    elapsed_s = (time_ns() - t0_ns) / 1e9
    println("\n[ADMM-hybrid/DC] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("rho_theta inicial: $(config.rho_theta) | rho_theta final: $(rho_theta[])")
    println("rho_flow  inicial: $(config.rho_flow)  | rho_flow  final: $(rho_flow[])")
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (hybrid/DC): $(_pretty_time(elapsed_s))")

    return total_real_cost, partitions, iter_done, elapsed_s
end
