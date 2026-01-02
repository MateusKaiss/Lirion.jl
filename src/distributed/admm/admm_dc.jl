# ============================================================
# ADMM-DC
# ============================================================

# ------------------------------------------------------------
# Partition initialization
# ------------------------------------------------------------

# ---------- dc bus partition ----------

function init_dc_bus_partition!(
    part::PartitionDC,
    sys_zone::SystemData,
    optimizer_factory,
    rho_theta::Float64,
)
    prepare_hotpaths!(part)
    model = build_model(sys_zone, DC; optimizer = optimizer_factory)
    part.model = model
    part.obj_base = model[:generation_cost]
    part.theta_var = model[:theta]
    part.objective_params = sys_zone.objective
    part.has_solution = false

    for I in keys(part.theta_var)
        set_start_value(part.theta_var[I], 0.0)
    end

    _build_quad_bus_expr!(part, rho_theta)
    @objective(model, Min, part.obj_base + part.quad_bus_expr)

    return nothing
end

function init_dc_bus_partition!(
    partitions::Partitions{PartitionDC},
    sys::ADMMSystem,
    optimizer_factory,
    rho_theta::Float64,
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n "ADMMSystem.parts must be aligned with partitions.parts"

    @inbounds for i in 1:n
        init_dc_bus_partition!(
            partitions.parts[i],
            sys.parts[i],
            optimizer_factory,
            rho_theta,
        )
    end
    return nothing
end

# ---------- dc branch partition ----------

function init_dc_branch_partition!(
    part::PartitionDC,
    sys_zone::SystemData,
    optimizer_factory,
    rho::Float64,
    interzonal_branches::Vector{Int},
)
    prepare_hotpaths!(part)
    part.zone_edges = _zone_cut_edges(part, interzonal_branches)

    model = build_model(sys_zone, DC; optimizer = optimizer_factory)
    part.model = model
    part.obj_base = model[:generation_cost]
    part.theta_var = model[:theta]
    part.objective_params = sys_zone.objective
    part.has_solution = false

    for I in keys(part.theta_var)
        set_start_value(part.theta_var[I], 0.0)
    end

    _build_quad_branch_expr!(part, rho, part.zone_edges)
    @objective(model, Min, part.obj_base + part.quad_branch_expr)

    return nothing
end

function init_dc_branch_partition!(
    partitions::Partitions{PartitionDC},
    sys::ADMMSystem,
    optimizer_factory,
    rho::Float64,
    interzonal_branches::Vector{Int},
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n "ADMMSystem.parts must be aligned with partitions.parts"

    @inbounds for i in 1:n
        init_dc_branch_partition!(
            partitions.parts[i],
            sys.parts[i],
            optimizer_factory,
            rho,
            interzonal_branches,
        )
    end
    return nothing
end

# ---------- dc hybrid partition ----------

function init_dc_hybrid_partition!(
    part::PartitionDC,
    sys_zone::SystemData,
    optimizer_factory,
    rho_theta::Float64,
    rho::Float64,
    interzonal_branches::Vector{Int},
)
    prepare_hotpaths!(part)
    part.zone_edges = _zone_cut_edges(part, interzonal_branches)

    model = build_model(sys_zone, DC; optimizer = optimizer_factory)
    part.model = model
    part.obj_base = model[:generation_cost]
    part.theta_var = model[:theta]
    part.objective_params = sys_zone.objective
    part.has_solution = false

    for I in keys(part.theta_var)
        set_start_value(part.theta_var[I], 0.0)
    end

    _build_quad_bus_expr!(part, rho_theta)
    _build_quad_branch_expr!(part, rho, part.zone_edges)

    @objective(model, Min, part.obj_base + part.quad_bus_expr + part.quad_branch_expr)

    return nothing
end

function init_dc_hybrid_partition!(
    partitions::Partitions{PartitionDC},
    sys::ADMMSystem,
    optimizer_factory,
    rho_theta::Float64,
    rho::Float64,
    interzonal_branches::Vector{Int},
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n "ADMMSystem.parts must be aligned with partitions.parts"

    @inbounds for i in 1:n
        init_dc_hybrid_partition!(
            partitions.parts[i],
            sys.parts[i],
            optimizer_factory,
            rho_theta,
            rho,
            interzonal_branches,
        )
    end
    return nothing
end

# ============================================================
# RUNNERS (ADMM-DC)
# ============================================================

# ------------------------------------------------------------
# Runner: ADMM bus/DC
# ------------------------------------------------------------

function run_admm(
    sys::ADMMSystem,
    config::DistributedOPFModel,
    ::ADMM{Bus, DCTag},
)::LirionRunResult
    t0_ns = time_ns()
    @assert config.base_model == DC "ADMM bus/DC requer base_model == DC."

    shared_buses = sys.shared_buses
    max_iter = config.max_iter

    rho = Ref(config.rho_theta)
    prev_rho = Ref(rho[])
    eps_abs = config.eps_abs
    eps_rel = config.eps_rel
    adaptive = _adaptive_on(config)

    optimizer_factory = resolve_optimizer(config)

    partitions = build_dc_partitions(sys)
    zone_ids = partitions.zones
    n_zones = length(zone_ids)

    shared_buses_sorted, shared_buses_to_pos = sorted_ids_and_id2idx(shared_buses)
    n_shared_buses = length(shared_buses_sorted)

    dual_theta = Dict{Int, Vector{Float64}}(
        zone_id => zeros(Float64, n_shared_buses) for zone_id in zone_ids
    )
    dual_theta_buf = Dict{Int, Vector{Float64}}(
        zone_id => zeros(Float64, n_shared_buses) for zone_id in zone_ids
    )

    theta_consensus = zeros(Float64, n_shared_buses)
    theta_consensus_prev = similar(theta_consensus)

    hist = ADMMHistory()
    hist.meta[:kind] = :dc_bus
    hist.meta[:tag] = "admm_dc_bus"

    theta_avg_history = Dict{Int, Vector{Float64}}()
    theta_local_history = Dict{Int, Dict{Int, Vector{Float64}}}()

    results = Vector{
        NamedTuple{(:gen_cost_expr, :real_cost, :theta_snap),
            Tuple{Float64, Float64, Dict{Int, Float64}}},
    }(
        undef,
        n_zones,
    )

    init_dc_bus_partition!(partitions, sys, optimizer_factory, rho[])

    iter_done = 0
    for iter in 1:max_iter
        iter_done = iter
        _FAST_USE_LOGS && println("\nIteração $iter (bus/DC)")

        refresh_rho_bus!(partitions, rho, prev_rho)

        Threads.@threads for i in 1:n_zones
            zone_id = zone_ids[i]
            part = get_partition(partitions, zone_id)
            model = part.model::JuMP.Model

            theta_dual_buf_z = dual_theta_buf[zone_id]
            theta_dual_buf_z .= dual_theta[zone_id]

            lin = _build_linear_bus_expr(
                part,
                shared_buses_sorted,
                shared_buses_to_pos,
                theta_consensus,
                theta_dual_buf_z,
                rho[],
            )

            @objective(model, Min, part.obj_base + part.quad_bus_expr + lin)
            optimize!(model)

            gen_cost_expr = value(model[:generation_cost])
            real_cost = compute_real_generation_cost(model, part.objective_params)

            thetas = snapshot_theta_shared!(part)
            shared_bus_ids = part.shared_buses_local
            theta_snap =
                Dict{Int, Float64}(b => thetas[j] for (j, b) in pairs(shared_bus_ids))

            results[i] = (
                gen_cost_expr = gen_cost_expr,
                real_cost = real_cost,
                theta_snap = theta_snap,
            )
        end

        for (i, zone_id) in enumerate(zone_ids)
            part = get_partition(partitions, zone_id)
            out = results[i]
            part.has_solution = true
            part.generation_cost_model_expr = out.gen_cost_expr
            part.real_generation_cost = out.real_cost

            @inbounds for (b, thetab) in out.theta_snap
                push!(
                    get!(
                        get!(theta_local_history, zone_id, Dict{Int, Vector{Float64}}()),
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
            shared_buses_to_pos,
        )

        @inbounds for b in shared_buses_sorted
            j = shared_buses_to_pos[b]
            push!(get!(theta_avg_history, b, Float64[]), theta_consensus[j])
        end

        update_duals_theta_all!(
            dual_theta,
            partitions,
            zone_ids,
            theta_consensus,
            rho,
            shared_buses_sorted,
            shared_buses_to_pos,
        )

        r = primal_residual_bus_dc(
            partitions,
            theta_consensus,
            shared_buses_sorted,
            shared_buses_to_pos,
        )
        s = dual_residual_bus_dc(theta_consensus, theta_consensus_prev)

        eps_r = eps_primal_bus_dc(
            partitions,
            theta_consensus,
            shared_buses_sorted,
            shared_buses_to_pos,
            eps_abs,
            eps_rel,
        )
        eps_s = eps_dual_bus_dc(theta_consensus, theta_consensus_prev, eps_abs, eps_rel)

        _record_iter!(hist, r, s, eps_r, eps_s; rhos = Dict(:rho_theta => rho[]))

        _FAST_USE_LOGS && println("‖r‖ = $r | eps_pri = $eps_r")
        _FAST_USE_LOGS && println("‖s‖ = $s | eps_du  = $eps_s")

        if stop_if_converged!(r, s, eps_r, eps_s)
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        maybe_adapt_rho!(adaptive, rho, r, s, config)
    end

    hist.heatmaps[:theta_avg] = theta_avg_history
    hist.heatmaps[:theta_local] = theta_local_history

    result = _finalize_admm_result(partitions, t0_ns, iter_done)
    _print_admm_summary!(
        "[ADMM-bus/DC]",
        result;
        rho_lines = ("rho_theta inicial: $(config.rho_theta) | rho_theta final: $(rho[])",),
    )

    return LirionRunResult(;
        result = result,
        model = nothing,
        system = sys,
        diagnostics = hist,
    )
end

# ------------------------------------------------------------
# Runner: ADMM branch/DC
# ------------------------------------------------------------

function run_admm(
    sys::ADMMSystem,
    config::DistributedOPFModel,
    ::ADMM{Branch, DCTag};
    use_legacy_rho::Bool = false,
)::LirionRunResult
    t0_ns = time_ns()
    @assert config.base_model == DC "ADMM-branch/DC requer base_model == DC."

    interzonal_branches = sys.interzonal_branches
    max_iter = config.max_iter

    rho = Ref(config.rho_flow_p)
    prev_rho = Ref(rho[])
    eps_abs = config.eps_abs
    eps_rel = config.eps_rel
    adaptive = _adaptive_on(config)

    optimizer_factory = resolve_optimizer(config)

    partitions = build_dc_partitions(sys)
    zone_ids = partitions.zones
    n_zones = length(zone_ids)

    flow_consensus = Dict{Int, Float64}(e => 0.0 for e in interzonal_branches)
    flow_consensus_prev = Dict{Int, Float64}()

    duals = Dict{Int, Dict{Int, Float64}}(
        zone_id => Dict{Int, Float64}() for zone_id in zone_ids
    )

    hist = ADMMHistory()
    hist.meta[:kind] = :dc_branch
    hist.meta[:tag] = "admm_dc_branch"

    results = Vector{NamedTuple{(:gen_cost_expr, :real_cost), Tuple{Float64, Float64}}}(
        undef,
        n_zones,
    )

    init_dc_branch_partition!(
        partitions,
        sys,
        optimizer_factory,
        rho[],
        interzonal_branches,
    )

    iter_done = 0
    for iter in 1:max_iter
        iter_done = iter
        _FAST_USE_LOGS && println("\n[ADMM-branch/DC] Iteração $iter")

        refresh_rho_branch!(partitions, rho, prev_rho)

        Threads.@threads for i in 1:n_zones
            zone_id = zone_ids[i]
            part = get_partition(partitions, zone_id)
            model = part.model::JuMP.Model

            lin = _build_linear_branch_expr(
                part,
                part.zone_edges,
                flow_consensus,
                duals[zone_id],
                rho[],
            )

            @objective(model, Min, part.obj_base + part.quad_branch_expr + lin)
            optimize!(model)

            results[i] = (
                gen_cost_expr = value(model[:generation_cost]),
                real_cost = compute_real_generation_cost(model, part.objective_params),
            )
        end

        for (i, zone_id) in enumerate(zone_ids)
            part = get_partition(partitions, zone_id)
            out = results[i]
            part.has_solution = true
            part.generation_cost_model_expr = out.gen_cost_expr
            part.real_generation_cost = out.real_cost
        end

        flow_consensus_prev = copy(flow_consensus)
        compute_flow_consensus!(flow_consensus, partitions, interzonal_branches)
        update_duals_flow!(duals, partitions, flow_consensus, interzonal_branches, rho[])

        r = primal_residual_branch_dc(partitions, flow_consensus, interzonal_branches)
        s = dual_residual_branch_dc(flow_consensus, flow_consensus_prev)

        eps_r = eps_primal_branch_dc(
            partitions,
            flow_consensus,
            interzonal_branches,
            eps_abs,
            eps_rel,
        )
        eps_s = eps_dual_branch_dc(flow_consensus, flow_consensus_prev, eps_abs, eps_rel)

        _record_iter!(hist, r, s, eps_r, eps_s; rhos = Dict(:rho_flow => rho[]))

        _FAST_USE_LOGS && println("‖r‖ = $r | eps_pri = $eps_r")
        _FAST_USE_LOGS && println("‖s‖ = $s | eps_du  = $eps_s")

        if stop_if_converged!(r, s, eps_r, eps_s)
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        if adaptive
            if use_legacy_rho
                if r > config.mu_incr * s
                    rho[] *= config.tau_incr
                elseif s > config.mu_decr * r
                    rho[] *= config.tau_decr
                end
            else
                adaptive_rho_dc!(
                    rho, r, s;
                    mu_incr = config.mu_incr,
                    mu_decr = config.mu_decr,
                    tau_incr = config.tau_incr,
                    tau_decr = config.tau_decr,
                )
            end
        end
    end

    result = _finalize_admm_result(partitions, t0_ns, iter_done)
    _print_admm_summary!(
        "[ADMM-branch/DC]",
        result;
        rho_lines = ("rho_flow inicial: $(config.rho_flow_p) | rho_flow final: $(rho[])",),
    )

    return LirionRunResult(;
        result = result,
        model = nothing,
        system = sys,
        diagnostics = hist,
    )
end

# ------------------------------------------------------------
# Runner: ADMM hybrid/DC
# ------------------------------------------------------------

function run_admm(
    sys::ADMMSystem,
    config::DistributedOPFModel,
    ::ADMM{Hybrid, DCTag},
)::LirionRunResult
    t0_ns = time_ns()
    @assert config.base_model == DC "Hybrid ADMM/DC requer DC."

    shared_buses = sys.shared_buses
    interzonal_branches = sys.interzonal_branches
    max_iter = config.max_iter

    rho_theta = Ref(config.rho_theta)
    rho_flow = Ref(config.rho_flow_p)

    prev_rho_theta = Ref(rho_theta[])
    prev_rho_flow = Ref(rho_flow[])

    eps_abs = config.eps_abs
    eps_rel = config.eps_rel
    adaptive = _adaptive_on(config)

    optimizer_factory = resolve_optimizer(config)

    partitions = build_dc_partitions(sys)
    zone_ids = partitions.zones
    n_zones = length(zone_ids)

    shared_buses_sorted, shared_buses_to_pos = sorted_ids_and_id2idx(shared_buses)
    n_shared_buses = length(shared_buses_sorted)

    theta_consensus = zeros(Float64, n_shared_buses)
    theta_consensus_prev = similar(theta_consensus)

    flow_consensus = Dict{Int, Float64}(e => 0.0 for e in interzonal_branches)
    flow_consensus_prev = Dict{Int, Float64}()

    duals_theta = Dict{Int, Vector{Float64}}(
        zone_id => zeros(Float64, n_shared_buses) for zone_id in zone_ids
    )
    dual_theta_buf = Dict{Int, Vector{Float64}}(
        zone_id => zeros(Float64, n_shared_buses) for zone_id in zone_ids
    )
    duals_flow = Dict{Int, Dict{Int, Float64}}(
        zone_id => Dict{Int, Float64}() for zone_id in zone_ids
    )

    hist = ADMMHistory()
    hist.meta[:kind] = :dc_hybrid
    hist.meta[:tag] = "admm_dc_hybrid"

    theta_avg_history = Dict{Int, Vector{Float64}}()

    results = Vector{NamedTuple{(:gen_cost_expr, :real_cost), Tuple{Float64, Float64}}}(
        undef,
        n_zones,
    )

    init_dc_hybrid_partition!(
        partitions,
        sys,
        optimizer_factory,
        rho_theta[],
        rho_flow[],
        interzonal_branches,
    )

    iter_done = 0
    for iter in 1:max_iter
        iter_done = iter
        _FAST_USE_LOGS && println("\n[ADMM-hybrid/DC] Iteração $iter")

        refresh_rho_bus!(partitions, rho_theta, prev_rho_theta)
        refresh_rho_branch!(partitions, rho_flow, prev_rho_flow)

        Threads.@threads for i in 1:n_zones
            zone_id = zone_ids[i]
            part = get_partition(partitions, zone_id)
            model = part.model::JuMP.Model

            theta_dual_buf_z = dual_theta_buf[zone_id]
            theta_dual_buf_z .= duals_theta[zone_id]

            lin_bus = _build_linear_bus_expr(
                part,
                shared_buses_sorted,
                shared_buses_to_pos,
                theta_consensus,
                theta_dual_buf_z,
                rho_theta[],
            )
            lin_branch = _build_linear_branch_expr(
                part,
                part.zone_edges,
                flow_consensus,
                duals_flow[zone_id],
                rho_flow[],
            )

            @objective(
                model,
                Min,
                part.obj_base + part.quad_bus_expr + part.quad_branch_expr + lin_bus +
                lin_branch
            )
            optimize!(model)

            results[i] = (
                gen_cost_expr = value(model[:generation_cost]),
                real_cost = compute_real_generation_cost(model, part.objective_params),
            )
        end

        for (i, zone_id) in enumerate(zone_ids)
            part = get_partition(partitions, zone_id)
            out = results[i]
            part.has_solution = true
            part.generation_cost_model_expr = out.gen_cost_expr
            part.real_generation_cost = out.real_cost
        end

        theta_consensus_prev .= theta_consensus
        flow_consensus_prev = copy(flow_consensus)

        compute_theta_consensus_vec!(
            theta_consensus,
            partitions,
            shared_buses_sorted,
            shared_buses_to_pos,
        )
        @inbounds for b in shared_buses_sorted
            j = shared_buses_to_pos[b]
            push!(get!(theta_avg_history, b, Float64[]), theta_consensus[j])
        end

        compute_flow_consensus!(flow_consensus, partitions, interzonal_branches)

        update_duals_theta_all!(
            duals_theta,
            partitions,
            zone_ids,
            theta_consensus,
            rho_theta,
            shared_buses_sorted,
            shared_buses_to_pos,
        )
        update_duals_flow!(
            duals_flow,
            partitions,
            flow_consensus,
            interzonal_branches,
            rho_flow[],
        )

        r_bus = primal_residual_bus_dc(
            partitions,
            theta_consensus,
            shared_buses_sorted,
            shared_buses_to_pos,
        )
        r_br = primal_residual_branch_dc(partitions, flow_consensus, interzonal_branches)
        r_tot = hypot(r_bus, r_br)

        s_bus = dual_residual_bus_dc(theta_consensus, theta_consensus_prev)
        s_br = dual_residual_branch_dc(flow_consensus, flow_consensus_prev)
        s_tot = hypot(s_bus, s_br)

        eps_bus_r = eps_primal_bus_dc(
            partitions,
            theta_consensus,
            shared_buses_sorted,
            shared_buses_to_pos,
            eps_abs,
            eps_rel,
        )
        eps_br_r = eps_primal_branch_dc(
            partitions,
            flow_consensus,
            interzonal_branches,
            eps_abs,
            eps_rel,
        )
        eps_r_tot = hypot(eps_bus_r, eps_br_r)

        eps_bus_s = eps_dual_bus_dc(theta_consensus, theta_consensus_prev, eps_abs, eps_rel)
        eps_br_s = eps_dual_branch_dc(flow_consensus, flow_consensus_prev, eps_abs, eps_rel)
        eps_s_tot = hypot(eps_bus_s, eps_br_s)

        _record_iter!(
            hist, r_tot, s_tot, eps_r_tot, eps_s_tot;
            rhos = Dict(:rho_theta => rho_theta[], :rho_flow => rho_flow[]),
        )

        _FAST_USE_LOGS && println("‖r‖ = $r_tot | eps_pri = $eps_r_tot")
        _FAST_USE_LOGS && println("‖s‖ = $s_tot | eps_du  = $eps_s_tot")

        if stop_if_converged!(r_tot, s_tot, eps_r_tot, eps_s_tot)
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        maybe_adapt_rho!(
            adaptive,
            rho_theta, r_bus, s_bus,
            rho_flow, r_br, s_br,
            config,
        )
    end

    hist.heatmaps[:theta_avg] = theta_avg_history

    result = _finalize_admm_result(partitions, t0_ns, iter_done)
    _print_admm_summary!(
        "[ADMM-hybrid/DC]",
        result;
        rho_lines = (
            "rho_theta inicial: $(config.rho_theta) | rho_theta final: $(rho_theta[])",
            "rho_flow  inicial: $(config.rho_flow_p) | rho_flow  final: $(rho_flow[])",
        ),
    )

    return LirionRunResult(;
        result = result,
        model = nothing,
        system = sys,
        diagnostics = hist,
    )
end
