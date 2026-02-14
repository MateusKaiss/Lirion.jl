# ============================================================
# ADMM-DC
# ============================================================

# ------------------------------------------------------------
# Partition initialization
# ------------------------------------------------------------

# ---------- dc state coupling partition (ex-bus) ----------

function init_dc_state_partition!(
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

    build_quad_bus_expr!(part, rho_theta)
    @objective(model, Min, part.obj_base + part.quad_bus_expr)

    return nothing
end

function init_dc_state_partition!(
    partitions::Partitions{PartitionDC},
    sys::ADMMSystem,
    optimizer_factory,
    rho_theta::Float64,
)
    n_zones = length(partitions.parts)
    @assert length(sys.parts) == n_zones "ADMMSystem.parts must be aligned with partitions.parts"

    @inbounds for i in 1:n_zones
        init_dc_state_partition!(
            partitions.parts[i],
            sys.parts[i],
            optimizer_factory,
            rho_theta,
        )
    end
    return nothing
end

# ---------- dc power coupling partition (ex-branch) ----------

function init_dc_power_partition!(
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

    build_quad_branch_expr!(part, rho, part.zone_edges)
    @objective(model, Min, part.obj_base + part.quad_branch_expr)

    return nothing
end

function init_dc_power_partition!(
    partitions::Partitions{PartitionDC},
    sys::ADMMSystem,
    optimizer_factory,
    rho::Float64,
    interzonal_branches::Vector{Int},
)
    n_zones = length(partitions.parts)
    @assert length(sys.parts) == n_zones "ADMMSystem.parts must be aligned with partitions.parts"

    @inbounds for i in 1:n_zones
        init_dc_power_partition!(
            partitions.parts[i],
            sys.parts[i],
            optimizer_factory,
            rho,
            interzonal_branches,
        )
    end
    return nothing
end

# ---------- dc hybrid coupling partition ----------

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

    build_quad_bus_expr!(part, rho_theta)
    build_quad_branch_expr!(part, rho, part.zone_edges)

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
    n_zones = length(partitions.parts)
    @assert length(sys.parts) == n_zones "ADMMSystem.parts must be aligned with partitions.parts"

    @inbounds for i in 1:n_zones
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
# Runner: ADMM state/DC
# ------------------------------------------------------------
function run_admm(
    sys::ADMMSystem,
    config::DistributedOPFModel,
    ::ADMM{State, DCTag};
    optimizer_factory = nothing,
)::LirionRunResult
    t0_ns = time_ns()
    @assert config.base_model == DC "ADMM state/DC requer base_model == DC."
    @assert config.coupling == :state "Dispatch State/DC requer coupling=:state."

    shared_buses = sys.shared_buses
    max_iter = config.max_iter

    rho = Ref(config.rho_theta)
    prev_rho = Ref(rho[])
    eps_abs = config.eps_abs
    eps_rel = config.eps_rel
    adaptive = _adaptive_on(config)

    partitions = build_dc_partitions(sys)
    zone_ids = partitions.zones
    n_zones = length(zone_ids)

    shared_buses_sorted, shared_buses_to_pos = sorted_ids_and_id2idx(shared_buses)
    n_shared_buses = length(shared_buses_sorted)

    dual_theta =
        Dict{Int, Vector{Float64}}(z => zeros(Float64, n_shared_buses) for z in zone_ids)
    dual_theta_buf =
        Dict{Int, Vector{Float64}}(z => zeros(Float64, n_shared_buses) for z in zone_ids)

    theta_consensus = zeros(Float64, n_shared_buses)
    theta_consensus_prev = similar(theta_consensus)

    hist = ADMMHistory()
    hist.meta[:kind] = :dc_state
    hist.meta[:tag] = "admm_dc_state"

    theta_avg_history = Dict{Int, Vector{Float64}}()

    results = Vector{NamedTuple{(:gen_cost_expr, :real_cost), Tuple{Float64, Float64}}}(
        undef,
        n_zones,
    )

    init_dc_state_partition!(partitions, sys, optimizer_factory, rho[])

    iter_done = 0
    for iter in 1:max_iter
        iter_done = iter
        _FAST_USE_LOGS && println("\nIteração $iter (state/DC)")

        refresh_rho_bus_dc!(partitions, rho, prev_rho)

        Threads.@threads for i in 1:n_zones
            zone_id = zone_ids[i]
            part = get_partition(partitions, zone_id)
            model = part.model::JuMP.Model

            theta_dual_buf_z = dual_theta_buf[zone_id]
            theta_dual_buf_z .= dual_theta[zone_id]

            lin = build_linear_bus_term(
                part,
                shared_buses_sorted,
                shared_buses_to_pos,
                theta_consensus,
                theta_dual_buf_z,
                rho[],
            )

            @objective(model, Min, part.obj_base + part.quad_bus_expr + lin)
            optimize!(model)

            results[i] = (
                gen_cost_expr = value(model[:generation_cost]),
                real_cost = compute_real_generation_cost(model, part.objective_params),
            )
        end

        @inbounds for (i, zone_id) in enumerate(zone_ids)
            part = get_partition(partitions, zone_id)
            out = results[i]
            part.has_solution = true
            part.generation_cost_model_expr = out.gen_cost_expr
            part.real_generation_cost = out.real_cost
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

        s = dual_residual_bus_dc(theta_consensus, theta_consensus_prev, rho[])

        eps_r = eps_primal_bus_dc(
            partitions,
            theta_consensus,
            shared_buses_sorted,
            shared_buses_to_pos,
            eps_abs,
            eps_rel,
        )

        eps_s = eps_dual_bus_dc_unscaled(
            dual_theta, partitions, zone_ids, shared_buses_to_pos, eps_abs, eps_rel,
        )

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

    result = _finalize_admm_result(partitions, t0_ns, iter_done)
    _print_admm_summary!(
        "[ADMM-state/DC]",
        result;
        rho_lines = ("rho_theta inicial: $(config.rho_theta) | rho_theta final: $(rho[])",),
    )

    return LirionRunResult(; result, model = nothing, system = sys, diagnostics = hist)
end

# ------------------------------------------------------------
# Runner: ADMM power/DC
# ------------------------------------------------------------
function run_admm(
    sys::ADMMSystem,
    config::DistributedOPFModel,
    ::ADMM{Power, DCTag};
    optimizer_factory = nothing,
)::LirionRunResult
    t0_ns = time_ns()
    @assert config.base_model == DC "ADMM-power/DC requer base_model == DC."
    @assert config.coupling == :power "Dispatch Power/DC requer coupling=:power."

    interzonal_branches = sys.interzonal_branches
    max_iter = config.max_iter

    rho = Ref(config.rho_p)
    prev_rho = Ref(rho[])
    eps_abs = config.eps_abs
    eps_rel = config.eps_rel
    adaptive = _adaptive_on(config)

    partitions = build_dc_partitions(sys)
    zone_ids = partitions.zones
    n_zones = length(zone_ids)

    flow_consensus = Dict{Int, Float64}(e => 0.0 for e in interzonal_branches)
    flow_consensus_prev = Dict{Int, Float64}()

    duals_flow = Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in zone_ids)

    hist = ADMMHistory()
    hist.meta[:kind] = :dc_power
    hist.meta[:tag] = "admm_dc_power"

    results = Vector{NamedTuple{(:gen_cost_expr, :real_cost), Tuple{Float64, Float64}}}(
        undef,
        n_zones,
    )

    init_dc_power_partition!(
        partitions,
        sys,
        optimizer_factory,
        rho[],
        interzonal_branches,
    )

    iter_done = 0
    for iter in 1:max_iter
        iter_done = iter
        _FAST_USE_LOGS && println("\n[ADMM-power/DC] Iteração $iter")

        refresh_rho_branch_dc!(partitions, rho, prev_rho)

        Threads.@threads for i in 1:n_zones
            zone_id = zone_ids[i]
            part = get_partition(partitions, zone_id)
            model = part.model::JuMP.Model

            lin = build_linear_branch_term(
                part,
                part.zone_edges,
                flow_consensus,
                duals_flow[zone_id],
                rho[],
            )

            @objective(model, Min, part.obj_base + part.quad_branch_expr + lin)
            optimize!(model)
            _require_optimal!(model; iter = iter, zone_id = zone_id, tag = "ADMM-power/DC")

            results[i] = (
                gen_cost_expr = value(model[:generation_cost]),
                real_cost = compute_real_generation_cost(model, part.objective_params),
            )
        end

        @inbounds for (i, zone_id) in enumerate(zone_ids)
            part = get_partition(partitions, zone_id)
            out = results[i]
            part.has_solution = true
            part.generation_cost_model_expr = out.gen_cost_expr
            part.real_generation_cost = out.real_cost
        end

        flow_consensus_prev = copy(flow_consensus)

        compute_flow_consensus!(flow_consensus, partitions, interzonal_branches)

        update_duals_flow!(
            duals_flow,
            partitions,
            flow_consensus,
            interzonal_branches,
            rho[],
        )

        r = primal_residual_branch_dc(partitions, flow_consensus, interzonal_branches)
        s = dual_residual_branch_dc(flow_consensus, flow_consensus_prev, rho[])

        eps_r = eps_primal_branch_dc(
            partitions,
            flow_consensus,
            interzonal_branches,
            eps_abs,
            eps_rel,
        )

        eps_s = eps_dual_branch_dc_unscaled(
            duals_flow,
            partitions,
            zone_ids,
            interzonal_branches,
            eps_abs,
            eps_rel,
        )

        _record_iter!(hist, r, s, eps_r, eps_s; rhos = Dict(:rho_p => rho[]))

        _FAST_USE_LOGS && println("‖r‖ = $r | eps_pri = $eps_r")
        _FAST_USE_LOGS && println("‖s‖ = $s | eps_du  = $eps_s")

        if stop_if_converged!(r, s, eps_r, eps_s)
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        maybe_adapt_rho!(adaptive, rho, r, s, config)
    end

    result = _finalize_admm_result(partitions, t0_ns, iter_done)
    _print_admm_summary!(
        "[ADMM-power/DC]",
        result;
        rho_lines = ("rho_p inicial: $(config.rho_p) | rho_p final: $(rho[])",),
    )

    return LirionRunResult(; result, model = nothing, system = sys, diagnostics = hist)
end

# ------------------------------------------------------------
# Runner: ADMM hybrid/DC
# ------------------------------------------------------------
function run_admm(
    sys::ADMMSystem,
    config::DistributedOPFModel,
    ::ADMM{Hybrid, DCTag};
    optimizer_factory = nothing,
)::LirionRunResult
    t0_ns = time_ns()
    @assert config.base_model == DC "Hybrid ADMM/DC requer DC."
    @assert config.coupling == :hybrid "Dispatch Hybrid/DC requer coupling=:hybrid."

    shared_buses = sys.shared_buses
    interzonal_branches = sys.interzonal_branches
    max_iter = config.max_iter

    rho_theta = Ref(config.rho_theta)
    rho_p = Ref(config.rho_p)

    prev_rho_theta = Ref(rho_theta[])
    prev_rho_p = Ref(rho_p[])

    eps_abs = config.eps_abs
    eps_rel = config.eps_rel
    adaptive = _adaptive_on(config)

    partitions = build_dc_partitions(sys)
    zone_ids = partitions.zones
    n_zones = length(zone_ids)

    shared_buses_sorted, shared_buses_to_pos = sorted_ids_and_id2idx(shared_buses)
    n_shared_buses = length(shared_buses_sorted)

    theta_consensus = zeros(Float64, n_shared_buses)
    theta_consensus_prev = similar(theta_consensus)

    flow_consensus = Dict{Int, Float64}(e => 0.0 for e in interzonal_branches)
    flow_consensus_prev = Dict{Int, Float64}()

    duals_theta =
        Dict{Int, Vector{Float64}}(z => zeros(Float64, n_shared_buses) for z in zone_ids)
    dual_theta_buf =
        Dict{Int, Vector{Float64}}(z => zeros(Float64, n_shared_buses) for z in zone_ids)

    duals_flow = Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in zone_ids)

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
        rho_p[],
        interzonal_branches,
    )

    iter_done = 0
    for iter in 1:max_iter
        iter_done = iter
        _FAST_USE_LOGS && println("\n[ADMM-hybrid/DC] Iteração $iter")

        refresh_rho_bus_dc!(partitions, rho_theta, prev_rho_theta)
        refresh_rho_branch_dc!(partitions, rho_p, prev_rho_p)

        Threads.@threads for i in 1:n_zones
            zone_id = zone_ids[i]
            part = get_partition(partitions, zone_id)
            model = part.model::JuMP.Model

            theta_dual_buf_z = dual_theta_buf[zone_id]
            theta_dual_buf_z .= duals_theta[zone_id]

            lin_bus = build_linear_bus_term(
                part,
                shared_buses_sorted,
                shared_buses_to_pos,
                theta_consensus,
                theta_dual_buf_z,
                rho_theta[],
            )
            lin_branch = build_linear_branch_term(
                part,
                part.zone_edges,
                flow_consensus,
                duals_flow[zone_id],
                rho_p[],
            )

            @objective(
                model,
                Min,
                part.obj_base + part.quad_bus_expr + part.quad_branch_expr + lin_bus +
                lin_branch
            )
            optimize!(model)
            _require_optimal!(model; iter = iter, zone_id = zone_id, tag = "ADMM-hybrid/DC")

            results[i] = (
                gen_cost_expr = value(model[:generation_cost]),
                real_cost = compute_real_generation_cost(model, part.objective_params),
            )
        end

        @inbounds for (i, zone_id) in enumerate(zone_ids)
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
            rho_p[],
        )

        r_bus = primal_residual_bus_dc(
            partitions,
            theta_consensus,
            shared_buses_sorted,
            shared_buses_to_pos,
        )
        r_br = primal_residual_branch_dc(partitions, flow_consensus, interzonal_branches)
        r_tot = hypot(r_bus, r_br)

        s_bus = dual_residual_bus_dc(theta_consensus, theta_consensus_prev, rho_theta[])
        s_br = dual_residual_branch_dc(flow_consensus, flow_consensus_prev, rho_p[])
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

        eps_bus_s = eps_dual_bus_dc_unscaled(
            duals_theta, partitions, zone_ids, shared_buses_to_pos, eps_abs, eps_rel,
        )
        eps_br_s = eps_dual_branch_dc_unscaled(
            duals_flow, partitions, zone_ids, interzonal_branches, eps_abs, eps_rel,
        )
        eps_s_tot = hypot(eps_bus_s, eps_br_s)

        _record_iter!(
            hist, r_tot, s_tot, eps_r_tot, eps_s_tot;
            rhos = Dict(:rho_theta => rho_theta[], :rho_p => rho_p[]),
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
            rho_p, r_br, s_br,
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
            "rho_p     inicial: $(config.rho_p)     | rho_p     final: $(rho_p[])",
        ),
    )

    return LirionRunResult(; result, model = nothing, system = sys, diagnostics = hist)
end
