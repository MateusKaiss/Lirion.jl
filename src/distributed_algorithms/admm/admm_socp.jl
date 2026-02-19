# ============================================================
# ADMM-SOCP
# ============================================================

# ------------------------------------------------------------
# Partition initialization
# ------------------------------------------------------------

# ---------- socp state [V² + (R,T)] partition ----------

function init_socp_v2rt_partition!(
    part::PartitionSOCP,
    sys_zone::SystemData,
    optimizer_factory,
    rho_v2::Float64,
    rho_r::Float64,
    rho_t::Float64,
    interzonal_branches,
)
    prepare_hotpaths!(part)
    part.zone_edges = _zone_cut_edges(part, interzonal_branches)

    model = build_model(sys_zone, SOCP; optimizer = optimizer_factory)
    part.model = model
    part.obj_base = model[:generation_cost]
    part.v2_var = model[:v_magnitude_sq]
    part.objective_params = sys_zone.objective
    part.has_solution = false

    for I in keys(part.v2_var)
        set_start_value(part.v2_var[I], 1.0)
    end

    build_quad_bus_expr!(part, rho_v2)
    build_quad_branch_rt_expr!(part, rho_r, rho_t, part.zone_edges)

    part.obj_static = part.obj_base + part.quad_bus_expr_socp + part.quad_rt_expr_socp
    @objective(model, Min, part.obj_static)

    return nothing
end

function init_socp_v2rt_partition!(
    partitions::Partitions{PartitionSOCP},
    sys::DistributedSystem,
    optimizer_factory,
    rho_v2::Float64,
    rho_r::Float64,
    rho_t::Float64,
    interzonal_branches,
)
    n_zones = length(partitions.parts)
    @assert n_zones == length(sys.parts) == length(sys.zones)

    @inbounds for i in 1:n_zones
        init_socp_v2rt_partition!(
            partitions.parts[i],
            sys.parts[i],
            optimizer_factory,
            rho_v2,
            rho_r,
            rho_t,
            interzonal_branches,
        )
    end
    return nothing
end

# ---------- socp power partition (P/Q) ----------

function init_socp_pq_partition!(
    part::PartitionSOCP,
    sys_zone::SystemData,
    optimizer_factory,
    rho_p::Float64,
    rho_q::Float64,
    interzonal_branches,
)
    prepare_hotpaths!(part)
    part.zone_edges = _zone_cut_edges(part, interzonal_branches)

    model = build_model(sys_zone, SOCP; optimizer = optimizer_factory)
    part.model = model
    part.obj_base = model[:generation_cost]
    part.objective_params = sys_zone.objective
    part.has_solution = false

    if is_gurobi(model)
        soft_anchor_branch_expr!(
            part;
            κV2 = 1e-6,
            κPQ = 1e-8,
            use_bus = true,
            use_pq = true,
        )
    else
        part.soft_anchor_expr = zero(JuMP.QuadExpr)
    end

    build_quad_branch_pq_expr!(part, rho_p, rho_q, part.zone_edges)

    part.obj_static = part.obj_base + part.quad_pq_expr_socp + part.soft_anchor_expr
    @objective(model, Min, part.obj_static)

    return nothing
end

function init_socp_pq_partition!(
    partitions::Partitions{PartitionSOCP},
    sys::DistributedSystem,
    optimizer_factory,
    rho_p::Float64,
    rho_q::Float64,
    interzonal_branches,
)
    n_zones = length(partitions.parts)
    @assert n_zones == length(sys.parts) == length(sys.zones)

    @inbounds for i in 1:n_zones
        init_socp_pq_partition!(
            partitions.parts[i],
            sys.parts[i],
            optimizer_factory,
            rho_p,
            rho_q,
            interzonal_branches,
        )
    end
    return nothing
end

# ---------- socp hybrid [(V² + P/Q + R/T)] partition ----------

function init_socp_hybrid_partition!(
    part::PartitionSOCP,
    sys_zone::SystemData,
    optimizer_factory,
    rho_v2::Float64,
    rho_p::Float64,
    rho_q::Float64,
    rho_r::Float64,
    rho_t::Float64,
    interzonal_branches,
)
    prepare_hotpaths!(part)
    part.zone_edges = _zone_cut_edges(part, interzonal_branches)

    model = build_model(sys_zone, SOCP; optimizer = optimizer_factory)
    part.model = model
    part.obj_base = model[:generation_cost]
    part.v2_var = model[:v_magnitude_sq]
    part.objective_params = sys_zone.objective
    part.has_solution = false

    for I in keys(part.v2_var)
        set_start_value(part.v2_var[I], 1.0)
    end

    build_quad_bus_expr!(part, rho_v2)
    build_quad_branch_pq_expr!(part, rho_p, rho_q, part.zone_edges)
    build_quad_branch_rt_expr!(part, rho_r, rho_t, part.zone_edges)

    part.obj_static =
        part.obj_base + part.quad_bus_expr_socp + part.quad_pq_expr_socp +
        part.quad_rt_expr_socp
    @objective(model, Min, part.obj_static)

    return nothing
end

function init_socp_hybrid_partition!(
    partitions::Partitions{PartitionSOCP},
    sys::DistributedSystem,
    optimizer_factory,
    rho_v2::Float64,
    rho_p::Float64,
    rho_q::Float64,
    rho_r::Float64,
    rho_t::Float64,
    interzonal_branches,
)
    n_zones = length(partitions.parts)
    @assert n_zones == length(sys.parts) == length(sys.zones)

    @inbounds for i in 1:n_zones
        init_socp_hybrid_partition!(
            partitions.parts[i],
            sys.parts[i],
            optimizer_factory,
            rho_v2,
            rho_p,
            rho_q,
            rho_r,
            rho_t,
            interzonal_branches,
        )
    end
    return nothing
end

# ============================================================
# RUNNERS (ADMM-SOCP)
# ============================================================

# ------------------------------------------------------------
# Runner: ADMM SOCP state (V² + R/T)
# ------------------------------------------------------------

function run_admm(
    sys::DistributedSystem,
    config::DistributedOPFModel,
    ::ADMM{State, SOCPTag};
    optimizer_factory = nothing,
)::LirionRunResult
    t0_ns = time_ns()
    @assert config.base_model == SOCP "ADMM-SOCP: base_model deve ser SOCP."
    @assert config.coupling == :state "Dispatch State/SOCP requer coupling=:state."

    if optimizer_factory === nothing
        optimizer_factory = resolve_optimizer(config)
    end

    shared_buses = sys.shared_buses
    interzonal_branches = sys.interzonal_branches
    max_iter = config.max_iter

    rho_v2 = Ref(config.rho_v2)
    rho_r = Ref(config.rho_r)
    rho_t = Ref(config.rho_t)

    prev_rho_v2 = Ref(rho_v2[])
    prev_rho_r = Ref(rho_r[])
    prev_rho_t = Ref(rho_t[])

    eps_abs = config.eps_abs
    eps_rel = config.eps_rel
    adaptive = _adaptive_on(config)

    partitions = build_socp_partitions(sys)
    zone_ids = partitions.zones
    n_zones = length(zone_ids)

    shared_buses_sorted, shared_buses_to_pos = sorted_ids_and_id2idx(shared_buses)
    n_shared_buses = length(shared_buses_sorted)

    v2_consensus = ones(Float64, n_shared_buses)
    v2_consensus_prev = similar(v2_consensus)

    duals_v2 = Dict{Int, Vector{Float64}}(
        zone_id => zeros(Float64, n_shared_buses) for zone_id in zone_ids
    )
    dual_v2_buf = Dict{Int, Vector{Float64}}(
        zone_id => zeros(Float64, n_shared_buses) for zone_id in zone_ids
    )

    r_consensus = Dict{Int, Float64}(e => 0.0 for e in interzonal_branches)
    t_consensus = Dict{Int, Float64}(e => 0.0 for e in interzonal_branches)
    r_consensus_prev = Dict{Int, Float64}()
    t_consensus_prev = Dict{Int, Float64}()

    duals_r = Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in zone_ids)
    duals_t = Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in zone_ids)

    hist = ADMMHistory()
    hist.meta[:kind] = :socp_state_v2rt
    hist.meta[:tag] = "admm_socp_state_v2rt"
    v2_avg_history = Dict{Int, Vector{Float64}}()

    results = Vector{NamedTuple{(:gen_cost_expr, :real_cost), Tuple{Float64, Float64}}}(
        undef, n_zones,
    )

    init_socp_v2rt_partition!(
        partitions, sys, optimizer_factory,
        rho_v2[], rho_r[], rho_t[], interzonal_branches,
    )

    iter_done = 0
    for iter in 1:max_iter
        iter_done = iter
        _FAST_USE_LOGS && println("\n[ADMM-SOCP-state (V²+R,T)] Iteração $iter")

        refresh_rho_v2_socp!(partitions, rho_v2, prev_rho_v2)
        refresh_rho_rt_socp!(partitions, rho_r, prev_rho_r, rho_t, prev_rho_t)

        Threads.@threads for i in 1:n_zones
            zone_id = zone_ids[i]
            part = get_partition(partitions, zone_id)
            model = part.model::JuMP.Model

            v2_dual_buf_z = dual_v2_buf[zone_id]
            v2_dual_buf_z .= duals_v2[zone_id]

            lin_bus = build_linear_bus_term(
                part,
                shared_buses_sorted,
                shared_buses_to_pos,
                v2_consensus,
                v2_dual_buf_z,
                rho_v2[],
            )

            lin_rt = build_linear_branch_rt_term(
                part,
                part.zone_edges,
                r_consensus,
                t_consensus,
                duals_r[zone_id],
                duals_t[zone_id],
                rho_r[],
                rho_t[],
            )

            @objective(model, Min, part.obj_static + (lin_bus + lin_rt))
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

        v2_consensus_prev .= v2_consensus
        r_consensus_prev = copy(r_consensus)
        t_consensus_prev = copy(t_consensus)

        compute_v2_consensus_vec!(
            v2_consensus,
            partitions,
            shared_buses_sorted,
            shared_buses_to_pos,
        )
        compute_rt_consensus!(r_consensus, t_consensus, partitions, interzonal_branches)

        @inbounds for b in shared_buses_sorted
            j = shared_buses_to_pos[b]
            push!(get!(v2_avg_history, b, Float64[]), v2_consensus[j])
        end

        update_duals_v2_all!(
            duals_v2,
            partitions,
            zone_ids,
            v2_consensus,
            rho_v2[],
            shared_buses_sorted,
            shared_buses_to_pos,
        )
        update_duals_rt_all!(
            duals_r,
            duals_t,
            partitions,
            zone_ids,
            r_consensus,
            t_consensus,
            rho_r[],
            rho_t[],
        )

        r_bus = primal_residual_bus_socp_vec(
            partitions,
            v2_consensus,
            shared_buses_sorted,
            shared_buses_to_pos,
        )
        r_rt = primal_residual_branch_rt_socp(
            partitions,
            r_consensus,
            t_consensus,
            interzonal_branches,
        )
        r = hypot(r_bus, r_rt)

        s_bus = dual_residual_bus_socp_vec(v2_consensus, v2_consensus_prev, rho_v2[])
        s_rt = dual_residual_branch_rt_socp(
            r_consensus,
            t_consensus,
            r_consensus_prev,
            t_consensus_prev,
            interzonal_branches,
            rho_r[],
            rho_t[],
        )
        s = hypot(s_bus, s_rt)

        eps_r = begin
            epsb = eps_primal_bus_socp_vec(
                partitions,
                v2_consensus,
                shared_buses_sorted,
                shared_buses_to_pos,
                eps_abs,
                eps_rel,
            )
            epsrt = eps_primal_branch_rt_socp(
                partitions,
                r_consensus,
                t_consensus,
                interzonal_branches,
                eps_abs,
                eps_rel,
            )
            hypot(epsb, epsrt)
        end

        eps_s = begin
            epsb = eps_dual_bus_socp_unscaled(
                duals_v2,
                partitions,
                zone_ids,
                shared_buses_to_pos,
                eps_abs,
                eps_rel,
            )
            epsrt = eps_dual_branch_rt_socp_unscaled(
                duals_r,
                duals_t,
                zone_ids,
                interzonal_branches,
                eps_abs,
                eps_rel,
            )
            hypot(epsb, epsrt)
        end

        _record_iter!(
            hist,
            r,
            s,
            eps_r,
            eps_s;
            rhos = Dict(:rho_v2 => rho_v2[], :rho_r => rho_r[], :rho_t => rho_t[]),
        )

        _FAST_USE_LOGS && println("‖r‖ = $r | eps_pri = $eps_r")
        _FAST_USE_LOGS && println("‖s‖ = $s | eps_du  = $eps_s")

        if stop_if_converged!(r, s, eps_r, eps_s)
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        maybe_adapt_rho!(adaptive, rho_v2, r_bus, s_bus, rho_r, rho_t, r_rt, s_rt, config)
    end

    hist.heatmaps[:v2_avg] = v2_avg_history

    result = _finalize_admm_result(partitions, t0_ns, iter_done)
    _print_admm_summary!(
        "[ADMM-SOCP-state (V²+R,T)]",
        result;
        rho_lines = (
            "rho(V²) inicial: $(config.rho_v2) | rho(V²) final: $(rho_v2[])",
            "rho_r   inicial: $(config.rho_r)  | rho_r final: $(rho_r[])",
            "rho_t   inicial: $(config.rho_t)  | rho_t final: $(rho_t[])",
        ),
    )

    return LirionRunResult(; result, model = nothing, system = sys, diagnostics = hist)
end

# ------------------------------------------------------------
# Runner: ADMM SOCP power (branch P/Q)
# ------------------------------------------------------------

function run_admm(
    sys::DistributedSystem,
    config::DistributedOPFModel,
    ::ADMM{Power, SOCPTag};
    optimizer_factory = nothing,
)::LirionRunResult
    t0_ns = time_ns()
    @assert config.base_model == SOCP "ADMM-SOCP-power: base_model deve ser SOCP."
    @assert config.coupling == :power "Dispatch Power/SOCP requer coupling=:power."

    if optimizer_factory === nothing
        optimizer_factory = resolve_optimizer(config)
    end

    interzonal_branches = sys.interzonal_branches
    max_iter = config.max_iter

    rho_p = Ref(config.rho_p)
    rho_q = Ref(config.rho_q)
    prev_rho_p = Ref(rho_p[])
    prev_rho_q = Ref(rho_q[])

    eps_abs = config.eps_abs
    eps_rel = config.eps_rel
    adaptive = _adaptive_on(config)

    partitions = build_socp_partitions(sys)
    zone_ids = partitions.zones
    n_zones = length(zone_ids)

    p_flow_consensus = Dict{Int, Float64}(e => 0.0 for e in interzonal_branches)
    q_flow_consensus = Dict{Int, Float64}(e => 0.0 for e in interzonal_branches)
    p_flow_consensus_prev = Dict{Int, Float64}()
    q_flow_consensus_prev = Dict{Int, Float64}()

    dualsP = Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in zone_ids)
    dualsQ = Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in zone_ids)

    hist = ADMMHistory()
    hist.meta[:kind] = :socp_power_pq
    hist.meta[:tag] = "admm_socp_power_pq"

    init_socp_pq_partition!(
        partitions,
        sys,
        optimizer_factory,
        rho_p[],
        rho_q[],
        interzonal_branches,
    )

    iter_done = 0
    for iter in 1:max_iter
        iter_done = iter
        _FAST_USE_LOGS && println("\n[ADMM-SOCP-power (P,Q)] Iteração $iter")

        refresh_rho_pq_socp!(partitions, rho_p, prev_rho_p, rho_q, prev_rho_q)

        Threads.@threads for i in 1:n_zones
            zone_id = zone_ids[i]
            part = get_partition(partitions, zone_id)
            model = part.model::JuMP.Model

            lin = build_linear_branch_pq_term(
                part,
                part.zone_edges,
                p_flow_consensus,
                q_flow_consensus,
                dualsP[zone_id],
                dualsQ[zone_id],
                rho_p[],
                rho_q[],
            )

            @objective(model, Min, part.obj_static + lin)
            optimize!(model)

            part.has_solution = true
            part.generation_cost_model_expr = value(model[:generation_cost])
            part.real_generation_cost =
                compute_real_generation_cost(model, part.objective_params)
        end

        p_flow_consensus_prev = copy(p_flow_consensus)
        q_flow_consensus_prev = copy(q_flow_consensus)

        compute_pq_consensus!(
            p_flow_consensus,
            q_flow_consensus,
            partitions,
            interzonal_branches,
        )

        update_duals_pq_all!(
            dualsP,
            dualsQ,
            partitions,
            zone_ids,
            p_flow_consensus,
            q_flow_consensus,
            rho_p[],
            rho_q[],
        )

        r = primal_residual_branch_pq_socp(
            partitions,
            p_flow_consensus,
            q_flow_consensus,
            interzonal_branches,
        )
        s = dual_residual_branch_pq_socp(
            p_flow_consensus,
            q_flow_consensus,
            p_flow_consensus_prev,
            q_flow_consensus_prev,
            interzonal_branches,
            rho_p[],
            rho_q[],
        )

        eps_r = eps_primal_branch_pq_socp(
            partitions,
            p_flow_consensus,
            q_flow_consensus,
            interzonal_branches,
            eps_abs,
            eps_rel,
        )
        eps_s = eps_dual_branch_pq_socp_unscaled(
            dualsP,
            dualsQ,
            zone_ids,
            interzonal_branches,
            eps_abs,
            eps_rel,
        )

        _record_iter!(
            hist,
            r,
            s,
            eps_r,
            eps_s;
            rhos = Dict(:rho_p => rho_p[], :rho_q => rho_q[]),
        )

        _FAST_USE_LOGS && println("‖r‖ = $r | eps_pri = $eps_r")
        _FAST_USE_LOGS && println("‖s‖ = $s | eps_du  = $eps_s")

        if stop_if_converged!(r, s, eps_r, eps_s)
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        maybe_adapt_rho!(adaptive, rho_p, rho_q, r, s, config)
    end

    result = _finalize_admm_result(partitions, t0_ns, iter_done)
    _print_admm_summary!(
        "[ADMM-SOCP-power (P,Q)]",
        result;
        rho_lines = (
            "rho_p inicial: $(config.rho_p) | rho_p final: $(rho_p[])",
            "rho_q inicial: $(config.rho_q) | rho_q final: $(rho_q[])",
        ),
    )

    return LirionRunResult(; result, model = nothing, system = sys, diagnostics = hist)
end

# ------------------------------------------------------------
# Runner: ADMM SOCP hybrid (V² + P/Q + R/T)
# ------------------------------------------------------------

function run_admm(
    sys::DistributedSystem,
    config::DistributedOPFModel,
    ::ADMM{Hybrid, SOCPTag};
    optimizer_factory = nothing,
)::LirionRunResult
    t0_ns = time_ns()
    @assert config.base_model == SOCP "ADMM-SOCP-hybrid: base_model deve ser SOCP."
    @assert config.coupling == :hybrid "Dispatch Hybrid/SOCP requer coupling=:hybrid."

    if optimizer_factory === nothing
        optimizer_factory = resolve_optimizer(config)
    end

    shared_buses = sys.shared_buses
    interzonal_branches = sys.interzonal_branches
    max_iter = config.max_iter

    rho_v2 = Ref(config.rho_v2)
    rho_p = Ref(config.rho_p)
    rho_q = Ref(config.rho_q)
    rho_r = Ref(config.rho_r)
    rho_t = Ref(config.rho_t)

    prev_rho_v2 = Ref(rho_v2[])
    prev_rho_p = Ref(rho_p[])
    prev_rho_q = Ref(rho_q[])
    prev_rho_r = Ref(rho_r[])
    prev_rho_t = Ref(rho_t[])

    eps_abs = config.eps_abs
    eps_rel = config.eps_rel
    adaptive = _adaptive_on(config)

    partitions = build_socp_partitions(sys)
    zone_ids = partitions.zones
    n_zones = length(zone_ids)

    shared_buses_sorted, shared_buses_to_pos = sorted_ids_and_id2idx(shared_buses)
    n_shared_buses = length(shared_buses_sorted)

    v2_consensus = ones(Float64, n_shared_buses)
    v2_consensus_prev = similar(v2_consensus)

    duals_v2 =
        Dict{Int, Vector{Float64}}(z => zeros(Float64, n_shared_buses) for z in zone_ids)
    dual_v2_buf =
        Dict{Int, Vector{Float64}}(z => zeros(Float64, n_shared_buses) for z in zone_ids)

    p_flow_consensus = Dict{Int, Float64}(e => 0.0 for e in interzonal_branches)
    q_flow_consensus = Dict{Int, Float64}(e => 0.0 for e in interzonal_branches)
    r_consensus = Dict{Int, Float64}(e => 0.0 for e in interzonal_branches)
    t_consensus = Dict{Int, Float64}(e => 0.0 for e in interzonal_branches)

    p_flow_consensus_prev = Dict{Int, Float64}()
    q_flow_consensus_prev = Dict{Int, Float64}()
    r_consensus_prev = Dict{Int, Float64}()
    t_consensus_prev = Dict{Int, Float64}()

    dualsP = Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in zone_ids)
    dualsQ = Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in zone_ids)
    duals_r = Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in zone_ids)
    duals_t = Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in zone_ids)

    hist = ADMMHistory()
    hist.meta[:kind] = :socp_hybrid
    hist.meta[:tag] = "admm_socp_hybrid"
    v2_avg_history = Dict{Int, Vector{Float64}}()

    init_socp_hybrid_partition!(
        partitions, sys, optimizer_factory,
        rho_v2[], rho_p[], rho_q[], rho_r[], rho_t[], interzonal_branches,
    )

    iter_done = 0
    for iter in 1:max_iter
        iter_done = iter
        _FAST_USE_LOGS && println("\n[ADMM-SOCP-hybrid] Iteração $iter")

        refresh_rho_hybrid_socp!(
            partitions,
            rho_v2, prev_rho_v2,
            rho_p, prev_rho_p,
            rho_q, prev_rho_q,
            rho_r, prev_rho_r,
            rho_t, prev_rho_t,
        )

        Threads.@threads for i in 1:n_zones
            zone_id = zone_ids[i]
            part = get_partition(partitions, zone_id)
            model = part.model::JuMP.Model

            v2_dual_buf_z = dual_v2_buf[zone_id]
            v2_dual_buf_z .= duals_v2[zone_id]

            lin_bus = build_linear_bus_term(
                part,
                shared_buses_sorted,
                shared_buses_to_pos,
                v2_consensus,
                v2_dual_buf_z,
                rho_v2[],
            )

            lin_pq = build_linear_branch_pq_term(
                part,
                part.zone_edges,
                p_flow_consensus,
                q_flow_consensus,
                dualsP[zone_id],
                dualsQ[zone_id],
                rho_p[],
                rho_q[],
            )

            lin_rt = build_linear_branch_rt_term(
                part,
                part.zone_edges,
                r_consensus,
                t_consensus,
                duals_r[zone_id],
                duals_t[zone_id],
                rho_r[],
                rho_t[],
            )

            @objective(model, Min, part.obj_static + (lin_bus + lin_pq + lin_rt))
            optimize!(model)

            part.has_solution = true
            part.generation_cost_model_expr = value(model[:generation_cost])
            part.real_generation_cost =
                compute_real_generation_cost(model, part.objective_params)
        end

        v2_consensus_prev .= v2_consensus
        p_flow_consensus_prev = copy(p_flow_consensus)
        q_flow_consensus_prev = copy(q_flow_consensus)
        r_consensus_prev = copy(r_consensus)
        t_consensus_prev = copy(t_consensus)

        compute_v2_consensus_vec!(
            v2_consensus,
            partitions,
            shared_buses_sorted,
            shared_buses_to_pos,
        )
        compute_pq_consensus!(
            p_flow_consensus,
            q_flow_consensus,
            partitions,
            interzonal_branches,
        )
        compute_rt_consensus!(r_consensus, t_consensus, partitions, interzonal_branches)

        @inbounds for b in shared_buses_sorted
            j = shared_buses_to_pos[b]
            push!(get!(v2_avg_history, b, Float64[]), v2_consensus[j])
        end

        update_duals_v2_all!(
            duals_v2,
            partitions,
            zone_ids,
            v2_consensus,
            rho_v2[],
            shared_buses_sorted,
            shared_buses_to_pos,
        )
        update_duals_pq_all!(
            dualsP,
            dualsQ,
            partitions,
            zone_ids,
            p_flow_consensus,
            q_flow_consensus,
            rho_p[],
            rho_q[],
        )
        update_duals_rt_all!(
            duals_r,
            duals_t,
            partitions,
            zone_ids,
            r_consensus,
            t_consensus,
            rho_r[],
            rho_t[],
        )

        r_bus = primal_residual_bus_socp_vec(
            partitions,
            v2_consensus,
            shared_buses_sorted,
            shared_buses_to_pos,
        )
        r_pq = primal_residual_branch_pq_socp(
            partitions,
            p_flow_consensus,
            q_flow_consensus,
            interzonal_branches,
        )
        r_rt = primal_residual_branch_rt_socp(
            partitions,
            r_consensus,
            t_consensus,
            interzonal_branches,
        )
        r = sqrt(r_bus^2 + r_pq^2 + r_rt^2)

        s_bus = dual_residual_bus_socp_vec(v2_consensus, v2_consensus_prev, rho_v2[])
        s_pq = dual_residual_branch_pq_socp(
            p_flow_consensus,
            q_flow_consensus,
            p_flow_consensus_prev,
            q_flow_consensus_prev,
            interzonal_branches,
            rho_p[],
            rho_q[],
        )
        s_rt = dual_residual_branch_rt_socp(
            r_consensus,
            t_consensus,
            r_consensus_prev,
            t_consensus_prev,
            interzonal_branches,
            rho_r[],
            rho_t[],
        )
        s = sqrt(s_bus^2 + s_pq^2 + s_rt^2)

        eps_r = begin
            epsb = eps_primal_bus_socp_vec(
                partitions,
                v2_consensus,
                shared_buses_sorted,
                shared_buses_to_pos,
                eps_abs,
                eps_rel,
            )
            epspq = eps_primal_branch_pq_socp(
                partitions,
                p_flow_consensus,
                q_flow_consensus,
                interzonal_branches,
                eps_abs,
                eps_rel,
            )
            epsrt = eps_primal_branch_rt_socp(
                partitions,
                r_consensus,
                t_consensus,
                interzonal_branches,
                eps_abs,
                eps_rel,
            )
            sqrt(epsb^2 + epspq^2 + epsrt^2)
        end

        eps_s = begin
            epsb = eps_dual_bus_socp_unscaled(
                duals_v2,
                partitions,
                zone_ids,
                shared_buses_to_pos,
                eps_abs,
                eps_rel,
            )
            epspq = eps_dual_branch_pq_socp_unscaled(
                dualsP,
                dualsQ,
                zone_ids,
                interzonal_branches,
                eps_abs,
                eps_rel,
            )
            epsrt = eps_dual_branch_rt_socp_unscaled(
                duals_r,
                duals_t,
                zone_ids,
                interzonal_branches,
                eps_abs,
                eps_rel,
            )
            sqrt(epsb^2 + epspq^2 + epsrt^2)
        end

        _record_iter!(
            hist, r, s, eps_r, eps_s;
            rhos = Dict(
                :rho_v2 => rho_v2[],
                :rho_p => rho_p[],
                :rho_q => rho_q[],
                :rho_r => rho_r[],
                :rho_t => rho_t[],
            ),
        )

        _FAST_USE_LOGS && println("‖r‖ = $r | eps_pri = $eps_r")
        _FAST_USE_LOGS && println("‖s‖ = $s | eps_du  = $eps_s")

        if stop_if_converged!(r, s, eps_r, eps_s)
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        maybe_adapt_rho!(
            adaptive,
            rho_v2, r_bus, s_bus,
            rho_p, rho_q, r_pq, s_pq,
            rho_r, rho_t, r_rt, s_rt,
            config,
        )
    end

    hist.heatmaps[:v2_avg] = v2_avg_history

    result = _finalize_admm_result(partitions, t0_ns, iter_done)
    _print_admm_summary!(
        "[ADMM-SOCP-hybrid]",
        result;
        rho_lines = (
            "rho(V²) inicial: $(config.rho_v2) | rho(V²) final: $(rho_v2[])",
            "rho_p   inicial: $(config.rho_p)  | rho_p final: $(rho_p[])",
            "rho_q   inicial: $(config.rho_q)  | rho_q final: $(rho_q[])",
            "rho_r   inicial: $(config.rho_r)  | rho_r final: $(rho_r[])",
            "rho_t   inicial: $(config.rho_t)  | rho_t final: $(rho_t[])",
        ),
    )

    return LirionRunResult(; result, model = nothing, system = sys, diagnostics = hist)
end
