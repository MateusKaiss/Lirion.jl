# ============================================================
# ADMM-AC
# ============================================================

# ------------------------------------------------------------
# Partition initialization
# ------------------------------------------------------------

# ---------- ac state coupling partition ----------

function init_ac_state_partition!(
    part::PartitionAC,
    sys_zone::SystemData,
    optimizer_factory,
    rho_theta::Float64,
    rho_v_mag::Float64,
)
    prepare_hotpaths!(part)

    model = build_model(sys_zone, AC; optimizer = optimizer_factory)
    part.model = model
    part.obj_base = model[:generation_cost]
    part.theta_var = model[:theta]
    part.v_mag_var = model[:v_magnitude]
    part.objective_params = sys_zone.objective
    part.has_solution = false

    for I in keys(part.theta_var)
        set_start_value(part.theta_var[I], 0.0)
    end
    for I in keys(part.v_mag_var)
        set_start_value(part.v_mag_var[I], 1.0)
    end

    build_quad_bus_expr!(part, rho_theta, rho_v_mag)
    @objective(model, Min, part.obj_base + part.quad_bus_expr_ac)

    return nothing
end

function init_ac_state_partition!(
    partitions::Partitions{PartitionAC},
    sys::DistributedSystem,
    optimizer_factory,
    rho_theta::Float64,
    rho_v_mag::Float64,
)
    n_zones = length(partitions.parts)
    @assert n_zones == length(sys.parts) == length(sys.zones)

    @inbounds for i in 1:n_zones
        init_ac_state_partition!(
            partitions.parts[i],
            sys.parts[i],
            optimizer_factory,
            rho_theta,
            rho_v_mag,
        )
    end
    return nothing
end

# ---------- ac power coupling partition  ----------

function init_ac_power_partition!(
    part::PartitionAC,
    sys_zone::SystemData,
    optimizer_factory,
    rhoP::Float64,
    rhoQ::Float64,
    interzonal_branches,
)
    prepare_hotpaths!(part)
    part.zone_edges = _zone_cut_edges(part, interzonal_branches)

    model = build_model(sys_zone, AC; optimizer = optimizer_factory)
    part.model = model
    part.obj_base = model[:generation_cost]
    part.theta_var = model[:theta]
    part.v_mag_var = model[:v_magnitude]
    part.objective_params = sys_zone.objective
    part.has_solution = false

    for I in keys(part.theta_var)
        set_start_value(part.theta_var[I], 0.0)
    end
    for I in keys(part.v_mag_var)
        set_start_value(part.v_mag_var[I], 1.0)
    end

    build_quad_branch_expr!(part, rhoP, rhoQ, part.zone_edges)
    @objective(model, Min, part.obj_base + part.quad_branch_expr_ac)

    return nothing
end

function init_ac_power_partition!(
    partitions::Partitions{PartitionAC},
    sys::DistributedSystem,
    optimizer_factory,
    rhoP::Float64,
    rhoQ::Float64,
    interzonal_branches,
)
    n_zones = length(partitions.parts)
    @assert n_zones == length(sys.parts) == length(sys.zones)

    @inbounds for i in 1:n_zones
        init_ac_power_partition!(
            partitions.parts[i],
            sys.parts[i],
            optimizer_factory,
            rhoP,
            rhoQ,
            interzonal_branches,
        )
    end
    return nothing
end

# ---------- ac hybrid coupling partition ----------

function init_ac_hybrid_partition!(
    part::PartitionAC,
    sys_zone::SystemData,
    optimizer_factory,
    rho_theta::Float64,
    rho_v_mag::Float64,
    rhoP::Float64,
    rhoQ::Float64,
    interzonal_branches,
)
    prepare_hotpaths!(part)
    part.zone_edges = _zone_cut_edges(part, interzonal_branches)

    model = build_model(sys_zone, AC; optimizer = optimizer_factory)
    part.model = model
    part.obj_base = model[:generation_cost]
    part.theta_var = model[:theta]
    part.v_mag_var = model[:v_magnitude]
    part.objective_params = sys_zone.objective
    part.has_solution = false

    for I in keys(part.theta_var)
        set_start_value(part.theta_var[I], 0.0)
    end
    for I in keys(part.v_mag_var)
        set_start_value(part.v_mag_var[I], 1.0)
    end

    build_quad_bus_expr!(part, rho_theta, rho_v_mag)
    build_quad_branch_expr!(part, rhoP, rhoQ, part.zone_edges)

    @objective(
        model,
        Min,
        part.obj_base + part.quad_bus_expr_ac + part.quad_branch_expr_ac,
    )

    return nothing
end

function init_ac_hybrid_partition!(
    partitions::Partitions{PartitionAC},
    sys::DistributedSystem,
    optimizer_factory,
    rho_theta::Float64,
    rho_v_mag::Float64,
    rhoP::Float64,
    rhoQ::Float64,
    interzonal_branches,
)
    n_zones = length(partitions.parts)
    @assert n_zones == length(sys.parts) == length(sys.zones)

    @inbounds for i in 1:n_zones
        init_ac_hybrid_partition!(
            partitions.parts[i],
            sys.parts[i],
            optimizer_factory,
            rho_theta,
            rho_v_mag,
            rhoP,
            rhoQ,
            interzonal_branches,
        )
    end
    return nothing
end

# ============================================================
# RUNNERS (ADMM-AC)
# ============================================================

# ------------------------------------------------------------
# Runner: ADMM state/AC  (ex-bus)
# ------------------------------------------------------------
function run_admm(
    sys::DistributedSystem,
    config::DistributedOPFModel,
    ::ADMM{State, AC};
    optimizer_factory = nothing,
)::LirionRunResult
    t0_ns = time_ns()
    @assert config.base_model == AC "ADMM-AC: base_model deve ser AC."
    @assert config.coupling == :state "Dispatch State/AC requer coupling=:state."

    full_shared_bus_ids = sys.full_shared_bus_ids
    max_iter = config.max_iter

    rho_theta = Ref(config.rho_theta)
    rho_v = Ref(config.rho_v)
    prev_rho_theta = Ref(rho_theta[])
    prev_rho_v = Ref(rho_v[])
    eps_abs = config.eps_abs
    eps_rel = config.eps_rel
    adaptive = _adaptive_on(config)

    partitions = build_ac_partitions(sys)
    zone_ids = partitions.zones
    n_zones = length(zone_ids)

    shared_buses_sorted, shared_buses_to_pos = sorted_ids_and_id2idx(full_shared_bus_ids)
    n_shared_buses = length(shared_buses_sorted)

    duals_theta = Dict{Int, Vector{Float64}}(z => zeros(n_shared_buses) for z in zone_ids)
    duals_v_mag = Dict{Int, Vector{Float64}}(z => zeros(n_shared_buses) for z in zone_ids)

    dual_theta_buf =
        Dict{Int, Vector{Float64}}(z => zeros(n_shared_buses) for z in zone_ids)
    dual_v_mag_buf =
        Dict{Int, Vector{Float64}}(z => zeros(n_shared_buses) for z in zone_ids)

    theta_consensus = zeros(Float64, n_shared_buses)
    v_mag_consensus = ones(Float64, n_shared_buses)
    theta_consensus_prev = similar(theta_consensus)
    v_mag_consensus_prev = similar(v_mag_consensus)

    hist = ADMMHistory()
    hist.meta[:kind] = :ac_state
    hist.meta[:tag] = "admm_ac_state"

    theta_avg_hist = Dict{Int, Vector{Float64}}()
    v_mag_avg_hist = Dict{Int, Vector{Float64}}()

    results = Vector{NamedTuple{(:expr_cost, :real_cost), Tuple{Float64, Float64}}}(
        undef,
        n_zones,
    )

    init_ac_state_partition!(partitions, sys, optimizer_factory, rho_theta[], rho_v[])

    iter_done = 0
    for iter in 1:max_iter
        iter_done = iter
        _FAST_USE_LOGS && println("\n[ADMM-AC-state] Iteração $iter")

        refresh_rho_bus_ac!(partitions, rho_theta, prev_rho_theta, rho_v, prev_rho_v)

        Threads.@threads for i in 1:n_zones
            zone_id = zone_ids[i]
            part = partitions.parts[i]
            model = part.model::JuMP.Model

            theta_dual_buf_z = dual_theta_buf[zone_id]
            v_mag_dual_buf_z = dual_v_mag_buf[zone_id]
            theta_dual_buf_z .= duals_theta[zone_id]
            v_mag_dual_buf_z .= duals_v_mag[zone_id]

            lin = build_linear_bus_term(
                part,
                shared_buses_sorted,
                shared_buses_to_pos,
                theta_consensus,
                v_mag_consensus,
                theta_dual_buf_z,
                v_mag_dual_buf_z,
                rho_theta[],
                rho_v[],
            )

            @objective(model, Min, part.obj_base + part.quad_bus_expr_ac + lin)
            optimize!(model)

            results[i] = (
                expr_cost = value(model[:generation_cost]),
                real_cost = compute_real_generation_cost(model, part.objective_params),
            )
        end

        @inbounds for i in 1:n_zones
            part = partitions.parts[i]
            out = results[i]
            part.has_solution = true
            part.generation_cost_model_expr = out.expr_cost
            part.real_generation_cost = out.real_cost
        end

        theta_consensus_prev .= theta_consensus
        v_mag_consensus_prev .= v_mag_consensus

        compute_ac_bus_consensus_vec!(
            theta_consensus,
            v_mag_consensus,
            partitions,
            shared_buses_sorted,
            shared_buses_to_pos,
        )

        @inbounds for b in shared_buses_sorted
            j = shared_buses_to_pos[b]
            push!(get!(theta_avg_hist, b, Float64[]), theta_consensus[j])
            push!(get!(v_mag_avg_hist, b, Float64[]), v_mag_consensus[j])
        end

        update_duals_bus_ac_all!(
            duals_theta,
            duals_v_mag,
            partitions,
            zone_ids,
            theta_consensus,
            v_mag_consensus,
            rho_theta,
            rho_v,
            shared_buses_sorted,
            shared_buses_to_pos;
            theta_consensus_prev = theta_consensus_prev,
            v_mag_consensus_prev = v_mag_consensus_prev,
        )

        r = primal_residual_bus_ac_vec(
            partitions,
            theta_consensus,
            v_mag_consensus,
            shared_buses_sorted,
            shared_buses_to_pos,
        )

        s = dual_residual_bus_ac_vec(
            theta_consensus,
            v_mag_consensus,
            theta_consensus_prev,
            v_mag_consensus_prev,
            rho_theta[],
            rho_v[],
        )

        eps_r = eps_primal_bus_ac_vec(
            partitions,
            theta_consensus,
            v_mag_consensus,
            shared_buses_sorted,
            shared_buses_to_pos,
            eps_abs,
            eps_rel,
        )

        eps_s = eps_dual_bus_ac_unscaled(
            duals_theta,
            duals_v_mag,
            partitions,
            zone_ids,
            shared_buses_to_pos,
            eps_abs,
            eps_rel,
        )

        _record_iter!(
            hist,
            r,
            s,
            eps_r,
            eps_s;
            rhos = Dict(:rho_theta => rho_theta[], :rho_v => rho_v[]),
        )

        _FAST_USE_LOGS && println("‖r‖ = $r | eps_pri = $eps_r")
        _FAST_USE_LOGS && println("‖s‖ = $s | eps_du  = $eps_s")

        if stop_if_converged!(r, s, eps_r, eps_s)
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        maybe_adapt_rho!(
            adaptive,
            rho_theta,
            rho_v,
            partitions,
            theta_consensus,
            v_mag_consensus,
            theta_consensus_prev,
            v_mag_consensus_prev,
            shared_buses_sorted,
            shared_buses_to_pos,
            config,
        )
    end

    hist.heatmaps[:theta_avg] = theta_avg_hist
    hist.heatmaps[:v_magavg] = v_mag_avg_hist

    result = _finalize_admm_result(partitions, t0_ns, iter_done)
    _print_admm_summary!(
        "[ADMM-state/AC]",
        result;
        rho_lines = (
            "rho_theta inicial: $(config.rho_theta) | rho_theta final: $(rho_theta[])",
            "rho_v     inicial: $(config.rho_v)     | rho_v     final: $(rho_v[])",
        ),
    )

    return LirionRunResult(; result, model = nothing, system = sys, diagnostics = hist)
end

# ------------------------------------------------------------
# Runner: ADMM power/AC
# ------------------------------------------------------------
function run_admm(
    sys::DistributedSystem,
    config::DistributedOPFModel,
    ::ADMM{Power, AC};
    optimizer_factory = nothing,
)::LirionRunResult
    t0_ns = time_ns()
    @assert config.base_model == AC "ADMM-AC-power: base_model deve ser AC."
    @assert config.coupling == :power "Dispatch Power/AC requer coupling=:power."

    interzonal_branches = sys.interzonal_branches
    max_iter = config.max_iter

    rho_p = Ref(config.rho_p)
    rho_q = Ref(config.rho_q)
    prev_rho_p = Ref(rho_p[])
    prev_rho_q = Ref(rho_q[])

    eps_abs = config.eps_abs
    eps_rel = config.eps_rel
    adaptive = _adaptive_on(config)

    partitions = build_ac_partitions(sys)
    zone_ids = partitions.zones
    n_zones = length(zone_ids)

    duals_p_flow =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in zone_ids)
    duals_q_flow =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in zone_ids)

    p_flow_consensus = Dict{Int, Float64}(e => 0.0 for e in interzonal_branches)
    q_flow_consensus = Dict{Int, Float64}(e => 0.0 for e in interzonal_branches)
    p_flow_consensus_prev = Dict{Int, Float64}()
    q_flow_consensus_prev = Dict{Int, Float64}()

    hist = ADMMHistory()
    hist.meta[:kind] = :ac_power
    hist.meta[:tag] = "admm_ac_power"

    results = Vector{NamedTuple{(:expr_cost, :real_cost), Tuple{Float64, Float64}}}(
        undef,
        n_zones,
    )

    init_ac_power_partition!(
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
        _FAST_USE_LOGS && println("\n[ADMM-AC-power] Iteração $iter")

        refresh_rho_branch_ac!(partitions, rho_p, prev_rho_p, rho_q, prev_rho_q)

        Threads.@threads for i in 1:n_zones
            zone_id = zone_ids[i]
            part = partitions.parts[i]
            model = part.model::JuMP.Model

            lin = build_linear_branch_term(
                part,
                part.zone_edges,
                p_flow_consensus,
                q_flow_consensus,
                duals_p_flow[zone_id],
                duals_q_flow[zone_id],
                rho_p[],
                rho_q[],
            )

            @objective(model, Min, part.obj_base + part.quad_branch_expr_ac + lin)
            optimize!(model)

            results[i] = (
                expr_cost = value(model[:generation_cost]),
                real_cost = compute_real_generation_cost(model, part.objective_params),
            )
        end

        @inbounds for i in 1:n_zones
            part = partitions.parts[i]
            out = results[i]
            part.has_solution = true
            part.generation_cost_model_expr = out.expr_cost
            part.real_generation_cost = out.real_cost
        end

        p_flow_consensus_prev = copy(p_flow_consensus)
        q_flow_consensus_prev = copy(q_flow_consensus)

        compute_ac_branch_consensus!(
            p_flow_consensus,
            q_flow_consensus,
            partitions,
            interzonal_branches,
        )

        update_duals_branch_ac!(
            duals_p_flow,
            duals_q_flow,
            partitions,
            p_flow_consensus,
            q_flow_consensus,
            interzonal_branches,
            rho_p[],
            rho_q[];
            p_flow_consensus_prev = p_flow_consensus_prev,
            q_flow_consensus_prev = q_flow_consensus_prev,
        )

        r = primal_residual_branch_ac(
            partitions,
            p_flow_consensus,
            q_flow_consensus,
            interzonal_branches,
        )

        s = dual_residual_branch_ac(
            p_flow_consensus,
            q_flow_consensus,
            p_flow_consensus_prev,
            q_flow_consensus_prev,
            interzonal_branches,
            rho_p[],
            rho_q[],
        )

        eps_r = eps_primal_branch_ac(
            partitions,
            p_flow_consensus,
            q_flow_consensus,
            interzonal_branches,
            eps_abs,
            eps_rel,
        )

        eps_s = eps_dual_branch_ac_unscaled(
            duals_p_flow,
            duals_q_flow,
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

        maybe_adapt_rho_socp!(adaptive, rho_p, rho_q, r, s, config)
    end

    result = _finalize_admm_result(partitions, t0_ns, iter_done)
    _print_admm_summary!(
        "[ADMM-power/AC]",
        result;
        rho_lines = (
            "rhoP inicial: $(config.rho_p) | rhoP final: $(rho_p[])",
            "rhoQ inicial: $(config.rho_q) | rhoQ final: $(rho_q[])",
        ),
    )

    return LirionRunResult(; result, model = nothing, system = sys, diagnostics = hist)
end

# ------------------------------------------------------------
# Runner: ADMM hybrid/AC
# ------------------------------------------------------------
function run_admm(
    sys::DistributedSystem,
    config::DistributedOPFModel,
    ::ADMM{Hybrid, AC};
    optimizer_factory = nothing,
)::LirionRunResult
    t0_ns = time_ns()
    @assert config.base_model == AC "Hybrid AC: base_model deve ser AC."
    @assert config.coupling == :hybrid "Dispatch Hybrid/AC requer coupling=:hybrid."

    full_shared_bus_ids = sys.full_shared_bus_ids
    interzonal_branches = sys.interzonal_branches
    max_iter = config.max_iter

    rho_theta = Ref(config.rho_theta)
    rho_v = Ref(config.rho_v)
    rho_p = Ref(config.rho_p)
    rho_q = Ref(config.rho_q)

    prev_rho_theta = Ref(rho_theta[])
    prev_rho_v = Ref(rho_v[])
    prev_rho_p = Ref(rho_p[])
    prev_rho_q = Ref(rho_q[])

    eps_abs = config.eps_abs
    eps_rel = config.eps_rel
    adaptive = _adaptive_on(config)

    partitions = build_ac_partitions(sys)
    zone_ids = partitions.zones
    n_zones = length(zone_ids)

    shared_buses_sorted, shared_buses_to_pos = sorted_ids_and_id2idx(full_shared_bus_ids)
    n_shared_buses = length(shared_buses_sorted)

    duals_theta = Dict{Int, Vector{Float64}}(z => zeros(n_shared_buses) for z in zone_ids)
    duals_v_mag = Dict{Int, Vector{Float64}}(z => zeros(n_shared_buses) for z in zone_ids)

    dual_theta_buf =
        Dict{Int, Vector{Float64}}(z => zeros(n_shared_buses) for z in zone_ids)
    dual_v_mag_buf =
        Dict{Int, Vector{Float64}}(z => zeros(n_shared_buses) for z in zone_ids)

    duals_p_flow =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in zone_ids)
    duals_q_flow =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in zone_ids)

    theta_consensus = zeros(Float64, n_shared_buses)
    v_mag_consensus = ones(Float64, n_shared_buses)
    theta_consensus_prev = similar(theta_consensus)
    v_mag_consensus_prev = similar(v_mag_consensus)

    p_flow_consensus = Dict{Int, Float64}(e => 0.0 for e in interzonal_branches)
    q_flow_consensus = Dict{Int, Float64}(e => 0.0 for e in interzonal_branches)
    p_flow_consensus_prev = Dict{Int, Float64}()
    q_flow_consensus_prev = Dict{Int, Float64}()

    hist = ADMMHistory()
    hist.meta[:kind] = :ac_hybrid
    hist.meta[:tag] = "admm_ac_hybrid"

    theta_avg_hist = Dict{Int, Vector{Float64}}()
    v_mag_avg_hist = Dict{Int, Vector{Float64}}()

    results = Vector{NamedTuple{(:expr_cost, :real_cost), Tuple{Float64, Float64}}}(
        undef,
        n_zones,
    )

    init_ac_hybrid_partition!(
        partitions,
        sys,
        optimizer_factory,
        rho_theta[],
        rho_v[],
        rho_p[],
        rho_q[],
        interzonal_branches,
    )

    iter_done = 0
    for iter in 1:max_iter
        iter_done = iter
        _FAST_USE_LOGS && println("\n[ADMM-AC-hybrid] Iteração $iter")

        refresh_rho_bus_ac!(partitions, rho_theta, prev_rho_theta, rho_v, prev_rho_v)
        refresh_rho_branch_ac!(partitions, rho_p, prev_rho_p, rho_q, prev_rho_q)

        Threads.@threads for i in 1:n_zones
            zone_id = zone_ids[i]
            part = partitions.parts[i]
            model = part.model::JuMP.Model

            θbuf = dual_theta_buf[zone_id]
            θbuf .= duals_theta[zone_id]
            Vbuf = dual_v_mag_buf[zone_id]
            Vbuf .= duals_v_mag[zone_id]

            lin_bus = build_linear_bus_term(
                part,
                shared_buses_sorted,
                shared_buses_to_pos,
                theta_consensus,
                v_mag_consensus,
                θbuf,
                Vbuf,
                rho_theta[],
                rho_v[],
            )

            lin_br = build_linear_branch_term(
                part,
                part.zone_edges,
                p_flow_consensus,
                q_flow_consensus,
                duals_p_flow[zone_id],
                duals_q_flow[zone_id],
                rho_p[],
                rho_q[],
            )

            @objective(
                model,
                Min,
                part.obj_base + part.quad_bus_expr_ac + part.quad_branch_expr_ac + lin_bus +
                lin_br
            )
            optimize!(model)

            results[i] = (
                expr_cost = value(model[:generation_cost]),
                real_cost = compute_real_generation_cost(model, part.objective_params),
            )
        end

        @inbounds for i in 1:n_zones
            part = partitions.parts[i]
            out = results[i]
            part.has_solution = true
            part.generation_cost_model_expr = out.expr_cost
            part.real_generation_cost = out.real_cost
        end

        theta_consensus_prev .= theta_consensus
        v_mag_consensus_prev .= v_mag_consensus
        p_flow_consensus_prev = copy(p_flow_consensus)
        q_flow_consensus_prev = copy(q_flow_consensus)

        compute_ac_bus_consensus_vec!(
            theta_consensus,
            v_mag_consensus,
            partitions,
            shared_buses_sorted,
            shared_buses_to_pos,
        )
        compute_ac_branch_consensus!(
            p_flow_consensus,
            q_flow_consensus,
            partitions,
            interzonal_branches,
        )

        @inbounds for b in shared_buses_sorted
            j = shared_buses_to_pos[b]
            push!(get!(theta_avg_hist, b, Float64[]), theta_consensus[j])
            push!(get!(v_mag_avg_hist, b, Float64[]), v_mag_consensus[j])
        end

        update_duals_bus_ac_all!(
            duals_theta,
            duals_v_mag,
            partitions,
            zone_ids,
            theta_consensus,
            v_mag_consensus,
            rho_theta,
            rho_v,
            shared_buses_sorted,
            shared_buses_to_pos;
            theta_consensus_prev = theta_consensus_prev,
            v_mag_consensus_prev = v_mag_consensus_prev,
        )

        update_duals_branch_ac!(
            duals_p_flow,
            duals_q_flow,
            partitions,
            p_flow_consensus,
            q_flow_consensus,
            interzonal_branches,
            rho_p[],
            rho_q[];
            p_flow_consensus_prev = p_flow_consensus_prev,
            q_flow_consensus_prev = q_flow_consensus_prev,
        )

        r_bus = primal_residual_bus_ac_vec(
            partitions,
            theta_consensus,
            v_mag_consensus,
            shared_buses_sorted,
            shared_buses_to_pos,
        )
        r_br = primal_residual_branch_ac(
            partitions,
            p_flow_consensus,
            q_flow_consensus,
            interzonal_branches,
        )
        r_tot = hypot(r_bus, r_br)

        s_bus = dual_residual_bus_ac_vec(
            theta_consensus, v_mag_consensus,
            theta_consensus_prev, v_mag_consensus_prev,
            rho_theta[], rho_v[],
        )
        s_br = dual_residual_branch_ac(
            p_flow_consensus, q_flow_consensus,
            p_flow_consensus_prev, q_flow_consensus_prev,
            interzonal_branches,
            rho_p[], rho_q[],
        )
        s_tot = hypot(s_bus, s_br)

        eps_bus_r = eps_primal_bus_ac_vec(
            partitions,
            theta_consensus,
            v_mag_consensus,
            shared_buses_sorted,
            shared_buses_to_pos,
            eps_abs,
            eps_rel,
        )
        eps_br_r = eps_primal_branch_ac(
            partitions,
            p_flow_consensus,
            q_flow_consensus,
            interzonal_branches,
            eps_abs,
            eps_rel,
        )
        eps_r_tot = hypot(eps_bus_r, eps_br_r)

        eps_bus_s = eps_dual_bus_ac_unscaled(
            duals_theta,
            duals_v_mag,
            partitions,
            zone_ids,
            shared_buses_to_pos,
            eps_abs,
            eps_rel,
        )
        eps_br_s = eps_dual_branch_ac_unscaled(
            duals_p_flow,
            duals_q_flow,
            zone_ids,
            interzonal_branches,
            eps_abs,
            eps_rel,
        )
        eps_s_tot = hypot(eps_bus_s, eps_br_s)

        _record_iter!(
            hist, r_tot, s_tot, eps_r_tot, eps_s_tot;
            rhos = Dict(
                :rho_theta => rho_theta[],
                :rho_v => rho_v[],
                :rho_p => rho_p[],
                :rho_q => rho_q[],
            ),
        )

        _FAST_USE_LOGS && println("‖r‖ = $r_tot | eps_pri = $eps_r_tot")
        _FAST_USE_LOGS && println("‖s‖ = $s_tot | eps_du  = $eps_s_tot")

        if stop_if_converged!(r_tot, s_tot, eps_r_tot, eps_s_tot)
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        maybe_adapt_rho!(
            adaptive,
            rho_theta,
            rho_v,
            partitions,
            theta_consensus,
            v_mag_consensus,
            theta_consensus_prev,
            v_mag_consensus_prev,
            shared_buses_sorted,
            shared_buses_to_pos,
            rho_p,
            rho_q,
            r_br,
            s_br,
            config,
        )
    end

    hist.heatmaps[:theta_avg] = theta_avg_hist
    hist.heatmaps[:v_magavg] = v_mag_avg_hist

    result = _finalize_admm_result(partitions, t0_ns, iter_done)
    _print_admm_summary!(
        "[ADMM-hybrid/AC]",
        result;
        rho_lines = (
            "rho_theta inicial: $(config.rho_theta) | rho_theta final: $(rho_theta[])",
            "rho_v     inicial: $(config.rho_v)     | rho_v     final: $(rho_v[])",
            "rhoP      inicial: $(config.rho_p)     | rhoP      final: $(rho_p[])",
            "rhoQ      inicial: $(config.rho_q)     | rhoQ      final: $(rho_q[])",
        ),
    )

    return LirionRunResult(; result, model = nothing, system = sys, diagnostics = hist)
end
