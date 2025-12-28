# ============================================================
# ADMM-SOCP
# ============================================================

# ------------------------------------------------------------
# Partition initialization
# ------------------------------------------------------------

# ----------------------------
# SOCP / V² + (R,T) partition
# ----------------------------

function init_socp_v2rt_partition!(
    sub::PartitionSOCP,
    sys_zone,
    optimizer_factory,
    rhoV2::Float64,
    rhoR::Float64,
    rhoT::Float64,
    cut_edges,
)
    prepare_hotpaths!(sub)
    sub.zone_edges = _zone_cut_edges(sub, cut_edges)

    m = build_model(sys_zone, SOCP; optimizer = optimizer_factory)
    sub.model = m
    sub.obj_base = m[:generation_cost]
    sub.v2_var = m[:v_magnitude_sq]
    sub.objective_params = sys_zone.objective
    sub.has_solution = false
    _set_model_silent!(m)

    for I in keys(sub.v2_var)
        set_start_value(sub.v2_var[I], 1.0)
    end

    _build_quad_bus_expr_socp!(sub, rhoV2)
    _build_quad_branch_rt_expr!(sub, rhoR, rhoT, sub.zone_edges)

    sub.obj_static = sub.obj_base + sub.quad_bus_expr_socp + sub.quad_rt_expr_socp
    @objective(m, Min, sub.obj_static)
    optimize!(m)

    return nothing
end

function init_socp_v2rt_partition!(
    partitions::Partitions{PartitionSOCP},
    sys::ADMMSystem,
    zone_ids::Vector{Int},
    optimizer_factory,
    rhoV2::Float64,
    rhoR::Float64,
    rhoT::Float64,
    cut_edges,
)
    for zone_id in zone_ids
        init_socp_v2rt_partition!(
            getzone(partitions, zone_id),
            sys.parts[zone_id],
            optimizer_factory,
            rhoV2,
            rhoR,
            rhoT,
            cut_edges,
        )
    end
    return nothing
end

# ----------------------------
# SOCP / (P,Q) partition
# ----------------------------

function init_socp_pq_partition!(
    sub::PartitionSOCP,
    sys_zone,
    optimizer_factory,
    rhoP::Float64,
    rhoQ::Float64,
    cut_edges,
)
    prepare_hotpaths!(sub)
    sub.zone_edges = _zone_cut_edges(sub, cut_edges)

    m = build_model(sys_zone, SOCP; optimizer = optimizer_factory)
    sub.model = m
    sub.obj_base = m[:generation_cost]
    sub.objective_params = sys_zone.objective
    sub.has_solution = false
    _set_model_silent!(m)

    if is_gurobi(m)
        _soft_anchor_branch_expr!(
            sub;
            κV2 = 1e-6,
            κPQ = 1e-8,
            use_bus = true,
            use_pq = true,
        )
    else
        sub.soft_anchor_expr = zero(JuMP.QuadExpr)
    end

    _build_quad_branch_pq_expr!(sub, rhoP, rhoQ, sub.zone_edges)
    sub.obj_static = sub.obj_base + sub.quad_pq_expr_socp + sub.soft_anchor_expr
    @objective(m, Min, sub.obj_static)
    optimize!(m)

    return nothing
end

function init_socp_pq_partition!(
    partitions::Partitions{PartitionSOCP},
    sys::ADMMSystem,
    zone_ids::Vector{Int},
    optimizer_factory,
    rhoP::Float64,
    rhoQ::Float64,
    cut_edges,
)
    for zone_id in zone_ids
        init_socp_pq_partition!(
            getzone(partitions, zone_id),
            sys.parts[zone_id],
            optimizer_factory,
            rhoP,
            rhoQ,
            cut_edges,
        )
    end
    return nothing
end

# ----------------------------
# SOCP / hybrid partition (V² + P/Q + R/T)
# ----------------------------

function init_socp_hybrid_partition!(
    sub::PartitionSOCP,
    sys_zone,
    optimizer_factory,
    rhoV2::Float64,
    rhoP::Float64,
    rhoQ::Float64,
    rhoR::Float64,
    rhoT::Float64,
    cut_edges,
)
    prepare_hotpaths!(sub)
    sub.zone_edges = _zone_cut_edges(sub, cut_edges)

    m = build_model(sys_zone, SOCP; optimizer = optimizer_factory)
    sub.model = m
    sub.obj_base = m[:generation_cost]
    sub.v2_var = m[:v_magnitude_sq]
    sub.objective_params = sys_zone.objective
    sub.has_solution = false
    _set_model_silent!(m)

    for I in keys(sub.v2_var)
        set_start_value(sub.v2_var[I], 1.0)
    end

    _build_quad_bus_expr_socp!(sub, rhoV2)
    _build_quad_branch_pq_expr!(sub, rhoP, rhoQ, sub.zone_edges)
    _build_quad_branch_rt_expr!(sub, rhoR, rhoT, sub.zone_edges)

    sub.obj_static =
        sub.obj_base + sub.quad_bus_expr_socp + sub.quad_pq_expr_socp +
        sub.quad_rt_expr_socp
    @objective(m, Min, sub.obj_static)
    optimize!(m)

    return nothing
end

function init_socp_hybrid_partition!(
    partitions::Partitions{PartitionSOCP},
    sys::ADMMSystem,
    zone_ids::Vector{Int},
    optimizer_factory,
    rhoV2::Float64,
    rhoP::Float64,
    rhoQ::Float64,
    rhoR::Float64,
    rhoT::Float64,
    cut_edges,
)
    for zone_id in zone_ids
        init_socp_hybrid_partition!(
            getzone(partitions, zone_id),
            sys.parts[zone_id],
            optimizer_factory,
            rhoV2,
            rhoP,
            rhoQ,
            rhoR,
            rhoT,
            cut_edges,
        )
    end
    return nothing
end

# ============================================================
# RUNNERS (ADMM-SOCP)
# ============================================================

# ------------------------------------------------------------
# Runner: ADMM SOCP (V² + R/T)
# ------------------------------------------------------------

function run_admm_socp(sys::ADMMSystem, config::DistributedOPFModel)::LirionRunResult
    t0_ns = time_ns()
    @assert config.base_model == SOCP "ADMM-SOCP: base_model deve ser SOCP."

    shared_buses = sys.shared_buses
    cut_edges = sys.cut_edges
    max_iter = config.max_iter

    rhoV2 = Ref(config.rho_v)
    rhoR = Ref(config.rho_flow_p)
    rhoT = Ref(config.rho_flow_q)

    prev_rhoV2 = Ref(rhoV2[])
    prev_rhoR = Ref(rhoR[])
    prev_rhoT = Ref(rhoT[])

    eps_abs = config.eps_abs
    eps_rel = config.eps_rel
    adaptive = _adaptive_on(config)

    optimizer_factory = resolve_socp_optimizer(config)

    partitions = build_socp_partitions(sys)     # ::Partitions{PartitionSOCP}
    zone_ids = partitions.zones
    n_zones = length(zone_ids)

    shared_buses_sorted, shared_buses_to_pos = sorted_ids_and_id2idx(shared_buses)
    n_shared_buses = length(shared_buses_sorted)

    v2bar = ones(Float64, n_shared_buses)
    v2bar_prev = similar(v2bar)

    duals_v2 = Dict{Int, Vector{Float64}}(
        zone_id => zeros(Float64, n_shared_buses) for zone_id in zone_ids
    )
    lambda_zone_buf = Dict{Int, Vector{Float64}}(
        zone_id => zeros(Float64, n_shared_buses) for zone_id in zone_ids
    )

    zR = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zT = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zR_prev = Dict{Int, Float64}()
    zT_prev = Dict{Int, Float64}()

    dualsR = Dict{Int, Dict{Int, Float64}}(
        zone_id => Dict{Int, Float64}() for zone_id in zone_ids
    )
    dualsT = Dict{Int, Dict{Int, Float64}}(
        zone_id => Dict{Int, Float64}() for zone_id in zone_ids
    )

    # ---------- diagnostics ----------
    hist = ADMMHistory()
    hist.meta[:kind] = :socp_v2rt
    hist.meta[:tag] = "admm_socp_v2rt"

    # opcional: histórico por barra (consenso de V²)
    v2_avg_history = Dict{Int, Vector{Float64}}()

    results = Vector{NamedTuple{(:gen_cost_expr, :real_cost), Tuple{Float64, Float64}}}(
        undef, n_zones,
    )

    init_socp_v2rt_partition!(
        partitions, sys, zone_ids, optimizer_factory,
        rhoV2[], rhoR[], rhoT[], cut_edges,
    )

    iter_done = 0
    for iter in 1:max_iter
        iter_done = iter
        _FAST_USE_LOGS && println("\n[ADMM-SOCP (V²+R,T)] Iteração $iter")

        refresh_rho_v2_socp!(partitions, rhoV2, prev_rhoV2)
        refresh_rho_rt_socp!(partitions, rhoR, prev_rhoR, rhoT, prev_rhoT)

        α_used = nothing
        v2_target = v2bar

        Threads.@threads for i in 1:n_zones
            zone_id = zone_ids[i]
            sub = getzone(partitions, zone_id)
            m = sub.model::JuMP.Model

            lambda_z = (lambda_zone_buf[zone_id] .= duals_v2[zone_id])

            lin_bus = _build_linear_bus_expr_socp(
                sub, shared_buses_sorted, shared_buses_to_pos, v2_target, lambda_z,
                rhoV2[],
            )
            lin_rt = _build_linear_branch_rt_expr(
                sub, sub.zone_edges, zR, zT, dualsR[zone_id], dualsT[zone_id], rhoR[],
                rhoT[],
            )

            @objective(m, Min, sub.obj_static + (lin_bus + lin_rt))
            optimize!(m)

            results[i] = (
                gen_cost_expr = value(m[:generation_cost]),
                real_cost = compute_real_generation_cost(m, sub.objective_params),
            )
        end

        for (i, zone_id) in enumerate(zone_ids)
            sub = getzone(partitions, zone_id)
            out = results[i]
            sub.has_solution = true
            sub.generation_cost_model_expr = out.gen_cost_expr
            sub.real_generation_cost = out.real_cost
        end

        v2bar_prev .= v2bar
        zR_prev = copy(zR)
        zT_prev = copy(zT)

        compute_v2_consensus_vec!(
            v2bar,
            partitions,
            shared_buses_sorted,
            shared_buses_to_pos,
        )
        compute_rt_consensus!(zR, zT, partitions, cut_edges)

        # histórico por barra (opcional)
        @inbounds for b in shared_buses_sorted
            j = shared_buses_to_pos[b]
            push!(get!(v2_avg_history, b, Float64[]), v2bar[j])
        end

        update_duals_v2_all!(
            duals_v2,
            partitions,
            zone_ids,
            v2bar,
            v2bar_prev,
            rhoV2[],
            α_used,
            shared_buses_sorted,
            shared_buses_to_pos,
        )
        update_duals_rt_all!(dualsR, dualsT, partitions, zone_ids, zR, zT, rhoR[], rhoT[])

        r_bus = primal_residual_bus_socp_vec(
            partitions, v2bar, shared_buses_sorted, shared_buses_to_pos,
        )
        r_rt = primal_residual_branch_rt_socp(partitions, zR, zT, cut_edges)
        r = hypot(r_bus, r_rt)

        s_bus = dual_residual_bus_socp_vec(v2bar, v2bar_prev)
        s_rt = dual_residual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges)
        s = hypot(s_bus, s_rt)

        eps_r = begin
            epsb = eps_primal_bus_socp_vec(
                partitions, v2bar, shared_buses_sorted, shared_buses_to_pos,
                eps_abs,
                eps_rel,
            )
            epsrt = eps_primal_branch_rt_socp(
                partitions, zR, zT, cut_edges, eps_abs, eps_rel,
            )
            hypot(epsb, epsrt)
        end
        eps_s = begin
            epsb = eps_dual_bus_socp_vec(v2bar, v2bar_prev, eps_abs, eps_rel)
            epsrt = eps_dual_branch_rt_socp(
                zR, zT, zR_prev, zT_prev, cut_edges, eps_abs, eps_rel,
            )
            hypot(epsb, epsrt)
        end

        _record_iter!(
            hist, r, s, eps_r, eps_s;
            rhos = Dict(
                :rho_v2 => rhoV2[],
                :rho_r => rhoR[],
                :rho_t => rhoT[],
            ),
        )

        _FAST_USE_LOGS && println("‖r‖ = $r | eps_pri = $eps_r")
        _FAST_USE_LOGS && println("‖s‖ = $s | eps_du  = $eps_s")

        if stop_if_converged!(r, s, eps_r, eps_s)
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        maybe_adapt_rho!(adaptive, rhoV2, r_bus, s_bus, rhoR, rhoT, r_rt, s_rt, config)
    end

    hist.heatmaps[:v2_avg] = v2_avg_history

    total_real_cost = sum(p.real_generation_cost for p in partitions.parts)
    total_expr_cost = sum(p.generation_cost_model_expr for p in partitions.parts)

    elapsed_s = (time_ns() - t0_ns) / 1e9
    println("\n[ADMM-SOCP (V²+R,T)] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("rho(V²) inicial: $(config.rho_v) | rho(V²) final: $(rhoV2[])")
    println("rhoR    inicial: $(config.rho_flow_p) | rhoR final: $(rhoR[])")
    println("rhoT    inicial: $(config.rho_flow_q) | rhoT final: $(rhoT[])")
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (SOCP V²+RT): $(_pretty_time(elapsed_s))")

    result = (
        total_real_cost = total_real_cost,
        total_expr_cost = total_expr_cost,
        partitions = partitions,
        iter_done = iter_done,
        elapsed_s = elapsed_s,
    )
    return LirionRunResult(;
        result = result,
        model = nothing,
        system = sys,
        diagnostics = hist,
    )
end

# ------------------------------------------------------------
# Runner: ADMM SOCP (branch P/Q)
# ------------------------------------------------------------

function run_admm_socp_branch(sys::ADMMSystem, config::DistributedOPFModel)::LirionRunResult
    t0_ns = time_ns()
    @assert config.base_model == SOCP "ADMM-SOCP-branch: base_model deve ser SOCP."

    cut_edges = sys.cut_edges
    max_iter = config.max_iter

    rhoP = Ref(config.rho_flow_p)
    rhoQ = Ref(config.rho_flow_q)
    prev_rhoP = Ref(rhoP[])
    prev_rhoQ = Ref(rhoQ[])

    eps_abs = config.eps_abs
    eps_rel = config.eps_rel
    adaptive = _adaptive_on(config)

    optimizer_factory = resolve_socp_optimizer(config)

    partitions = build_socp_partitions(sys)   # ::Partitions{PartitionSOCP}
    zone_ids = partitions.zones
    n_zones = length(zone_ids)

    zP = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zQ = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zP_prev = Dict{Int, Float64}()
    zQ_prev = Dict{Int, Float64}()

    dualsP = Dict{Int, Dict{Int, Float64}}(
        zone_id => Dict{Int, Float64}() for zone_id in zone_ids
    )
    dualsQ = Dict{Int, Dict{Int, Float64}}(
        zone_id => Dict{Int, Float64}() for zone_id in zone_ids
    )

    # ---------- diagnostics ----------
    hist = ADMMHistory()
    hist.meta[:kind] = :socp_pq
    hist.meta[:tag] = "admm_socp_pq"

    init_socp_pq_partition!(
        partitions, sys, zone_ids, optimizer_factory, rhoP[], rhoQ[], cut_edges,
    )

    iter_done = 0
    for iter in 1:max_iter
        iter_done = iter
        _FAST_USE_LOGS && println("\n[ADMM-SOCP-branch (P,Q)] Iteração $iter")

        refresh_rho_pq_socp!(partitions, rhoP, prev_rhoP, rhoQ, prev_rhoQ)

        Threads.@threads for i in 1:n_zones
            zone_id = zone_ids[i]
            sub = getzone(partitions, zone_id)
            m = sub.model::JuMP.Model

            lin = _build_linear_branch_pq_expr(
                sub, sub.zone_edges, zP, zQ, dualsP[zone_id], dualsQ[zone_id], rhoP[],
                rhoQ[],
            )

            @objective(m, Min, sub.obj_static + lin)
            optimize!(m)

            sub.has_solution = true
            sub.generation_cost_model_expr = value(m[:generation_cost])
            sub.real_generation_cost = compute_real_generation_cost(m, sub.objective_params)
        end

        zP_prev = copy(zP)
        zQ_prev = copy(zQ)

        compute_pq_consensus!(zP, zQ, partitions, cut_edges)

        update_duals_pq_all!(dualsP, dualsQ, partitions, zone_ids, zP, zQ, rhoP[], rhoQ[])

        r = primal_residual_branch_pq_socp(partitions, zP, zQ, cut_edges)
        s = dual_residual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges)

        eps_r = eps_primal_branch_pq_socp(partitions, zP, zQ, cut_edges, eps_abs, eps_rel)
        eps_s =
            eps_dual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges, eps_abs, eps_rel)

        _record_iter!(
            hist, r, s, eps_r, eps_s;
            rhos = Dict(:rho_p => rhoP[], :rho_q => rhoQ[]),
        )

        _FAST_USE_LOGS && println("‖r‖ = $r | eps_pri = $eps_r")
        _FAST_USE_LOGS && println("‖s‖ = $s | eps_du  = $eps_s")

        if stop_if_converged!(r, s, eps_r, eps_s)
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        maybe_adapt_rho!(adaptive, rhoP, rhoQ, r, s, config)
    end

    total_real_cost = sum(p.real_generation_cost for p in partitions.parts)
    total_expr_cost = sum(p.generation_cost_model_expr for p in partitions.parts)

    elapsed_s = (time_ns() - t0_ns) / 1e9
    println("\n[ADMM-SOCP-branch (P,Q)] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("rhoP inicial: $(config.rho_flow_p) | rhoP final: $(rhoP[])")
    println("rhoQ inicial: $(config.rho_flow_q) | rhoQ final: $(rhoQ[])")
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (SOCP P/Q): $(_pretty_time(elapsed_s))")

    result = (
        total_real_cost = total_real_cost,
        total_expr_cost = total_expr_cost,
        partitions = partitions,
        iter_done = iter_done,
        elapsed_s = elapsed_s,
    )
    return LirionRunResult(;
        result = result,
        model = nothing,
        system = sys,
        diagnostics = hist,
    )
end

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

@inline function _get_rho_rt_initial(config::DistributedOPFModel)
    # fallback: se não existir rho_flow_r/t, usa p/q
    ρr = try
        getfield(config, :rho_flow_r)
    catch
        config.rho_flow_p
    end
    ρt = try
        getfield(config, :rho_flow_t)
    catch
        config.rho_flow_q
    end
    return ρr, ρt
end

# ------------------------------------------------------------
# Runner: ADMM SOCP hybrid (V² + P/Q + R/T)
# ------------------------------------------------------------

function run_admm_socp_hybrid(sys::ADMMSystem, config::DistributedOPFModel)::LirionRunResult
    t0_ns = time_ns()
    @assert config.base_model == SOCP "ADMM-SOCP-hybrid: base_model deve ser SOCP."

    shared_buses = sys.shared_buses
    cut_edges = sys.cut_edges
    max_iter = config.max_iter

    rhoV2 = Ref(config.rho_v)

    rhoP = Ref(config.rho_flow_p)
    rhoQ = Ref(config.rho_flow_q)

    ρr0, ρt0 = _get_rho_rt_initial(config)
    rhoR = Ref(ρr0)
    rhoT = Ref(ρt0)

    prev_rhoV2 = Ref(rhoV2[])
    prev_rhoP = Ref(rhoP[])
    prev_rhoQ = Ref(rhoQ[])
    prev_rhoR = Ref(rhoR[])
    prev_rhoT = Ref(rhoT[])

    eps_abs = config.eps_abs
    eps_rel = config.eps_rel
    adaptive = _adaptive_on(config)

    optimizer_factory = resolve_socp_optimizer(config)

    partitions = build_socp_partitions(sys)   # ::Partitions{PartitionSOCP}
    zone_ids = partitions.zones
    n_zones = length(zone_ids)

    shared_buses_sorted, shared_buses_to_pos = sorted_ids_and_id2idx(shared_buses)
    n_shared_buses = length(shared_buses_sorted)

    v2bar = ones(Float64, n_shared_buses)
    v2bar_prev = similar(v2bar)

    duals_v2 = Dict{Int, Vector{Float64}}(
        zone_id => zeros(Float64, n_shared_buses) for zone_id in zone_ids
    )
    lambda_zone_buf = Dict{Int, Vector{Float64}}(
        zone_id => zeros(Float64, n_shared_buses) for zone_id in zone_ids
    )

    zP = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zQ = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zR = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zT = Dict{Int, Float64}(e => 0.0 for e in cut_edges)

    zP_prev = Dict{Int, Float64}()
    zQ_prev = Dict{Int, Float64}()
    zR_prev = Dict{Int, Float64}()
    zT_prev = Dict{Int, Float64}()

    dualsP = Dict{Int, Dict{Int, Float64}}(
        zone_id => Dict{Int, Float64}() for zone_id in zone_ids
    )
    dualsQ = Dict{Int, Dict{Int, Float64}}(
        zone_id => Dict{Int, Float64}() for zone_id in zone_ids
    )
    dualsR = Dict{Int, Dict{Int, Float64}}(
        zone_id => Dict{Int, Float64}() for zone_id in zone_ids
    )
    dualsT = Dict{Int, Dict{Int, Float64}}(
        zone_id => Dict{Int, Float64}() for zone_id in zone_ids
    )

    # ---------- diagnostics ----------
    hist = ADMMHistory()
    hist.meta[:kind] = :socp_hybrid
    hist.meta[:tag] = "admm_socp_hybrid"
    v2_avg_history = Dict{Int, Vector{Float64}}()

    init_socp_hybrid_partition!(
        partitions, sys, zone_ids, optimizer_factory,
        rhoV2[], rhoP[], rhoQ[], rhoR[], rhoT[], cut_edges,
    )

    iter_done = 0
    for iter in 1:max_iter
        iter_done = iter
        _FAST_USE_LOGS && println("\n[ADMM-SOCP-hybrid] Iteração $iter")

        refresh_rho_hybrid_socp!(
            partitions,
            rhoV2, prev_rhoV2,
            rhoP, prev_rhoP,
            rhoQ, prev_rhoQ,
            rhoR, prev_rhoR,
            rhoT, prev_rhoT,
        )

        α_used = nothing
        v2_target = v2bar

        Threads.@threads for i in 1:n_zones
            zone_id = zone_ids[i]
            sub = getzone(partitions, zone_id)
            m = sub.model::JuMP.Model

            lambda_z = (lambda_zone_buf[zone_id] .= duals_v2[zone_id])

            lin_bus = _build_linear_bus_expr_socp(
                sub, shared_buses_sorted, shared_buses_to_pos, v2_target, lambda_z,
                rhoV2[],
            )
            lin_pq = _build_linear_branch_pq_expr(
                sub, sub.zone_edges, zP, zQ, dualsP[zone_id], dualsQ[zone_id], rhoP[],
                rhoQ[],
            )
            lin_rt = _build_linear_branch_rt_expr(
                sub, sub.zone_edges, zR, zT, dualsR[zone_id], dualsT[zone_id], rhoR[],
                rhoT[],
            )

            @objective(m, Min, sub.obj_static + (lin_bus + lin_pq + lin_rt))
            optimize!(m)

            sub.has_solution = true
            sub.generation_cost_model_expr = value(m[:generation_cost])
            sub.real_generation_cost = compute_real_generation_cost(m, sub.objective_params)
        end

        v2bar_prev .= v2bar
        zP_prev = copy(zP)
        zQ_prev = copy(zQ)
        zR_prev = copy(zR)
        zT_prev = copy(zT)

        compute_v2_consensus_vec!(
            v2bar,
            partitions,
            shared_buses_sorted,
            shared_buses_to_pos,
        )
        compute_pq_consensus!(zP, zQ, partitions, cut_edges)
        compute_rt_consensus!(zR, zT, partitions, cut_edges)

        @inbounds for b in shared_buses_sorted
            j = shared_buses_to_pos[b]
            push!(get!(v2_avg_history, b, Float64[]), v2bar[j])
        end

        update_duals_v2_all!(
            duals_v2, partitions, zone_ids, v2bar, v2bar_prev, rhoV2[],
            α_used, shared_buses_sorted, shared_buses_to_pos,
        )
        update_duals_pq_all!(dualsP, dualsQ, partitions, zone_ids, zP, zQ, rhoP[], rhoQ[])
        update_duals_rt_all!(dualsR, dualsT, partitions, zone_ids, zR, zT, rhoR[], rhoT[])

        r_bus = primal_residual_bus_socp_vec(
            partitions, v2bar, shared_buses_sorted, shared_buses_to_pos,
        )
        r_pq = primal_residual_branch_pq_socp(partitions, zP, zQ, cut_edges)
        r_rt = primal_residual_branch_rt_socp(partitions, zR, zT, cut_edges)
        r = sqrt(r_bus^2 + r_pq^2 + r_rt^2)

        s_bus = dual_residual_bus_socp_vec(v2bar, v2bar_prev)
        s_pq = dual_residual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges)
        s_rt = dual_residual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges)
        s = sqrt(s_bus^2 + s_pq^2 + s_rt^2)

        eps_r = begin
            epsb = eps_primal_bus_socp_vec(
                partitions, v2bar, shared_buses_sorted, shared_buses_to_pos,
                eps_abs,
                eps_rel,
            )
            epspq = eps_primal_branch_pq_socp(
                partitions,
                zP,
                zQ,
                cut_edges,
                eps_abs,
                eps_rel,
            )
            epsrt = eps_primal_branch_rt_socp(
                partitions,
                zR,
                zT,
                cut_edges,
                eps_abs,
                eps_rel,
            )
            sqrt(epsb^2 + epspq^2 + epsrt^2)
        end
        eps_s = begin
            epsb = eps_dual_bus_socp_vec(v2bar, v2bar_prev, eps_abs, eps_rel)
            epspq = eps_dual_branch_pq_socp(
                zP,
                zQ,
                zP_prev,
                zQ_prev,
                cut_edges,
                eps_abs,
                eps_rel,
            )
            epsrt = eps_dual_branch_rt_socp(
                zR,
                zT,
                zR_prev,
                zT_prev,
                cut_edges,
                eps_abs,
                eps_rel,
            )
            sqrt(epsb^2 + epspq^2 + epsrt^2)
        end

        _record_iter!(
            hist, r, s, eps_r, eps_s;
            rhos = Dict(
                :rho_v2 => rhoV2[],
                :rho_p => rhoP[],
                :rho_q => rhoQ[],
                :rho_r => rhoR[],
                :rho_t => rhoT[],
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
            rhoV2, r_bus, s_bus,
            rhoP, rhoQ, r_pq, s_pq,
            rhoR, rhoT, r_rt, s_rt,
            config,
        )
    end

    hist.heatmaps[:v2_avg] = v2_avg_history

    total_real_cost = sum(p.real_generation_cost for p in partitions.parts)
    total_expr_cost = sum(p.generation_cost_model_expr for p in partitions.parts)

    elapsed_s = (time_ns() - t0_ns) / 1e9
    println("\n[ADMM-SOCP-hybrid] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("rho(V²) inicial: $(config.rho_v) | rho(V²) final: $(rhoV2[])")
    println("rhoP    inicial: $(config.rho_flow_p) | rhoP final: $(rhoP[])")
    println("rhoQ    inicial: $(config.rho_flow_q) | rhoQ final: $(rhoQ[])")
    println("rhoR    inicial: $(ρr0) | rhoR final: $(rhoR[])")
    println("rhoT    inicial: $(ρt0) | rhoT final: $(rhoT[])")
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (SOCP híbrido): $(_pretty_time(elapsed_s))")

    result = (
        total_real_cost = total_real_cost,
        total_expr_cost = total_expr_cost,
        partitions = partitions,
        iter_done = iter_done,
        elapsed_s = elapsed_s,
    )
    return LirionRunResult(;
        result = result,
        model = nothing,
        system = sys,
        diagnostics = hist,
    )
end
