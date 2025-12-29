# ============================================================
# ADMM-AC
# ============================================================

# ------------------------------------------------------------
# Partition initialization
# ------------------------------------------------------------

# ----------------------------
# AC / bus partition
# ----------------------------

function init_ac_bus_partition!(
    part::PartitionAC,
    sys_zone,
    optimizer_factory,
    rho_theta::Float64,
    rhoV::Float64,
)
    prepare_hotpaths!(part)

    model = build_model(sys_zone, AC; optimizer = optimizer_factory)
    part.model = model
    part.obj_base = model[:generation_cost]
    part.theta_var = model[:theta]
    part.vmag_var = model[:v_magnitude]
    part.objective_params = sys_zone.objective
    part.has_solution = false

    for I in keys(part.theta_var)
        set_start_value(part.theta_var[I], 0.0)
    end
    for I in keys(part.vmag_var)
        set_start_value(part.vmag_var[I], 1.0)
    end

    _build_quad_bus_expr_ac!(part, rho_theta, rhoV)
    @objective(model, Min, part.obj_base + part.quad_bus_expr_ac)

    return nothing
end

function init_ac_bus_partition!(
    partitions::Partitions{PartitionAC},
    sys::ADMMSystem,
    zone_ids::Vector{Int},
    optimizer_factory,
    rho_theta::Float64,
    rhoV::Float64,
)
    for zone_id in zone_ids
        init_ac_bus_partition!(
            get_partition(partitions, zone_id),
            sys.parts[zone_id],
            optimizer_factory,
            rho_theta,
            rhoV,
        )
    end
    return nothing
end

# ----------------------------
# AC / branch partition
# ----------------------------

function init_ac_branch_partition!(
    part::PartitionAC,
    sys_zone,
    optimizer_factory,
    rhoP::Float64,
    rhoQ::Float64,
    cut_edges,
)
    prepare_hotpaths!(part)
    part.zone_edges = _zone_cut_edges(part, cut_edges)

    model = build_model(sys_zone, AC; optimizer = optimizer_factory)
    part.model = model
    part.obj_base = model[:generation_cost]
    part.theta_var = model[:theta]
    part.vmag_var = model[:v_magnitude]
    part.objective_params = sys_zone.objective
    part.has_solution = false

    for I in keys(part.theta_var)
        set_start_value(part.theta_var[I], 0.0)
    end
    for I in keys(part.vmag_var)
        set_start_value(part.vmag_var[I], 1.0)
    end

    _build_quad_branch_expr_ac!(part, rhoP, rhoQ, part.zone_edges)
    @objective(model, Min, part.obj_base + part.quad_branch_expr_ac)

    return nothing
end

function init_ac_branch_partition!(
    partitions::Partitions{PartitionAC},
    sys::ADMMSystem,
    zone_ids::Vector{Int},
    optimizer_factory,
    rhoP::Float64,
    rhoQ::Float64,
    cut_edges,
)
    for zone_id in zone_ids
        init_ac_branch_partition!(
            get_partition(partitions, zone_id),
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
# AC / hybrid partition
# ----------------------------

function init_ac_hybrid_partition!(
    part::PartitionAC,
    sys_zone,
    optimizer_factory,
    rho_theta::Float64,
    rhoV::Float64,
    rhoP::Float64,
    rhoQ::Float64,
    cut_edges,
)
    prepare_hotpaths!(part)
    part.zone_edges = _zone_cut_edges(part, cut_edges)

    model = build_model(sys_zone, AC; optimizer = optimizer_factory)
    part.model = model
    part.obj_base = model[:generation_cost]
    part.theta_var = model[:theta]
    part.vmag_var = model[:v_magnitude]
    part.objective_params = sys_zone.objective
    part.has_solution = false

    for I in keys(part.theta_var)
        set_start_value(part.theta_var[I], 0.0)
    end
    for I in keys(part.vmag_var)
        set_start_value(part.vmag_var[I], 1.0)
    end

    _build_quad_bus_expr_ac!(part, rho_theta, rhoV)
    _build_quad_branch_expr_ac!(part, rhoP, rhoQ, part.zone_edges)

    @objective(model, Min, part.obj_base + part.quad_bus_expr_ac + part.quad_branch_expr_ac)

    return nothing
end

function init_ac_hybrid_partition!(
    partitions::Partitions{PartitionAC},
    sys::ADMMSystem,
    zone_ids::Vector{Int},
    optimizer_factory,
    rho_theta::Float64,
    rhoV::Float64,
    rhoP::Float64,
    rhoQ::Float64,
    cut_edges,
)
    for zone_id in zone_ids
        init_ac_hybrid_partition!(
            get_partition(partitions, zone_id),
            sys.parts[zone_id],
            optimizer_factory,
            rho_theta,
            rhoV,
            rhoP,
            rhoQ,
            cut_edges,
        )
    end
    return nothing
end

# ============================================================
# RUNNERS (ADMM-AC)
# ============================================================

# ------------------------------------------------------------
# Runner: ADMM bus/AC
# ------------------------------------------------------------

function run_admm_ac(sys::ADMMSystem, config::DistributedOPFModel)::LirionRunResult
    t0_ns = time_ns()
    @assert config.base_model == AC "ADMM-AC: base_model deve ser AC."

    shared_buses = sys.shared_buses
    max_iter = config.max_iter

    rho_theta = Ref(config.rho_theta)
    rho_v = Ref(config.rho_v)
    prev_rho_theta = Ref(rho_theta[])
    prev_rho_v = Ref(rho_v[])

    eps_abs = config.eps_abs
    eps_rel = config.eps_rel
    adaptive = _adaptive_on(config)

    optimizer_factory = resolve_optimizer(config)

    partitions = build_ac_partitions(sys)
    zone_ids = partitions.zones
    n_zones = length(zone_ids)

    shared_buses_sorted, shared_buses_to_pos = sorted_ids_and_id2idx(shared_buses)
    n_shared_buses = length(shared_buses_sorted)

    duals_theta = Dict{Int, Vector{Float64}}(
        zone_id => zeros(Float64, n_shared_buses) for zone_id in zone_ids
    )
    duals_v = Dict{Int, Vector{Float64}}(
        zone_id => zeros(Float64, n_shared_buses) for zone_id in zone_ids
    )

    dual_theta_buf = Dict{Int, Vector{Float64}}(
        zone_id => zeros(Float64, n_shared_buses) for zone_id in zone_ids
    )
    dual_v_buf = Dict{Int, Vector{Float64}}(
        zone_id => zeros(Float64, n_shared_buses) for zone_id in zone_ids
    )

    theta_consensus = zeros(Float64, n_shared_buses)
    vmag_consensus = ones(Float64, n_shared_buses)
    theta_consensus_prev = similar(theta_consensus)
    vmag_consensus_prev = similar(vmag_consensus)

    hist = ADMMHistory()
    hist.meta[:kind] = :ac_bus
    hist.meta[:tag] = "admm_ac_bus"

    theta_avg_hist = Dict{Int, Vector{Float64}}()
    vmag_avg_hist = Dict{Int, Vector{Float64}}()

    results = Vector{NamedTuple{(:expr_cost, :real_cost), Tuple{Float64, Float64}}}(
        undef,
        n_zones,
    )

    init_ac_bus_partition!(
        partitions,
        sys,
        zone_ids,
        optimizer_factory,
        rho_theta[],
        rho_v[],
    )

    iter_done = 0
    for iter in 1:max_iter
        iter_done = iter
        _FAST_USE_LOGS && println("\n[ADMM-AC] Iteração $iter")

        refresh_rho_bus_ac!(partitions, rho_theta, prev_rho_theta, rho_v, prev_rho_v)

        alpha_used = nothing
        theta_target = theta_consensus
        vmag_target = vmag_consensus

        Threads.@threads for i in 1:n_zones
            zone_id = zone_ids[i]
            part = get_partition(partitions, zone_id)
            model = part.model::JuMP.Model

            lambda_theta = (dual_theta_buf[zone_id] .= duals_theta[zone_id])
            lambda_v = (dual_v_buf[zone_id] .= duals_v[zone_id])

            lin = _build_linear_bus_expr_ac!(
                part,
                shared_buses_sorted,
                shared_buses_to_pos,
                theta_target,
                vmag_target,
                lambda_theta,
                lambda_v,
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

        for (i, zone_id) in enumerate(zone_ids)
            part = get_partition(partitions, zone_id)
            out = results[i]
            part.has_solution = true
            part.generation_cost_model_expr = out.expr_cost
            part.real_generation_cost = out.real_cost
        end

        theta_consensus_prev .= theta_consensus
        vmag_consensus_prev .= vmag_consensus

        compute_ac_bus_consensus_vec!(
            theta_consensus,
            vmag_consensus,
            partitions,
            shared_buses_sorted,
            shared_buses_to_pos,
        )

        @inbounds for b in shared_buses_sorted
            j = shared_buses_to_pos[b]
            push!(get!(theta_avg_hist, b, Float64[]), theta_consensus[j])
            push!(get!(vmag_avg_hist, b, Float64[]), vmag_consensus[j])
        end

        update_duals_bus_ac_all!(
            duals_theta,
            duals_v,
            partitions,
            zone_ids,
            theta_consensus,
            vmag_consensus,
            rho_theta,
            rho_v,
            shared_buses_sorted,
            shared_buses_to_pos;
            alpha = alpha_used,
            theta_consensus_prev = theta_consensus_prev,
            vmag_consensus_prev = vmag_consensus_prev,
        )

        r = primal_residual_bus_ac_vec(
            partitions,
            theta_consensus,
            vmag_consensus,
            shared_buses_sorted,
            shared_buses_to_pos,
        )
        s = dual_residual_bus_ac_vec(
            theta_consensus,
            vmag_consensus,
            theta_consensus_prev,
            vmag_consensus_prev,
        )

        eps_r = eps_primal_bus_ac_vec(
            partitions,
            theta_consensus,
            vmag_consensus,
            shared_buses_sorted,
            shared_buses_to_pos,
            eps_abs,
            eps_rel,
        )
        eps_s = eps_dual_bus_ac_vec(
            theta_consensus,
            vmag_consensus,
            theta_consensus_prev,
            vmag_consensus_prev,
            eps_abs,
            eps_rel,
        )

        _record_iter!(
            hist, r, s, eps_r, eps_s;
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
            vmag_consensus,
            theta_consensus_prev,
            vmag_consensus_prev,
            shared_buses_sorted,
            shared_buses_to_pos,
            config,
        )
    end

    hist.heatmaps[:theta_avg] = theta_avg_hist
    hist.heatmaps[:vmag_avg] = vmag_avg_hist

    total_real_cost = sum(p.real_generation_cost for p in partitions.parts)
    total_expr_cost = sum(p.generation_cost_model_expr for p in partitions.parts)

    elapsed_s = (time_ns() - t0_ns) / 1e9
    println("\n[ADMM-AC] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("rho_theta inicial: $(config.rho_theta) | rho_theta final: $(rho_theta[])")
    println("rho_v     inicial: $(config.rho_v)     | rho_v     final: $(rho_v[])")
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (bus/AC): $(_pretty_time(elapsed_s))")

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
# Runner: ADMM branch/AC
# ------------------------------------------------------------

function run_admm_ac_branch(
    sys::ADMMSystem,
    config::DistributedOPFModel,
    warm_start::Bool = false,
)::LirionRunResult
    t0_ns = time_ns()
    @assert config.base_model == AC "ADMM-AC-branch: base_model deve ser AC."

    cut_edges = sys.cut_edges
    max_iter = config.max_iter

    rho_p = Ref(config.rho_flow_p)
    rho_q = Ref(config.rho_flow_q)
    prev_rho_p = Ref(rho_p[])
    prev_rho_q = Ref(rho_q[])

    eps_abs = config.eps_abs
    eps_rel = config.eps_rel
    adaptive = _adaptive_on(config)

    optimizer_factory = resolve_optimizer(config)

    partitions = build_ac_partitions(sys)
    zone_ids = partitions.zones
    n_zones = length(zone_ids)

    duals_p = Dict{Int, Dict{Int, Float64}}(
        zone_id => Dict{Int, Float64}() for zone_id in zone_ids
    )
    duals_q = Dict{Int, Dict{Int, Float64}}(
        zone_id => Dict{Int, Float64}() for zone_id in zone_ids
    )

    pflow_consensus = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    qflow_consensus = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    pflow_consensus_prev = Dict{Int, Float64}()
    qflow_consensus_prev = Dict{Int, Float64}()

    hist = ADMMHistory()
    hist.meta[:kind] = :ac_branch
    hist.meta[:tag] = "admm_ac_branch"

    results = Vector{NamedTuple{(:expr_cost, :real_cost), Tuple{Float64, Float64}}}(
        undef,
        n_zones,
    )

    init_ac_branch_partition!(
        partitions,
        sys,
        zone_ids,
        optimizer_factory,
        rho_p[],
        rho_q[],
        cut_edges,
    )

    iter_done = 0
    for iter in 1:max_iter
        iter_done = iter
        _FAST_USE_LOGS && println("\n[ADMM-AC-branch] Iteração $iter")

        refresh_rho_branch_ac!(partitions, rho_p, prev_rho_p, rho_q, prev_rho_q)

        alpha_used = nothing

        Threads.@threads for i in 1:n_zones
            zone_id = zone_ids[i]
            part = get_partition(partitions, zone_id)
            model = part.model::JuMP.Model

            lin = _build_linear_branch_expr_ac!(
                part,
                part.zone_edges,
                pflow_consensus,
                qflow_consensus,
                duals_p[zone_id],
                duals_q[zone_id],
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

        for (i, zone_id) in enumerate(zone_ids)
            part = get_partition(partitions, zone_id)
            out = results[i]
            part.has_solution = true
            part.generation_cost_model_expr = out.expr_cost
            part.real_generation_cost = out.real_cost
        end

        pflow_consensus_prev = copy(pflow_consensus)
        qflow_consensus_prev = copy(qflow_consensus)

        compute_ac_branch_consensus!(
            pflow_consensus,
            qflow_consensus,
            partitions,
            cut_edges,
        )

        update_duals_branch_ac!(
            duals_p,
            duals_q,
            partitions,
            pflow_consensus,
            qflow_consensus,
            cut_edges,
            rho_p[],
            rho_q[];
            alpha = alpha_used,
            pflow_consensus_prev = pflow_consensus_prev,
            qflow_consensus_prev = qflow_consensus_prev,
        )

        r = primal_residual_branch_ac(
            partitions,
            pflow_consensus,
            qflow_consensus,
            cut_edges,
        )
        s = dual_residual_branch_ac(
            pflow_consensus,
            qflow_consensus,
            pflow_consensus_prev,
            qflow_consensus_prev,
            cut_edges,
        )

        eps_r = eps_primal_branch_ac(
            partitions,
            pflow_consensus,
            qflow_consensus,
            cut_edges,
            eps_abs,
            eps_rel,
        )
        eps_s = eps_dual_branch_ac(
            pflow_consensus,
            qflow_consensus,
            pflow_consensus_prev,
            qflow_consensus_prev,
            cut_edges,
            eps_abs,
            eps_rel,
        )

        _record_iter!(
            hist, r, s, eps_r, eps_s;
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

    total_real_cost = sum(p.real_generation_cost for p in partitions.parts)
    total_expr_cost = sum(p.generation_cost_model_expr for p in partitions.parts)

    elapsed_s = (time_ns() - t0_ns) / 1e9
    println("\n[ADMM-AC-branch] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("rhoP inicial: $(config.rho_flow_p) | rhoP final: $(rho_p[])")
    println("rhoQ inicial: $(config.rho_flow_q) | rhoQ final: $(rho_q[])")
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (branch/AC): $(_pretty_time(elapsed_s))")

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
# Runner: ADMM hybrid/AC
# ------------------------------------------------------------

function run_admm_ac_hybrid(sys::ADMMSystem, config::DistributedOPFModel)::LirionRunResult
    t0_ns = time_ns()
    @assert config.base_model == AC "Hybrid AC: base_model deve ser AC."

    shared_buses = sys.shared_buses
    cut_edges = sys.cut_edges
    max_iter = config.max_iter

    rho_theta = Ref(config.rho_theta)
    rho_v = Ref(config.rho_v)
    rho_p = Ref(config.rho_flow_p)
    rho_q = Ref(config.rho_flow_q)

    prev_rho_theta = Ref(rho_theta[])
    prev_rho_v = Ref(rho_v[])
    prev_rho_p = Ref(rho_p[])
    prev_rho_q = Ref(rho_q[])

    eps_abs = config.eps_abs
    eps_rel = config.eps_rel
    adaptive = _adaptive_on(config)

    optimizer_factory = resolve_optimizer(config)

    partitions = build_ac_partitions(sys)
    zone_ids = partitions.zones
    n_zones = length(zone_ids)

    shared_buses_sorted, shared_buses_to_pos = sorted_ids_and_id2idx(shared_buses)
    n_shared_buses = length(shared_buses_sorted)

    duals_theta = Dict{Int, Vector{Float64}}(
        zone_id => zeros(Float64, n_shared_buses) for zone_id in zone_ids
    )
    duals_v = Dict{Int, Vector{Float64}}(
        zone_id => zeros(Float64, n_shared_buses) for zone_id in zone_ids
    )

    dual_theta_buf = Dict{Int, Vector{Float64}}(
        zone_id => zeros(Float64, n_shared_buses) for zone_id in zone_ids
    )
    dual_v_buf = Dict{Int, Vector{Float64}}(
        zone_id => zeros(Float64, n_shared_buses) for zone_id in zone_ids
    )

    duals_p = Dict{Int, Dict{Int, Float64}}(
        zone_id => Dict{Int, Float64}() for zone_id in zone_ids
    )
    duals_q = Dict{Int, Dict{Int, Float64}}(
        zone_id => Dict{Int, Float64}() for zone_id in zone_ids
    )

    theta_consensus = zeros(Float64, n_shared_buses)
    vmag_consensus = ones(Float64, n_shared_buses)
    theta_consensus_prev = similar(theta_consensus)
    vmag_consensus_prev = similar(vmag_consensus)

    pflow_consensus = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    qflow_consensus = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    pflow_consensus_prev = Dict{Int, Float64}()
    qflow_consensus_prev = Dict{Int, Float64}()

    hist = ADMMHistory()
    hist.meta[:kind] = :ac_hybrid
    hist.meta[:tag] = "admm_ac_hybrid"

    theta_avg_hist = Dict{Int, Vector{Float64}}()
    vmag_avg_hist = Dict{Int, Vector{Float64}}()

    results = Vector{NamedTuple{(:expr_cost, :real_cost), Tuple{Float64, Float64}}}(
        undef,
        n_zones,
    )

    init_ac_hybrid_partition!(
        partitions, sys, zone_ids, optimizer_factory,
        rho_theta[], rho_v[], rho_p[], rho_q[], cut_edges,
    )

    iter_done = 0
    for iter in 1:max_iter
        iter_done = iter
        _FAST_USE_LOGS && println("\n[ADMM-AC-hybrid] Iteração $iter")

        refresh_rho_bus_ac!(partitions, rho_theta, prev_rho_theta, rho_v, prev_rho_v)
        refresh_rho_branch_ac!(partitions, rho_p, prev_rho_p, rho_q, prev_rho_q)

        alpha_used = nothing
        theta_target = theta_consensus
        vmag_target = vmag_consensus

        Threads.@threads for i in 1:n_zones
            zone_id = zone_ids[i]
            part = get_partition(partitions, zone_id)
            model = part.model::JuMP.Model

            lambda_theta = (dual_theta_buf[zone_id] .= duals_theta[zone_id])
            lambda_v = (dual_v_buf[zone_id] .= duals_v[zone_id])

            lin_bus = _build_linear_bus_expr_ac!(
                part,
                shared_buses_sorted,
                shared_buses_to_pos,
                theta_target,
                vmag_target,
                lambda_theta,
                lambda_v,
                rho_theta[],
                rho_v[],
            )
            lin_br = _build_linear_branch_expr_ac!(
                part,
                part.zone_edges,
                pflow_consensus,
                qflow_consensus,
                duals_p[zone_id],
                duals_q[zone_id],
                rho_p[],
                rho_q[],
            )

            @objective(
                model,
                Min,
                part.obj_base + part.quad_bus_expr_ac + part.quad_branch_expr_ac +
                (lin_bus + lin_br)
            )
            optimize!(model)

            results[i] = (
                expr_cost = value(model[:generation_cost]),
                real_cost = compute_real_generation_cost(model, part.objective_params),
            )
        end

        for (i, zone_id) in enumerate(zone_ids)
            part = get_partition(partitions, zone_id)
            out = results[i]
            part.has_solution = true
            part.generation_cost_model_expr = out.expr_cost
            part.real_generation_cost = out.real_cost
        end

        theta_consensus_prev .= theta_consensus
        vmag_consensus_prev .= vmag_consensus
        pflow_consensus_prev = copy(pflow_consensus)
        qflow_consensus_prev = copy(qflow_consensus)

        compute_ac_bus_consensus_vec!(
            theta_consensus,
            vmag_consensus,
            partitions,
            shared_buses_sorted,
            shared_buses_to_pos,
        )

        @inbounds for b in shared_buses_sorted
            j = shared_buses_to_pos[b]
            push!(get!(theta_avg_hist, b, Float64[]), theta_consensus[j])
            push!(get!(vmag_avg_hist, b, Float64[]), vmag_consensus[j])
        end

        compute_ac_branch_consensus!(
            pflow_consensus,
            qflow_consensus,
            partitions,
            cut_edges,
        )

        update_duals_bus_ac_all!(
            duals_theta,
            duals_v,
            partitions,
            zone_ids,
            theta_consensus,
            vmag_consensus,
            rho_theta,
            rho_v,
            shared_buses_sorted,
            shared_buses_to_pos;
            alpha = alpha_used,
            theta_consensus_prev = theta_consensus_prev,
            vmag_consensus_prev = vmag_consensus_prev,
        )

        update_duals_branch_ac!(
            duals_p,
            duals_q,
            partitions,
            pflow_consensus,
            qflow_consensus,
            cut_edges,
            rho_p[],
            rho_q[];
            alpha = alpha_used,
            pflow_consensus_prev = pflow_consensus_prev,
            qflow_consensus_prev = qflow_consensus_prev,
        )

        r_bus = primal_residual_bus_ac_vec(
            partitions,
            theta_consensus,
            vmag_consensus,
            shared_buses_sorted,
            shared_buses_to_pos,
        )
        r_br = primal_residual_branch_ac(
            partitions,
            pflow_consensus,
            qflow_consensus,
            cut_edges,
        )
        r_tot = hypot(r_bus, r_br)

        s_bus = dual_residual_bus_ac_vec(
            theta_consensus,
            vmag_consensus,
            theta_consensus_prev,
            vmag_consensus_prev,
        )
        s_br = dual_residual_branch_ac(
            pflow_consensus,
            qflow_consensus,
            pflow_consensus_prev,
            qflow_consensus_prev,
            cut_edges,
        )
        s_tot = hypot(s_bus, s_br)

        eps_bus_r = eps_primal_bus_ac_vec(
            partitions,
            theta_consensus,
            vmag_consensus,
            shared_buses_sorted,
            shared_buses_to_pos,
            eps_abs,
            eps_rel,
        )
        eps_br_r = eps_primal_branch_ac(
            partitions,
            pflow_consensus,
            qflow_consensus,
            cut_edges,
            eps_abs,
            eps_rel,
        )
        eps_r_tot = hypot(eps_bus_r, eps_br_r)

        eps_bus_s = eps_dual_bus_ac_vec(
            theta_consensus,
            vmag_consensus,
            theta_consensus_prev,
            vmag_consensus_prev,
            eps_abs,
            eps_rel,
        )
        eps_br_s = eps_dual_branch_ac(
            pflow_consensus,
            qflow_consensus,
            pflow_consensus_prev,
            qflow_consensus_prev,
            cut_edges,
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
            vmag_consensus,
            theta_consensus_prev,
            vmag_consensus_prev,
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
    hist.heatmaps[:vmag_avg] = vmag_avg_hist

    total_real_cost = sum(p.real_generation_cost for p in partitions.parts)
    total_expr_cost = sum(p.generation_cost_model_expr for p in partitions.parts)

    elapsed_s = (time_ns() - t0_ns) / 1e9
    println("\n[ADMM-AC-hybrid] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("rho_theta inicial: $(config.rho_theta) | rho_theta final: $(rho_theta[])")
    println("rho_v     inicial: $(config.rho_v)     | rho_v     final: $(rho_v[])")
    println("rhoP      inicial: $(config.rho_flow_p)| rhoP      final: $(rho_p[])")
    println("rhoQ      inicial: $(config.rho_flow_q)| rhoQ      final: $(rho_q[])")
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (hybrid/AC): $(_pretty_time(elapsed_s))")

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
