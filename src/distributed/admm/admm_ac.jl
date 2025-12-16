function init_ac_bus_partition!(
    part::PartitionAC,
    sys_zone,
    optimizer_factory,
    rho_theta::Float64,
    rhoV::Float64,
)
    prepare_hotpaths_ac!(part)

    m = build_model(sys_zone, AC; optimizer = optimizer_factory)
    part.model = m
    part.obj_base = m[:generation_cost]
    part.theta_var = m[:theta]
    part.vmag_var = m[:v_magnitude]
    part.objective_params = sys_zone.objective
    part.has_solution = false

    _set_model_silent!(m)
    for I in keys(part.theta_var)
        set_start_value(part.theta_var[I], 0.0)
    end
    for I in keys(part.vmag_var)
        set_start_value(part.vmag_var[I], 1.0)
    end

    _build_quad_bus_expr_ac!(part, rho_theta, rhoV)
    @objective(m, Min, part.obj_base + part.quad_bus_expr_ac)

    return nothing
end

function init_ac_bus_partition!(
    partitions::Dict{Int, PartitionAC},
    sys::ADMMSystem,
    zone_ids::Vector{Int},
    optimizer_factory,
    rho_theta::Float64,
    rhoV::Float64,
)
    for z in zone_ids
        init_ac_bus_partition!(
            partitions[z],
            sys.parts[z],
            optimizer_factory,
            rho_theta,
            rhoV,
        )
    end
    return nothing
end

function init_ac_branch_partition!(
    part::PartitionAC,
    sys_zone,
    optimizer_factory,
    rhoP::Float64,
    rhoQ::Float64,
    cut_edges,
)
    prepare_hotpaths_ac!(part)
    part.zone_edges = _zone_cut_edges(part, cut_edges)

    m = build_model(sys_zone, AC; optimizer = optimizer_factory)
    part.model = m
    part.obj_base = m[:generation_cost]
    part.theta_var = m[:theta]
    part.vmag_var = m[:v_magnitude]
    part.objective_params = sys_zone.objective
    part.has_solution = false

    _set_model_silent!(m)
    for I in keys(part.theta_var)
        set_start_value(part.theta_var[I], 0.0)
    end
    for I in keys(part.vmag_var)
        set_start_value(part.vmag_var[I], 1.0)
    end

    _build_quad_branch_expr_ac!(part, rhoP, rhoQ, part.zone_edges)
    @objective(m, Min, part.obj_base + part.quad_branch_expr_ac)

    return nothing
end

# BRANCH/AC init — todas as partições
function init_ac_branch_partition!(
    partitions::Dict{Int, PartitionAC},
    sys::ADMMSystem,
    zone_ids::Vector{Int},
    optimizer_factory,
    rhoP::Float64,
    rhoQ::Float64,
    cut_edges,
)
    for z in zone_ids
        init_ac_branch_partition!(
            partitions[z],
            sys.parts[z],
            optimizer_factory,
            rhoP,
            rhoQ,
            cut_edges,
        )
    end
    return nothing
end

# ------------------------------------------------------------
# HYBRID/AC init (bus theta,V + branch P,Q) — uma partição
# ------------------------------------------------------------
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
    prepare_hotpaths_ac!(part)
    part.zone_edges = _zone_cut_edges(part, cut_edges)

    m = build_model(sys_zone, AC; optimizer = optimizer_factory)
    part.model = m
    part.obj_base = m[:generation_cost]
    part.theta_var = m[:theta]
    part.vmag_var = m[:v_magnitude]
    part.objective_params = sys_zone.objective
    part.has_solution = false

    _set_model_silent!(m)
    for I in keys(part.theta_var)
        set_start_value(part.theta_var[I], 0.0)
    end
    for I in keys(part.vmag_var)
        set_start_value(part.vmag_var[I], 1.0)
    end

    _build_quad_bus_expr_ac!(part, rho_theta, rhoV)
    _build_quad_branch_expr_ac!(part, rhoP, rhoQ, part.zone_edges)

    @objective(m, Min, part.obj_base + part.quad_bus_expr_ac + part.quad_branch_expr_ac)

    return nothing
end

# HYBRID/AC init — todas as partições
function init_ac_hybrid_partition!(
    partitions::Dict{Int, PartitionAC},
    sys::ADMMSystem,
    zone_ids::Vector{Int},
    optimizer_factory,
    rho_theta::Float64,
    rhoV::Float64,
    rhoP::Float64,
    rhoQ::Float64,
    cut_edges,
)
    for z in zone_ids
        init_ac_hybrid_partition!(
            partitions[z],
            sys.parts[z],
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
# RUNNERS — removendo os "for z in zone_ids init_..."
# ============================================================

function run_admm_ac(config::DistributedOPFModel, sys::ADMMSystem)
    t0_ns = time_ns()

    @assert config.base_model == AC "ADMM-AC: base_model deve ser AC."

    shared_buses = sys.shared_buses
    max_iter = config.max_iter

    rho_theta = Ref(config.rho_theta)
    rho_v = Ref(config.rho_v)

    eps_abs = config.eps_primal
    eps_rel = config.eps_dual
    adaptive = _adaptive_on(config)

    optimizer_factory = resolve_optimizer_ac(config)

    partitions = build_ac_partitions(sys)
    zone_ids = collect(keys(partitions))
    n_zones = length(zone_ids)

    shared_buses_sorted, shared_buses_idx = sorted_ids_and_pos(shared_buses)
    n_shared_buses = length(shared_buses_sorted)

    duals_theta =
        Dict{Int, Vector{Float64}}(z => zeros(Float64, n_shared_buses) for z in zone_ids)
    duals_v =
        Dict{Int, Vector{Float64}}(z => zeros(Float64, n_shared_buses) for z in zone_ids)

    lambdatheta_zone_buf =
        Dict{Int, Vector{Float64}}(z => zeros(Float64, n_shared_buses) for z in zone_ids)
    lambdaV_zone_buf =
        Dict{Int, Vector{Float64}}(z => zeros(Float64, n_shared_buses) for z in zone_ids)

    theta_consensus = zeros(Float64, n_shared_buses)
    vmag_consensus = ones(Float64, n_shared_buses)

    theta_consensus_prev = similar(theta_consensus)
    vmag_consensus_prev = similar(vmag_consensus)

    res_primal_hist = Float64[]
    res_dual_hist = Float64[]
    eps_primal_hist = Float64[]
    eps_dual_hist = Float64[]
    rho_theta_hist = Float64[]
    rho_v_hist = Float64[]
    theta_avg_hist = Dict{Int, Vector{Float64}}()

    prev_rho_theta = rho_theta[]
    prev_rho_v = rho_v[]

    results = Vector{NamedTuple{(:expr_cost, :real_cost),
        Tuple{Float64, Float64}}}(
        undef,
        n_zones,
    )

    # ✅ init (some o for do runner)
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

        if (rho_theta[] != prev_rho_theta) || (rho_v[] != prev_rho_v)
            for part in values(partitions)
                _build_quad_bus_expr_ac!(part, rho_theta[], rho_v[])
            end
            prev_rho_theta = rho_theta[]
            prev_rho_v = rho_v[]
        end

        α_used = nothing
        theta_target = theta_consensus
        vmag_target = vmag_consensus

        Threads.@threads for i in 1:n_zones
            z = zone_ids[i]
            part = partitions[z]
            m = part.model::JuMP.Model

            lambdatheta = (lambdatheta_zone_buf[z] .= duals_theta[z])
            lambdaV = (lambdaV_zone_buf[z] .= duals_v[z])

            lin = _build_linear_bus_expr_ac!(
                part,
                shared_buses_sorted,
                shared_buses_idx,
                theta_target,
                vmag_target,
                lambdatheta,
                lambdaV,
                rho_theta[],
                rho_v[],
            )

            @objective(m, Min, part.obj_base + part.quad_bus_expr_ac + lin)
            optimize!(m)

            results[i] = (
                expr_cost = value(m[:generation_cost]),
                real_cost = compute_real_generation_cost(m, part.objective_params),
            )
        end

        for (i, z) in enumerate(zone_ids)
            part = partitions[z]
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
            shared_buses_idx,
        )

        @inbounds for b in shared_buses_sorted
            j = shared_buses_idx[b]
            push!(get!(theta_avg_hist, b, Float64[]), theta_consensus[j])
        end

        for z in zone_ids
            update_duals_bus_ac_vec!(
                duals_theta[z],
                duals_v[z],
                partitions[z],
                theta_consensus,
                vmag_consensus,
                rho_theta[],
                rho_v[],
                shared_buses_sorted,
                shared_buses_idx;
                α = α_used,
                theta_consensus_prev = theta_consensus_prev,
                vmag_consensus_prev = vmag_consensus_prev,
            )
        end

        r = primal_residual_bus_ac_vec(
            partitions,
            theta_consensus,
            vmag_consensus,
            shared_buses_sorted,
            shared_buses_idx,
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
            shared_buses_idx,
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

        push!(res_primal_hist, r)
        push!(res_dual_hist, s)
        push!(eps_primal_hist, eps_r)
        push!(eps_dual_hist, eps_s)
        push!(rho_theta_hist, rho_theta[])
        push!(rho_v_hist, rho_v[])

        _FAST_USE_LOGS && println("‖r‖ = $r | eps_pri = $eps_r")
        _FAST_USE_LOGS && println("‖s‖ = $s | eps_du  = $eps_s")

        if r < eps_r && s < eps_s
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        if adaptive
            rtheta, rV = primal_residual_bus_ac_components(
                partitions,
                theta_consensus,
                vmag_consensus,
                shared_buses_sorted,
                shared_buses_idx,
            )
            stheta, sV = dual_residual_bus_ac_components(
                theta_consensus,
                vmag_consensus,
                theta_consensus_prev,
                vmag_consensus_prev,
            )

            adaptive_rho_ac!(
                rho_theta,
                rtheta,
                stheta;
                mu_incr = config.mu_incr,
                mu_decr = config.mu_decr,
                tau_incr = config.tau_incr,
                tau_decr = config.tau_decr,
            )
            adaptive_rho_ac!(
                rho_v,
                rV,
                sV;
                mu_incr = config.mu_incr,
                mu_decr = config.mu_decr,
                tau_incr = config.tau_incr,
                tau_decr = config.tau_decr,
            )
        end
    end

    save_ac_bus_plots(
        res_primal_hist, res_dual_hist, eps_primal_hist, eps_dual_hist,
        rho_theta_hist, rho_v_hist, theta_avg_hist;
        stem = "admm_ac",
    )

    total_real_cost = sum(p.real_generation_cost for p in values(partitions))
    total_expr_cost = sum(p.generation_cost_model_expr for p in values(partitions))

    elapsed_s = (time_ns() - t0_ns) / 1e9
    println("\n[ADMM-AC] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("rho_theta inicial: $(config.rho_theta) | rho_theta final: $(rho_theta[])")
    println("rho_v     inicial: $(config.rho_v)     | rho_v     final: $(rho_v[])")
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (bus/AC): $(_pretty_time(elapsed_s))")

    return total_real_cost, partitions, iter_done, elapsed_s
end

function run_admm_ac_branch(
    config::DistributedOPFModel,
    sys::ADMMSystem;
    warm_start::Bool = false,
)
    t0_ns = time_ns()
    @assert config.base_model == AC "ADMM-AC-branch: base_model deve ser AC."

    cut_edges = sys.cut_edges
    max_iter = config.max_iter

    rho_p = Ref(config.rho_flow_p)
    rho_q = Ref(config.rho_flow_q)

    eps_abs = config.eps_primal
    eps_rel = config.eps_dual
    adaptive = _adaptive_on(config)

    optimizer_factory = resolve_optimizer_ac(config)

    partitions = build_ac_partitions(sys)
    zone_ids = collect(keys(partitions))
    n_zones = length(zone_ids)

    duals_p = Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in zone_ids)
    duals_q = Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in zone_ids)

    zP = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zQ = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zP_prev = Dict{Int, Float64}()
    zQ_prev = Dict{Int, Float64}()

    res_primal_hist = Float64[]
    res_dual_hist = Float64[]
    eps_primal_hist = Float64[]
    eps_dual_hist = Float64[]
    rho_p_hist = Float64[]
    rho_q_hist = Float64[]

    prev_rho_p = rho_p[]
    prev_rho_q = rho_q[]

    results = Vector{NamedTuple{(:expr_cost, :real_cost),
        Tuple{Float64, Float64}}}(
        undef,
        n_zones,
    )

    # ✅ init (some o for do runner)
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

        if (rho_p[] != prev_rho_p) || (rho_q[] != prev_rho_q)
            for part in values(partitions)
                _build_quad_branch_expr_ac!(part, rho_p[], rho_q[], part.zone_edges)
            end
            prev_rho_p = rho_p[]
            prev_rho_q = rho_q[]
        end

        α_used = nothing

        Threads.@threads for i in 1:n_zones
            z = zone_ids[i]
            part = partitions[z]
            m = part.model::JuMP.Model

            lin = _build_linear_branch_expr_ac!(
                part,
                part.zone_edges,
                zP,
                zQ,
                duals_p[z],
                duals_q[z],
                rho_p[],
                rho_q[],
            )

            @objective(m, Min, part.obj_base + part.quad_branch_expr_ac + lin)
            optimize!(m)

            results[i] = (
                expr_cost = value(m[:generation_cost]),
                real_cost = compute_real_generation_cost(m, part.objective_params),
            )
        end

        for (i, z) in enumerate(zone_ids)
            part = partitions[z]
            out = results[i]
            part.has_solution = true
            part.generation_cost_model_expr = out.expr_cost
            part.real_generation_cost = out.real_cost
        end

        zP_prev = copy(zP)
        zQ_prev = copy(zQ)
        compute_ac_branch_consensus!(zP, zQ, partitions, cut_edges)

        update_duals_branch_ac!(
            duals_p,
            duals_q,
            partitions,
            zP,
            zQ,
            cut_edges,
            rho_p[],
            rho_q[];
            α = α_used,
            zP_prev = zP_prev,
            zQ_prev = zQ_prev,
        )

        r = primal_residual_branch_ac(partitions, zP, zQ, cut_edges)
        s = dual_residual_branch_ac(zP, zQ, zP_prev, zQ_prev, cut_edges)

        eps_r = eps_primal_branch_ac(partitions, zP, zQ, cut_edges, eps_abs, eps_rel)
        eps_s = eps_dual_branch_ac(zP, zQ, zP_prev, zQ_prev, cut_edges, eps_abs, eps_rel)

        push!(res_primal_hist, r)
        push!(res_dual_hist, s)
        push!(eps_primal_hist, eps_r)
        push!(eps_dual_hist, eps_s)
        push!(rho_p_hist, rho_p[])
        push!(rho_q_hist, rho_q[])

        _FAST_USE_LOGS && println("‖r‖ = $r | eps_pri = $eps_r")
        _FAST_USE_LOGS && println("‖s‖ = $s | eps_du  = $eps_s")

        if r < eps_r && s < eps_s
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        if adaptive
            adaptive_rho_ac!(
                rho_p,
                r,
                s;
                mu_incr = config.mu_incr,
                mu_decr = config.mu_decr,
                tau_incr = config.tau_incr,
                tau_decr = config.tau_decr,
            )
            adaptive_rho_ac!(
                rho_q,
                r,
                s;
                mu_incr = config.mu_incr,
                mu_decr = config.mu_decr,
                tau_incr = config.tau_incr,
                tau_decr = config.tau_decr,
            )
        end
    end

    save_ac_branch_plots(
        res_primal_hist, res_dual_hist, eps_primal_hist, eps_dual_hist,
        rho_p_hist, rho_q_hist;
        stem = "admm_ac_branch",
    )

    total_real_cost = sum(p.real_generation_cost for p in values(partitions))
    total_expr_cost = sum(p.generation_cost_model_expr for p in values(partitions))

    elapsed_s = (time_ns() - t0_ns) / 1e9
    println("\n[ADMM-AC-branch] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("rhoP inicial: $(config.rho_flow_p) | rhoP final: $(rho_p[])")
    println("rhoQ inicial: $(config.rho_flow_q) | rhoQ final: $(rho_q[])")
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (branch/AC): $(_pretty_time(elapsed_s))")

    return total_real_cost, partitions, iter_done, elapsed_s
end

function run_admm_ac_hybrid(config::DistributedOPFModel, sys::ADMMSystem)
    t0_ns = time_ns()
    @assert config.base_model == AC "Hybrid AC: base_model deve ser AC."

    shared_buses = sys.shared_buses
    cut_edges = sys.cut_edges
    max_iter = config.max_iter

    rho_theta = Ref(config.rho_theta)
    rho_v = Ref(config.rho_v)
    rho_p = Ref(config.rho_flow_p)
    rho_q = Ref(config.rho_flow_q)

    eps_abs = config.eps_primal
    eps_rel = config.eps_dual
    adaptive = _adaptive_on(config)

    optimizer_factory = resolve_optimizer_ac(config)

    partitions = build_ac_partitions(sys)
    zone_ids = collect(keys(partitions))
    n_zones = length(zone_ids)

    shared_buses_sorted, shared_buses_idx = sorted_ids_and_pos(shared_buses)
    n_shared_buses = length(shared_buses_sorted)

    duals_theta =
        Dict{Int, Vector{Float64}}(z => zeros(Float64, n_shared_buses) for z in zone_ids)
    duals_v =
        Dict{Int, Vector{Float64}}(z => zeros(Float64, n_shared_buses) for z in zone_ids)

    lambdatheta_zone_buf =
        Dict{Int, Vector{Float64}}(z => zeros(Float64, n_shared_buses) for z in zone_ids)
    lambdaV_zone_buf =
        Dict{Int, Vector{Float64}}(z => zeros(Float64, n_shared_buses) for z in zone_ids)

    duals_p = Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in zone_ids)
    duals_q = Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in zone_ids)

    theta_consensus = zeros(Float64, n_shared_buses)
    vmag_consensus = ones(Float64, n_shared_buses)
    theta_consensus_prev = similar(theta_consensus)
    vmag_consensus_prev = similar(vmag_consensus)

    zP = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zQ = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zP_prev = Dict{Int, Float64}()
    zQ_prev = Dict{Int, Float64}()

    res_primal_hist = Float64[]
    res_dual_hist = Float64[]
    eps_primal_hist = Float64[]
    eps_dual_hist = Float64[]

    rho_theta_hist = Float64[]
    rho_v_hist = Float64[]
    rho_p_hist = Float64[]
    rho_q_hist = Float64[]

    theta_avg_hist = Dict{Int, Vector{Float64}}()

    prev_rho_theta = rho_theta[]
    prev_rho_v = rho_v[]
    prev_rho_p = rho_p[]
    prev_rho_q = rho_q[]

    results = Vector{NamedTuple{(:expr_cost, :real_cost),
        Tuple{Float64, Float64}}}(
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

        if (rho_theta[] != prev_rho_theta) || (rho_v[] != prev_rho_v)
            for part in values(partitions)
                _build_quad_bus_expr_ac!(part, rho_theta[], rho_v[])
            end
            prev_rho_theta = rho_theta[]
            prev_rho_v = rho_v[]
        end

        if (rho_p[] != prev_rho_p) || (rho_q[] != prev_rho_q)
            for part in values(partitions)
                _build_quad_branch_expr_ac!(part, rho_p[], rho_q[], part.zone_edges)
            end
            prev_rho_p = rho_p[]
            prev_rho_q = rho_q[]
        end

        α_used = nothing
        theta_target = theta_consensus
        vmag_target = vmag_consensus

        Threads.@threads for i in 1:n_zones
            z = zone_ids[i]
            part = partitions[z]
            m = part.model::JuMP.Model

            lambdatheta = (lambdatheta_zone_buf[z] .= duals_theta[z])
            lambdaV = (lambdaV_zone_buf[z] .= duals_v[z])

            lin_bus = _build_linear_bus_expr_ac!(
                part,
                shared_buses_sorted,
                shared_buses_idx,
                theta_target,
                vmag_target,
                lambdatheta,
                lambdaV,
                rho_theta[],
                rho_v[],
            )
            lin_br = _build_linear_branch_expr_ac!(
                part,
                part.zone_edges,
                zP,
                zQ,
                duals_p[z],
                duals_q[z],
                rho_p[],
                rho_q[],
            )

            @objective(
                m,
                Min,
                part.obj_base + part.quad_bus_expr_ac + part.quad_branch_expr_ac +
                (lin_bus + lin_br)
            )
            optimize!(m)

            results[i] = (
                expr_cost = value(m[:generation_cost]),
                real_cost = compute_real_generation_cost(m, part.objective_params),
            )
        end

        for (i, z) in enumerate(zone_ids)
            part = partitions[z]
            out = results[i]
            part.has_solution = true
            part.generation_cost_model_expr = out.expr_cost
            part.real_generation_cost = out.real_cost
        end

        theta_consensus_prev .= theta_consensus
        vmag_consensus_prev .= vmag_consensus
        zP_prev = copy(zP)
        zQ_prev = copy(zQ)

        compute_ac_bus_consensus_vec!(
            theta_consensus,
            vmag_consensus,
            partitions,
            shared_buses_sorted,
            shared_buses_idx,
        )
        @inbounds for b in shared_buses_sorted
            j = shared_buses_idx[b]
            push!(get!(theta_avg_hist, b, Float64[]), theta_consensus[j])
        end
        compute_ac_branch_consensus!(zP, zQ, partitions, cut_edges)

        for z in zone_ids
            update_duals_bus_ac_vec!(
                duals_theta[z],
                duals_v[z],
                partitions[z],
                theta_consensus,
                vmag_consensus,
                rho_theta[],
                rho_v[],
                shared_buses_sorted,
                shared_buses_idx;
                α = α_used,
                theta_consensus_prev = theta_consensus_prev,
                vmag_consensus_prev = vmag_consensus_prev,
            )
        end
        update_duals_branch_ac!(
            duals_p,
            duals_q,
            partitions,
            zP,
            zQ,
            cut_edges,
            rho_p[],
            rho_q[];
            α = α_used,
            zP_prev = zP_prev,
            zQ_prev = zQ_prev,
        )

        r_bus = primal_residual_bus_ac_vec(
            partitions,
            theta_consensus,
            vmag_consensus,
            shared_buses_sorted,
            shared_buses_idx,
        )
        r_br = primal_residual_branch_ac(partitions, zP, zQ, cut_edges)
        r_tot = hypot(r_bus, r_br)

        s_bus = dual_residual_bus_ac_vec(
            theta_consensus,
            vmag_consensus,
            theta_consensus_prev,
            vmag_consensus_prev,
        )
        s_br = dual_residual_branch_ac(zP, zQ, zP_prev, zQ_prev, cut_edges)
        s_tot = hypot(s_bus, s_br)

        eps_bus_r = eps_primal_bus_ac_vec(
            partitions,
            theta_consensus,
            vmag_consensus,
            shared_buses_sorted,
            shared_buses_idx,
            eps_abs,
            eps_rel,
        )
        eps_br_r = eps_primal_branch_ac(partitions, zP, zQ, cut_edges, eps_abs, eps_rel)
        eps_r_tot = hypot(eps_bus_r, eps_br_r)

        eps_bus_s = eps_dual_bus_ac_vec(
            theta_consensus,
            vmag_consensus,
            theta_consensus_prev,
            vmag_consensus_prev,
            eps_abs,
            eps_rel,
        )
        eps_br_s = eps_dual_branch_ac(zP, zQ, zP_prev, zQ_prev, cut_edges, eps_abs, eps_rel)
        eps_s_tot = hypot(eps_bus_s, eps_br_s)

        push!(res_primal_hist, r_tot)
        push!(res_dual_hist, s_tot)
        push!(eps_primal_hist, eps_r_tot)
        push!(eps_dual_hist, eps_s_tot)

        push!(rho_theta_hist, rho_theta[])
        push!(rho_v_hist, rho_v[])
        push!(rho_p_hist, rho_p[])
        push!(rho_q_hist, rho_q[])

        _FAST_USE_LOGS && println("‖r‖ = $r_tot | eps_pri = $eps_r_tot")
        _FAST_USE_LOGS && println("‖s‖ = $s_tot | eps_du  = $eps_s_tot")

        if r_tot < eps_r_tot && s_tot < eps_s_tot
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        if adaptive
            rtheta, rV = primal_residual_bus_ac_components(
                partitions,
                theta_consensus,
                vmag_consensus,
                shared_buses_sorted,
                shared_buses_idx,
            )
            stheta, sV = dual_residual_bus_ac_components(
                theta_consensus,
                vmag_consensus,
                theta_consensus_prev,
                vmag_consensus_prev,
            )

            adaptive_rho_ac!(
                rho_theta,
                rtheta,
                stheta;
                mu_incr = config.mu_incr,
                mu_decr = config.mu_decr,
                tau_incr = config.tau_incr,
                tau_decr = config.tau_decr,
            )
            adaptive_rho_ac!(
                rho_v,
                rV,
                sV;
                mu_incr = config.mu_incr,
                mu_decr = config.mu_decr,
                tau_incr = config.tau_incr,
                tau_decr = config.tau_decr,
            )

            adaptive_rho_ac!(
                rho_p,
                r_br,
                s_br;
                mu_incr = config.mu_incr,
                mu_decr = config.mu_decr,
                tau_incr = config.tau_incr,
                tau_decr = config.tau_decr,
            )
            adaptive_rho_ac!(
                rho_q,
                r_br,
                s_br;
                mu_incr = config.mu_incr,
                mu_decr = config.mu_decr,
                tau_incr = config.tau_incr,
                tau_decr = config.tau_decr,
            )
        end
    end

    total_real_cost = sum(p.real_generation_cost for p in values(partitions))
    total_expr_cost = sum(p.generation_cost_model_expr for p in values(partitions))

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

    return total_real_cost, partitions, iter_done, elapsed_s
end
