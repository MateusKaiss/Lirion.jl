############################
# Helpers: init partitions #
############################

# BUS/AC init (theta + V mag)
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

# BRANCH/AC init (P,Q on cut edges)
function init_ac_branch_partition!(
    part::PartitionAC,
    sys_zone,
    optimizer_factory,
    rhoP::Float64,
    rhoQ::Float64,
    cut_edges,
)
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

# HYBRID/AC init (bus theta,V + branch P,Q)
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


############################
# ADMM: BUS/AC              #
############################

function run_admm_ac(config::DistributedOPFModel, sys::ADMMSystem)
    t0 = time_ns()
    @assert config.base_model == AC "ADMM-AC: base_model deve ser AC."

    partitions = build_ac_partitions(sys)
    shared_buses = sys.shared_buses

    rho_theta = Ref(config.rho_theta)
    rhoV = Ref(config.rho_v)
    eps_abs, eps_rel = config.eps_primal, config.eps_dual
    adaptive = _adaptive_on(config)

    shared_buses_sorted, shared_buses_idx = prepare_dense_maps(shared_buses)

    duals_theta_dense = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )
    dualsV_dense = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )

    theta_consensus_vec = zeros(Float64, length(shared_buses_sorted))
    vmag_consensus_vec = ones(Float64, length(shared_buses_sorted))
    theta_consensus_prev = similar(theta_consensus_vec)
    vmag_consensus_prev = similar(vmag_consensus_vec)

    res_pri_hist = Float64[]
    res_du_hist = Float64[]
    eps_pri_hist = Float64[]
    eps_du_hist = Float64[]
    rho_theta_hist = Float64[]
    rhoV_hist = Float64[]
    theta_avg_hist = Dict{Int, Vector{Float64}}()

    optimizer_factory = resolve_optimizer_ac(config)

    # init por zona (hotpaths + model + quad bus)
    for (zone_id, part) in partitions
        init_ac_bus_partition!(
            part,
            sys.parts[zone_id],
            optimizer_factory,
            rho_theta[],
            rhoV[],
        )
    end

    prev_rho_theta, prev_rhoV = rho_theta[], rhoV[]
    zone_ids = collect(keys(partitions))
    lambdatheta_zone_buf = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in zone_ids
    )
    lambdaV_zone_buf = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in zone_ids
    )

    iter_done = 0

    for k in 1:config.max_iter
        iter_done = k
        _FAST_USE_LOGS && println("\n[ADMM-AC] Iteração $k")

        if rho_theta[] != prev_rho_theta || rhoV[] != prev_rhoV
            for (_, part) in partitions
                _build_quad_bus_expr_ac!(part, rho_theta[], rhoV[])
            end
            prev_rho_theta, prev_rhoV = rho_theta[], rhoV[]
        end

        α_used = nothing
        theta_consensus_target = theta_consensus_vec
        vmag_consensus_target = vmag_consensus_vec

        nz = length(zone_ids)
        results = Vector{Any}(undef, nz)

        Threads.@threads for i in 1:nz
            z = zone_ids[i]
            part = partitions[z]

            lambdatheta = (lambdatheta_zone_buf[z] .= duals_theta_dense[z])
            lambdaV = (lambdaV_zone_buf[z] .= dualsV_dense[z])

            lin = _build_linear_bus_expr_ac!(
                part,
                shared_buses_sorted,
                shared_buses_idx,
                theta_consensus_target,
                vmag_consensus_target,
                lambdatheta,
                lambdaV,
                rho_theta[],
                rhoV[],
            )

            m = part.model::JuMP.Model
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

        theta_consensus_prev .= theta_consensus_vec
        vmag_consensus_prev .= vmag_consensus_vec
        compute_ac_bus_consensus_vec!(theta_consensus_vec, vmag_consensus_vec, partitions, shared_buses_sorted, shared_buses_idx)

        @inbounds for b in shared_buses_sorted
            j = shared_buses_idx[b]
            push!(get!(theta_avg_hist, b, Float64[]), theta_consensus_vec[j])
        end

        for z in zone_ids
            update_duals_bus_ac_vec!(
                duals_theta_dense[z],
                dualsV_dense[z],
                partitions[z],
                theta_consensus_vec,
                vmag_consensus_vec,
                rho_theta[],
                rhoV[],
                shared_buses_sorted,
                shared_buses_idx;
                α = α_used,
                theta_consensus_prev = theta_consensus_prev,
                vmag_consensus_prev = vmag_consensus_prev,
            )
        end

        r = primal_residual_bus_ac_vec(partitions, theta_consensus_vec, vmag_consensus_vec, shared_buses_sorted, shared_buses_idx)
        s = dual_residual_bus_ac_vec(theta_consensus_vec, vmag_consensus_vec, theta_consensus_prev, vmag_consensus_prev)

        epsp = eps_primal_bus_ac_vec(
            partitions, theta_consensus_vec, vmag_consensus_vec, shared_buses_sorted, shared_buses_idx, eps_abs, eps_rel
        )
        epsd = eps_dual_bus_ac_vec(theta_consensus_vec, vmag_consensus_vec, theta_consensus_prev, vmag_consensus_prev, eps_abs, eps_rel)

        push!(res_pri_hist, r)
        push!(res_du_hist, s)
        push!(eps_pri_hist, epsp)
        push!(eps_du_hist, epsd)

        _FAST_USE_LOGS && println("‖r‖ = $r | eps_pri = $epsp")
        _FAST_USE_LOGS && println("‖s‖ = $s | eps_du  = $epsd")

        push!(rho_theta_hist, rho_theta[])
        push!(rhoV_hist, rhoV[])

        rtheta, rV = primal_residual_bus_ac_components(partitions, theta_consensus_vec, vmag_consensus_vec, shared_buses_sorted, shared_buses_idx)
        stheta, sV = dual_residual_bus_ac_components(theta_consensus_vec, vmag_consensus_vec, theta_consensus_prev, vmag_consensus_prev)

        if r < epsp && s < epsd
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        if adaptive
            adaptive_rho_ac!(rho_theta, rtheta, stheta; mu_incr=config.mu_incr, mu_decr=config.mu_decr, tau_incr=config.tau_incr, tau_decr=config.tau_decr)
            adaptive_rho_ac!(rhoV, rV, sV; mu_incr=config.mu_incr, mu_decr=config.mu_decr, tau_incr=config.tau_incr, tau_decr=config.tau_decr)
        end
    end

    save_ac_bus_plots(
        res_pri_hist, res_du_hist, eps_pri_hist, eps_du_hist,
        rho_theta_hist, rhoV_hist, theta_avg_hist; stem = "admm_ac"
    )

    total_real_cost = sum(p.real_generation_cost for p in values(partitions))
    total_expr_cost = sum(p.generation_cost_model_expr for p in values(partitions))

    println("\n[ADMM-AC] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("rho_theta_final = $(rho_theta[]) | rhoV_final = $(rhoV[])")
    elapsed_s = (time_ns() - t0) / 1e9
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (bus/AC): $(_pretty_time(elapsed_s))")

    return total_real_cost, partitions, iter_done, elapsed_s
end


#################################
# ADMM: BRANCH/AC               #
#################################

function run_admm_ac_branch(
    config::DistributedOPFModel,
    sys::ADMMSystem;
    warm_start::Bool = false,  # mantido só pra compatibilidade (não usado aqui)
)
    t0 = time_ns()
    @assert config.base_model == AC "ADMM-AC-branch: base_model deve ser AC."

    partitions = build_ac_partitions(sys)
    cut_edges = sys.cut_edges

    rhoP = Ref(config.rho_flow_p)
    rhoQ = Ref(config.rho_flow_q)
    eps_abs, eps_rel = config.eps_primal, config.eps_dual
    adaptive = _adaptive_on(config)

    dualsP = Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))
    dualsQ = Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))

    zP = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zQ = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zP_prev = Dict{Int, Float64}()
    zQ_prev = Dict{Int, Float64}()

    res_pri_hist = Float64[]
    res_du_hist = Float64[]
    eps_pri_hist = Float64[]
    eps_du_hist = Float64[]
    rhoP_hist = Float64[]
    rhoQ_hist = Float64[]

    optimizer_factory = resolve_optimizer_ac(config)

    # init por zona (model + quad branch)
    for (zone_id, part) in partitions
        init_ac_branch_partition!(
            part,
            sys.parts[zone_id],
            optimizer_factory,
            rhoP[],
            rhoQ[],
            cut_edges,
        )
    end

    prev_rhoP, prev_rhoQ = rhoP[], rhoQ[]
    zone_ids = collect(keys(partitions))
    iter_done = 0

    for k in 1:config.max_iter
        iter_done = k
        _FAST_USE_LOGS && println("\n[ADMM-AC-branch] Iteração $k")

        if rhoP[] != prev_rhoP || rhoQ[] != prev_rhoQ
            for (_, part) in partitions
                _build_quad_branch_expr_ac!(part, rhoP[], rhoQ[], part.zone_edges)
            end
            prev_rhoP, prev_rhoQ = rhoP[], rhoQ[]
        end

        α_used = nothing

        nz = length(zone_ids)
        results = Vector{Any}(undef, nz)

        Threads.@threads for i in 1:nz
            z = zone_ids[i]
            part = partitions[z]

            lin = _build_linear_branch_expr_ac!(
                part,
                part.zone_edges,
                zP,
                zQ,
                dualsP[z],
                dualsQ[z],
                rhoP[],
                rhoQ[],
            )

            m = part.model::JuMP.Model
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
            dualsP,
            dualsQ,
            partitions,
            zP,
            zQ,
            cut_edges,
            rhoP[],
            rhoQ[];
            α = α_used,
            zP_prev = zP_prev,
            zQ_prev = zQ_prev,
        )

        r = primal_residual_branch_ac(partitions, zP, zQ, cut_edges)
        s = dual_residual_branch_ac(zP, zQ, zP_prev, zQ_prev, cut_edges)
        epsp = eps_primal_branch_ac(partitions, zP, zQ, cut_edges, eps_abs, eps_rel)
        epsd = eps_dual_branch_ac(zP, zQ, zP_prev, zQ_prev, cut_edges, eps_abs, eps_rel)

        push!(res_pri_hist, r)
        push!(res_du_hist, s)
        push!(eps_pri_hist, epsp)
        push!(eps_du_hist, epsd)

        _FAST_USE_LOGS && println("‖r‖ = $r | eps_pri = $epsp")
        _FAST_USE_LOGS && println("‖s‖ = $s | eps_du  = $epsd")

        push!(rhoP_hist, rhoP[])
        push!(rhoQ_hist, rhoQ[])

        if r < epsp && s < epsd
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        if adaptive
            adaptive_rho_ac!(rhoP, r, s; mu_incr=config.mu_incr, mu_decr=config.mu_decr, tau_incr=config.tau_incr, tau_decr=config.tau_decr)
            adaptive_rho_ac!(rhoQ, r, s; mu_incr=config.mu_incr, mu_decr=config.mu_decr, tau_incr=config.tau_incr, tau_decr=config.tau_decr)
        end
    end

    save_ac_branch_plots(
        res_pri_hist, res_du_hist, eps_pri_hist, eps_du_hist,
        rhoP_hist, rhoQ_hist; stem = "admm_ac_branch"
    )

    total_real_cost = sum(p.real_generation_cost for p in values(partitions))
    total_expr_cost = sum(p.generation_cost_model_expr for p in values(partitions))

    println("\n[ADMM-AC-branch] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("rhoP_final = $(rhoP[]) | rhoQ_final = $(rhoQ[])")
    elapsed_s = (time_ns() - t0) / 1e9
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (branch/AC): $(_pretty_time(elapsed_s))")

    return total_real_cost, partitions, iter_done, elapsed_s
end


#################################
# ADMM: HÍBRIDO/AC              #
#################################

function run_admm_ac_hybrid(config::DistributedOPFModel, sys::ADMMSystem)
    t0 = time_ns()
    @assert config.base_model == AC "Hybrid AC: base_model deve ser AC."

    partitions = build_ac_partitions(sys)
    shared_buses = sys.shared_buses
    cut_edges = sys.cut_edges

    rho_theta = Ref(config.rho_theta)
    rhoV = Ref(config.rho_v)
    rhoP = Ref(config.rho_flow_p)
    rhoQ = Ref(config.rho_flow_q)
    eps_abs, eps_rel = config.eps_primal, config.eps_dual
    adaptive = _adaptive_on(config)

    shared_buses_sorted, shared_buses_idx = prepare_dense_maps(shared_buses)

    duals_theta_dense = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )
    dualsV_dense = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )
    dualsP = Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))
    dualsQ = Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))

    theta_consensus_vec = zeros(Float64, length(shared_buses_sorted))
    vmag_consensus_vec = ones(Float64, length(shared_buses_sorted))
    theta_consensus_prev = similar(theta_consensus_vec)
    vmag_consensus_prev = similar(vmag_consensus_vec)

    zP = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zQ = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zP_prev = Dict{Int, Float64}()
    zQ_prev = Dict{Int, Float64}()

    res_pri_hist = Float64[]
    res_du_hist = Float64[]
    eps_pri_hist = Float64[]
    eps_du_hist = Float64[]
    rho_theta_hist = Float64[]
    rhoV_hist = Float64[]
    rhoP_hist = Float64[]
    rhoQ_hist = Float64[]
    theta_avg_hist = Dict{Int, Vector{Float64}}()

    optimizer_factory = resolve_optimizer_ac(config)

    # init por zona (hotpaths + zone_edges + model + quads bus/branch)
    for (zone_id, part) in partitions
        init_ac_hybrid_partition!(
            part,
            sys.parts[zone_id],
            optimizer_factory,
            rho_theta[],
            rhoV[],
            rhoP[],
            rhoQ[],
            cut_edges,
        )
    end

    prev_rho_theta, prev_rhoV = rho_theta[], rhoV[]
    prev_rhoP, prev_rhoQ = rhoP[], rhoQ[]
    zone_ids = collect(keys(partitions))
    lambdatheta_zone_buf = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in zone_ids
    )
    lambdaV_zone_buf = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in zone_ids
    )

    iter_done = 0

    for k in 1:config.max_iter
        iter_done = k
        _FAST_USE_LOGS && println("\n[ADMM-AC-hybrid] Iteração $k")

        if rho_theta[] != prev_rho_theta || rhoV[] != prev_rhoV
            for (_, part) in partitions
                _build_quad_bus_expr_ac!(part, rho_theta[], rhoV[])
            end
            prev_rho_theta, prev_rhoV = rho_theta[], rhoV[]
        end
        if rhoP[] != prev_rhoP || rhoQ[] != prev_rhoQ
            for (_, part) in partitions
                _build_quad_branch_expr_ac!(part, rhoP[], rhoQ[], part.zone_edges)
            end
            prev_rhoP, prev_rhoQ = rhoP[], rhoQ[]
        end

        α_used = nothing
        theta_consensus_target = theta_consensus_vec
        vmag_consensus_target = vmag_consensus_vec

        nz = length(zone_ids)
        results = Vector{Any}(undef, nz)

        Threads.@threads for i in 1:nz
            z = zone_ids[i]
            part = partitions[z]

            lambdatheta = (lambdatheta_zone_buf[z] .= duals_theta_dense[z])
            lambdaV = (lambdaV_zone_buf[z] .= dualsV_dense[z])

            lin_bus = _build_linear_bus_expr_ac!(
                part,
                shared_buses_sorted,
                shared_buses_idx,
                theta_consensus_target,
                vmag_consensus_target,
                lambdatheta,
                lambdaV,
                rho_theta[],
                rhoV[],
            )
            lin_br = _build_linear_branch_expr_ac!(
                part,
                part.zone_edges,
                zP,
                zQ,
                dualsP[z],
                dualsQ[z],
                rhoP[],
                rhoQ[],
            )

            m = part.model::JuMP.Model
            @objective(
                m,
                Min,
                part.obj_base + part.quad_bus_expr_ac + part.quad_branch_expr_ac +
                (lin_bus + lin_br),
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

        theta_consensus_prev .= theta_consensus_vec
        vmag_consensus_prev .= vmag_consensus_vec
        zP_prev = copy(zP)
        zQ_prev = copy(zQ)

        compute_ac_bus_consensus_vec!(theta_consensus_vec, vmag_consensus_vec, partitions, shared_buses_sorted, shared_buses_idx)
        @inbounds for b in shared_buses_sorted
            j = shared_buses_idx[b]
            push!(get!(theta_avg_hist, b, Float64[]), theta_consensus_vec[j])
        end
        compute_ac_branch_consensus!(zP, zQ, partitions, cut_edges)

        for z in zone_ids
            update_duals_bus_ac_vec!(
                duals_theta_dense[z],
                dualsV_dense[z],
                partitions[z],
                theta_consensus_vec,
                vmag_consensus_vec,
                rho_theta[],
                rhoV[],
                shared_buses_sorted,
                shared_buses_idx;
                α = α_used,
                theta_consensus_prev = theta_consensus_prev,
                vmag_consensus_prev = vmag_consensus_prev,
            )
        end
        update_duals_branch_ac!(
            dualsP,
            dualsQ,
            partitions,
            zP,
            zQ,
            cut_edges,
            rhoP[],
            rhoQ[];
            α = α_used,
            zP_prev = zP_prev,
            zQ_prev = zQ_prev,
        )

        r_bus = primal_residual_bus_ac_vec(partitions, theta_consensus_vec, vmag_consensus_vec, shared_buses_sorted, shared_buses_idx)
        r_br  = primal_residual_branch_ac(partitions, zP, zQ, cut_edges)
        r_tot = hypot(r_bus, r_br)

        s_bus = dual_residual_bus_ac_vec(theta_consensus_vec, vmag_consensus_vec, theta_consensus_prev, vmag_consensus_prev)
        s_br  = dual_residual_branch_ac(zP, zQ, zP_prev, zQ_prev, cut_edges)
        s_tot = hypot(s_bus, s_br)

        epsb_pri = eps_primal_bus_ac_vec(
            partitions, theta_consensus_vec, vmag_consensus_vec, shared_buses_sorted, shared_buses_idx, eps_abs, eps_rel
        )
        epse_pri = eps_primal_branch_ac(partitions, zP, zQ, cut_edges, eps_abs, eps_rel)
        eps_pri  = hypot(epsb_pri, epse_pri)

        epsb_du = eps_dual_bus_ac_vec(theta_consensus_vec, vmag_consensus_vec, theta_consensus_prev, vmag_consensus_prev, eps_abs, eps_rel)
        epse_du = eps_dual_branch_ac(zP, zQ, zP_prev, zQ_prev, cut_edges, eps_abs, eps_rel)
        eps_du  = hypot(epsb_du, epse_du)

        push!(res_pri_hist, r_tot)
        push!(res_du_hist, s_tot)
        push!(eps_pri_hist, eps_pri)
        push!(eps_du_hist, eps_du)

        _FAST_USE_LOGS && println("‖r‖ = $r_tot | eps_pri = $eps_pri")
        _FAST_USE_LOGS && println("‖s‖ = $s_tot | eps_du  = $eps_du")

        push!(rho_theta_hist, rho_theta[])
        push!(rhoV_hist, rhoV[])
        push!(rhoP_hist, rhoP[])
        push!(rhoQ_hist, rhoQ[])

        if r_tot < eps_pri && s_tot < eps_du
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        if adaptive
            rtheta, rV = primal_residual_bus_ac_components(partitions, theta_consensus_vec, vmag_consensus_vec, shared_buses_sorted, shared_buses_idx)
            stheta, sV = dual_residual_bus_ac_components(theta_consensus_vec, vmag_consensus_vec, theta_consensus_prev, vmag_consensus_prev)

            adaptive_rho_ac!(rho_theta, rtheta, stheta; mu_incr=config.mu_incr, mu_decr=config.mu_decr, tau_incr=config.tau_incr, tau_decr=config.tau_decr)
            adaptive_rho_ac!(rhoV, rV, sV; mu_incr=config.mu_incr, mu_decr=config.mu_decr, tau_incr=config.tau_incr, tau_decr=config.tau_decr)

            adaptive_rho_ac!(rhoP, r_br, s_br; mu_incr=config.mu_incr, mu_decr=config.mu_decr, tau_incr=config.tau_incr, tau_decr=config.tau_decr)
            adaptive_rho_ac!(rhoQ, r_br, s_br; mu_incr=config.mu_incr, mu_decr=config.mu_decr, tau_incr=config.tau_incr, tau_decr=config.tau_decr)
        end
    end

    # mantém seus plots/saves exatamente como você já tinha (não re-colei aqui)
    # (cole aqui embaixo o mesmo bloco de plots que você mandou)

    total_real_cost = sum(p.real_generation_cost for p in values(partitions))
    total_expr_cost = sum(p.generation_cost_model_expr for p in values(partitions))

    println("\n[ADMM-AC-hybrid] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("rho_theta_final = $(rho_theta[]) | rhoV_final = $(rhoV[]) | rhoP_final = $(rhoP[]) | rhoQ_final = $(rhoQ[])")
    elapsed_s = (time_ns() - t0) / 1e9
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (hybrid/AC): $(_pretty_time(elapsed_s))")

    return total_real_cost, partitions, iter_done, elapsed_s
end