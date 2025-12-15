function run_admm_ac(config::DistributedOPFModel, sys::ADMMSystem)
    t0 = time_ns()
    @assert config.base_model == AC "ADMM-AC: base_model deve ser AC."

    partitions = build_ac_partitions(sys)
    shared_buses = sys.shared_buses

    ρθ = Ref(config.rho_theta)
    ρV = Ref(config.rho_v)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual
    adaptive = _adaptive_on(config)

    for (_, part) in partitions
        prepare_hotpaths_ac!(part)
    end
    shared_buses_sorted, b2pos = prepare_dense_maps(shared_buses)

    dualsθ_dense = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )
    dualsV_dense = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )

    θ̄_vec = zeros(Float64, length(shared_buses_sorted))
    V̄_vec = ones(Float64, length(shared_buses_sorted))
    θ̄_prev = similar(θ̄_vec)
    V̄_prev = similar(V̄_vec)

    res_pri_hist = Float64[]
    res_du_hist = Float64[]
    eps_pri_hist = Float64[]
    eps_du_hist = Float64[]
    rhoθ_hist = Float64[]
    rhoV_hist = Float64[]
    θ_avg_hist = Dict{Int, Vector{Float64}}()

    optimizer_factory = resolve_optimizer_ac(config)

    # build inicial + termo quadrático em bus
    for (zone_id, part) in partitions
        sys_zone = sys.parts[zone_id]

        m = build_model(
            sys_zone,
            AC;
            optimizer = optimizer_factory,
        )

        part.model = m
        part.obj_base = m[:generation_cost]
        part.theta_var = m[:theta]
        part.vmag_var = m[:v_magnitude]
        part.has_solution = false

        _set_model_silent!(m)
        for I in keys(part.theta_var)
            set_start_value(part.theta_var[I], 0.0)
        end
        for I in keys(part.vmag_var)
            set_start_value(part.vmag_var[I], 1.0)
        end

        _build_quad_bus_expr_ac!(part, ρθ[], ρV[])
        @objective(m, Min, part.obj_base + part.quad_bus_expr_ac)
    end

    prev_ρθ, prev_ρV = ρθ[], ρV[]
    zone_ids = collect(keys(partitions))
    λθ_zone_buf = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in zone_ids
    )
    λV_zone_buf = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in zone_ids
    )

    iter_done = 0

    for k in 1:config.max_iter
        iter_done = k
        _FAST_USE_LOGS && println("\n[ADMM-AC] Iteração $k")

        if ρθ[] != prev_ρθ || ρV[] != prev_ρV
            for (_, part) in partitions
                _build_quad_bus_expr_ac!(part, ρθ[], ρV[])
            end
            prev_ρθ, prev_ρV = ρθ[], ρV[]
        end

        # sem over-relaxation / FISTA
        α_used = nothing
        θ̄_target = θ̄_vec
        V̄_target = V̄_vec

        nz = length(zone_ids)
        results = Vector{Any}(undef, nz)

        Threads.@threads for i in 1:nz
            z = zone_ids[i]
            part = partitions[z]

            λθ = (λθ_zone_buf[z] .= dualsθ_dense[z])
            λV = (λV_zone_buf[z] .= dualsV_dense[z])

            lin = _build_linear_bus_expr_ac!(
                part,
                shared_buses_sorted,
                b2pos,
                θ̄_target,
                V̄_target,
                λθ,
                λV,
                ρθ[],
                ρV[],
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

        θ̄_prev .= θ̄_vec
        V̄_prev .= V̄_vec
        compute_ac_bus_consensus_vec!(
            θ̄_vec,
            V̄_vec,
            partitions,
            shared_buses_sorted,
            b2pos,
        )
        @inbounds for b in shared_buses_sorted
            j = b2pos[b]
            push!(get!(θ_avg_hist, b, Float64[]), θ̄_vec[j])
        end

        for z in zone_ids
            update_duals_bus_ac_vec!(
                dualsθ_dense[z],
                dualsV_dense[z],
                partitions[z],
                θ̄_vec,
                V̄_vec,
                ρθ[],
                ρV[],
                shared_buses_sorted,
                b2pos;
                α = α_used,
                θ̄_prev = θ̄_prev,
                V̄_prev = V̄_prev,
            )
        end

        r = primal_residual_bus_ac_vec(
            partitions,
            θ̄_vec,
            V̄_vec,
            shared_buses_sorted,
            b2pos,
        )
        s = dual_residual_bus_ac_vec(θ̄_vec, V̄_vec, θ̄_prev, V̄_prev)
        εp = eps_primal_bus_ac_vec(
            partitions,
            θ̄_vec,
            V̄_vec,
            shared_buses_sorted,
            b2pos,
            ϵ_abs,
            ϵ_rel,
        )
        εd = eps_dual_bus_ac_vec(
            θ̄_vec,
            V̄_vec,
            θ̄_prev,
            V̄_prev,
            ϵ_abs,
            ϵ_rel,
        )

        push!(res_pri_hist, r)
        push!(res_du_hist, s)
        push!(eps_pri_hist, εp)
        push!(eps_du_hist, εd)

        _FAST_USE_LOGS && println("‖r‖ = $r | ε_pri = $εp")
        _FAST_USE_LOGS && println("‖s‖ = $s | ε_du  = $εd")

        push!(rhoθ_hist, ρθ[])
        push!(rhoV_hist, ρV[])

        rθ, rV = primal_residual_bus_ac_components(
            partitions,
            θ̄_vec,
            V̄_vec,
            shared_buses_sorted,
            b2pos,
        )
        sθ, sV = dual_residual_bus_ac_components(
            θ̄_vec,
            V̄_vec,
            θ̄_prev,
            V̄_prev,
        )

        if r < εp && s < εd
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        if adaptive
            adaptive_rho_ac!(
                ρθ,
                rθ,
                sθ;
                μ_incr = config.mu_incr,
                μ_decr = config.mu_decr,
                τ_incr = config.tau_incr,
                τ_decr = config.tau_decr,
            )
            adaptive_rho_ac!(
                ρV,
                rV,
                sV;
                μ_incr = config.mu_incr,
                μ_decr = config.mu_decr,
                τ_incr = config.tau_incr,
                τ_decr = config.tau_decr,
            )
        end
    end

    save_ac_bus_plots(
        res_pri_hist,
        res_du_hist,
        eps_pri_hist,
        eps_du_hist,
        rhoθ_hist,
        rhoV_hist,
        θ_avg_hist;
        stem = "admm_ac",
    )

    total_real_cost = sum(p.real_generation_cost for p in values(partitions))
    total_expr_cost = sum(p.generation_cost_model_expr for p in values(partitions))

    println("\n[ADMM-AC] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("ρθ_final = $(ρθ[]) | ρV_final = $(ρV[])")
    elapsed_s = (time_ns() - t0) / 1e9
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (bus/AC): $(_pretty_time(elapsed_s))")

    return total_real_cost, partitions, iter_done, elapsed_s
end

#################################
# ADMM: BRANCH (AC) – sys       #
#################################

function run_admm_ac_branch(
    config::DistributedOPFModel,
    sys::ADMMSystem;
    warm_start::Bool = false,
)
    t0 = time_ns()
    @assert config.base_model == AC "ADMM-AC-branch: base_model deve ser AC."

    partitions = build_ac_partitions(sys)
    cut_edges = sys.cut_edges

    for (_, part) in partitions
        part.zone_edges = _zone_cut_edges(part, cut_edges)
    end

    ρP = Ref(config.rho_flow_p)
    ρQ = Ref(config.rho_flow_q)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual
    adaptive = _adaptive_on(config)

    dualsP =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))
    dualsQ =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))
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

    # build inicial + termo quadrático em branch
    for (zone_id, part) in partitions
        sys_zone = sys.parts[zone_id]

        m = build_model(
            sys_zone,
            AC;
            optimizer = optimizer_factory,
        )
        part.model = m
        part.obj_base = m[:generation_cost]
        part.theta_var = m[:theta]
        part.vmag_var = m[:v_magnitude]
        part.has_solution = false

        _set_model_silent!(m)
        for I in keys(part.theta_var)
            set_start_value(part.theta_var[I], 0.0)
        end
        for I in keys(part.vmag_var)
            set_start_value(part.vmag_var[I], 1.0)
        end

        _build_quad_branch_expr_ac!(part, ρP[], ρQ[], part.zone_edges)
        @objective(m, Min, part.obj_base + part.quad_branch_expr_ac)
    end

    prev_ρP, prev_ρQ = ρP[], ρQ[]
    zone_ids = collect(keys(partitions))
    iter_done = 0

    for k in 1:config.max_iter
        iter_done = k
        _FAST_USE_LOGS && println("\n[ADMM-AC-branch] Iteração $k")

        if ρP[] != prev_ρP || ρQ[] != prev_ρQ
            for (_, part) in partitions
                _build_quad_branch_expr_ac!(part, ρP[], ρQ[], part.zone_edges)
            end
            prev_ρP, prev_ρQ = ρP[], ρQ[]
        end

        # sem over-relaxation
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
                ρP[],
                ρQ[],
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
            ρP[],
            ρQ[];
            α = α_used,
            zP_prev = zP_prev,
            zQ_prev = zQ_prev,
        )

        r = primal_residual_branch_ac(partitions, zP, zQ, cut_edges)
        s = dual_residual_branch_ac(zP, zQ, zP_prev, zQ_prev, cut_edges)
        εp = eps_primal_branch_ac(partitions, zP, zQ, cut_edges, ϵ_abs, ϵ_rel)
        εd = eps_dual_branch_ac(zP, zQ, zP_prev, zQ_prev, cut_edges, ϵ_abs, ϵ_rel)

        push!(res_pri_hist, r)
        push!(res_du_hist, s)
        push!(eps_pri_hist, εp)
        push!(eps_du_hist, εd)

        _FAST_USE_LOGS && println("‖r‖ = $r | ε_pri = $εp")
        _FAST_USE_LOGS && println("‖s‖ = $s | ε_du  = $εd")

        push!(rhoP_hist, ρP[])
        push!(rhoQ_hist, ρQ[])

        if r < εp && s < εd
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        if adaptive
            adaptive_rho_ac!(
                ρP,
                r,
                s;
                μ_incr = config.mu_incr,
                μ_decr = config.mu_decr,
                τ_incr = config.tau_incr,
                τ_decr = config.tau_decr,
            )
            adaptive_rho_ac!(
                ρQ,
                r,
                s;
                μ_incr = config.mu_incr,
                μ_decr = config.mu_decr,
                τ_incr = config.tau_incr,
                τ_decr = config.tau_decr,
            )
        end
    end

    save_ac_branch_plots(
        res_pri_hist,
        res_du_hist,
        eps_pri_hist,
        eps_du_hist,
        rhoP_hist,
        rhoQ_hist;
        stem = "admm_ac_branch",
    )

    total_real_cost = sum(p.real_generation_cost for p in values(partitions))
    total_expr_cost = sum(p.generation_cost_model_expr for p in values(partitions))

    println("\n[ADMM-AC-branch] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("ρP_final = $(ρP[]) | ρQ_final = $(ρQ[])")
    elapsed_s = (time_ns() - t0) / 1e9
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (branch/AC): $(_pretty_time(elapsed_s))")

    return total_real_cost, partitions, iter_done, elapsed_s
end

#################################
# ADMM: HÍBRIDO (AC) – sys      #
#################################

function run_admm_ac_hybrid(config::DistributedOPFModel, sys::ADMMSystem)
    t0 = time_ns()
    @assert config.base_model == AC "Hybrid AC: base_model deve ser AC."

    partitions = build_ac_partitions(sys)
    shared_buses = sys.shared_buses
    cut_edges = sys.cut_edges

    for (_, part) in partitions
        part.zone_edges = _zone_cut_edges(part, cut_edges)
        prepare_hotpaths_ac!(part)
    end

    ρθ = Ref(config.rho_theta)
    ρV = Ref(config.rho_v)
    ρP = Ref(config.rho_flow_p)
    ρQ = Ref(config.rho_flow_q)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual
    adaptive = _adaptive_on(config)

    adaptive_bus = adaptive
    adaptive_br = adaptive

    shared_buses_sorted, b2pos = prepare_dense_maps(shared_buses)

    dualsθ_dense = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )
    dualsV_dense = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )
    dualsP =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))
    dualsQ =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))

    θ̄_vec = zeros(Float64, length(shared_buses_sorted))
    V̄_vec = ones(Float64, length(shared_buses_sorted))
    θ̄_prev = similar(θ̄_vec)
    V̄_prev = similar(V̄_vec)

    zP = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zQ = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zP_prev = Dict{Int, Float64}()
    zQ_prev = Dict{Int, Float64}()

    res_pri_hist = Float64[]
    res_du_hist = Float64[]
    eps_pri_hist = Float64[]
    eps_du_hist = Float64[]
    rhoθ_hist = Float64[]
    rhoV_hist = Float64[]
    rhoP_hist = Float64[]
    rhoQ_hist = Float64[]
    θ_avg_hist = Dict{Int, Vector{Float64}}()

    optimizer_factory = resolve_optimizer_ac(config)

    # build inicial + termos quadráticos
    for (zone_id, part) in partitions
        sys_zone = sys.parts[zone_id]
        m = build_model(
            sys_zone,
            AC;
            optimizer = optimizer_factory,
        )
        part.model = m
        part.obj_base = m[:generation_cost]
        part.theta_var = m[:theta]
        part.vmag_var = m[:v_magnitude]
        part.has_solution = false

        _set_model_silent!(m)
        for I in keys(part.theta_var)
            set_start_value(part.theta_var[I], 0.0)
        end
        for I in keys(part.vmag_var)
            set_start_value(part.vmag_var[I], 1.0)
        end

        _build_quad_bus_expr_ac!(part, ρθ[], ρV[])
        _build_quad_branch_expr_ac!(part, ρP[], ρQ[], part.zone_edges)
        @objective(
            m,
            Min,
            part.obj_base + part.quad_bus_expr_ac + part.quad_branch_expr_ac,
        )
    end

    prev_ρθ, prev_ρV = ρθ[], ρV[]
    prev_ρP, prev_ρQ = ρP[], ρQ[]
    zone_ids = collect(keys(partitions))
    λθ_zone_buf = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in zone_ids
    )
    λV_zone_buf = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in zone_ids
    )

    iter_done = 0

    for k in 1:config.max_iter
        iter_done = k
        _FAST_USE_LOGS && println("\n[ADMM-AC-hybrid] Iteração $k")

        if ρθ[] != prev_ρθ || ρV[] != prev_ρV
            for (_, part) in partitions
                _build_quad_bus_expr_ac!(part, ρθ[], ρV[])
            end
            prev_ρθ, prev_ρV = ρθ[], ρV[]
        end
        if ρP[] != prev_ρP || ρQ[] != prev_ρQ
            for (_, part) in partitions
                _build_quad_branch_expr_ac!(part, ρP[], ρQ[], part.zone_edges)
            end
            prev_ρP, prev_ρQ = ρP[], ρQ[]
        end

        # sem over-relaxation / aceleração
        α_used = nothing
        θ̄_target = θ̄_vec
        V̄_target = V̄_vec

        nz = length(zone_ids)
        results = Vector{Any}(undef, nz)
        Threads.@threads for i in 1:nz
            z = zone_ids[i]
            part = partitions[z]

            λθ = (λθ_zone_buf[z] .= dualsθ_dense[z])
            λV = (λV_zone_buf[z] .= dualsV_dense[z])

            lin_bus = _build_linear_bus_expr_ac!(
                part,
                shared_buses_sorted,
                b2pos,
                θ̄_target,
                V̄_target,
                λθ,
                λV,
                ρθ[],
                ρV[],
            )
            lin_br = _build_linear_branch_expr_ac!(
                part,
                part.zone_edges,
                zP,
                zQ,
                dualsP[z],
                dualsQ[z],
                ρP[],
                ρQ[],
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

        θ̄_prev .= θ̄_vec
        V̄_prev .= V̄_vec
        zP_prev = copy(zP)
        zQ_prev = copy(zQ)

        compute_ac_bus_consensus_vec!(
            θ̄_vec,
            V̄_vec,
            partitions,
            shared_buses_sorted,
            b2pos,
        )
        for b in shared_buses_sorted
            j = b2pos[b]
            push!(get!(θ_avg_hist, b, Float64[]), θ̄_vec[j])
        end
        compute_ac_branch_consensus!(zP, zQ, partitions, cut_edges)

        for z in zone_ids
            update_duals_bus_ac_vec!(
                dualsθ_dense[z],
                dualsV_dense[z],
                partitions[z],
                θ̄_vec,
                V̄_vec,
                ρθ[],
                ρV[],
                shared_buses_sorted,
                b2pos;
                α = α_used,
                θ̄_prev = θ̄_prev,
                V̄_prev = V̄_prev,
            )
        end
        update_duals_branch_ac!(
            dualsP,
            dualsQ,
            partitions,
            zP,
            zQ,
            cut_edges,
            ρP[],
            ρQ[];
            α = α_used,
            zP_prev = zP_prev,
            zQ_prev = zQ_prev,
        )

        r_bus = primal_residual_bus_ac_vec(
            partitions,
            θ̄_vec,
            V̄_vec,
            shared_buses_sorted,
            b2pos,
        )
        r_br = primal_residual_branch_ac(partitions, zP, zQ, cut_edges)
        r_tot = hypot(r_bus, r_br)

        s_bus = dual_residual_bus_ac_vec(θ̄_vec, V̄_vec, θ̄_prev, V̄_prev)
        s_br = dual_residual_branch_ac(zP, zQ, zP_prev, zQ_prev, cut_edges)
        s_tot = hypot(s_bus, s_br)

        εb_pri = eps_primal_bus_ac_vec(
            partitions,
            θ̄_vec,
            V̄_vec,
            shared_buses_sorted,
            b2pos,
            ϵ_abs,
            ϵ_rel,
        )
        εe_pri = eps_primal_branch_ac(partitions, zP, zQ, cut_edges, ϵ_abs, ϵ_rel)
        ε_pri = hypot(εb_pri, εe_pri)

        εb_du = eps_dual_bus_ac_vec(θ̄_vec, V̄_vec, θ̄_prev, V̄_prev, ϵ_abs, ϵ_rel)
        εe_du = eps_dual_branch_ac(zP, zQ, zP_prev, zQ_prev, cut_edges, ϵ_abs, ϵ_rel)
        ε_du = hypot(εb_du, εe_du)

        push!(res_pri_hist, r_tot)
        push!(res_du_hist, s_tot)
        push!(eps_pri_hist, ε_pri)
        push!(eps_du_hist, ε_du)

        _FAST_USE_LOGS && println("‖r‖ = $r_tot | ε_pri = $ε_pri")
        _FAST_USE_LOGS && println("‖s‖ = $s_tot | ε_du  = $ε_du")

        push!(rhoθ_hist, ρθ[])
        push!(rhoV_hist, ρV[])
        push!(rhoP_hist, ρP[])
        push!(rhoQ_hist, ρQ[])

        if r_tot < ε_pri && s_tot < ε_du
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        if adaptive_bus
            rθ, rV = primal_residual_bus_ac_components(
                partitions,
                θ̄_vec,
                V̄_vec,
                shared_buses_sorted,
                b2pos,
            )
            sθ, sV = dual_residual_bus_ac_components(
                θ̄_vec,
                V̄_vec,
                θ̄_prev,
                V̄_prev,
            )
            adaptive_rho_ac!(
                ρθ,
                rθ,
                sθ;
                μ_incr = config.mu_incr,
                μ_decr = config.mu_decr,
                τ_incr = config.tau_incr,
                τ_decr = config.tau_decr,
            )
            adaptive_rho_ac!(
                ρV,
                rV,
                sV;
                μ_incr = config.mu_incr,
                μ_decr = config.mu_decr,
                τ_incr = config.tau_incr,
                τ_decr = config.tau_decr,
            )
        end
        if adaptive_br
            adaptive_rho_ac!(
                ρP,
                r_br,
                s_br;
                μ_incr = config.mu_incr,
                μ_decr = config.mu_decr,
                τ_incr = config.tau_incr,
                τ_decr = config.tau_decr,
            )
            adaptive_rho_ac!(
                ρQ,
                r_br,
                s_br;
                μ_incr = config.mu_incr,
                μ_decr = config.mu_decr,
                τ_incr = config.tau_incr,
                τ_decr = config.tau_decr,
            )
        end
    end

    plot_admm_convergence(
        res_pri_hist,
        res_du_hist,
        eps_pri_hist,
        eps_du_hist,
        "plot_convergencia_admm_ac_hybrid.html",
    )
    PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(;
                    x = 1:length(rhoθ_hist),
                    y = rhoθ_hist,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρθ",
                );
                p.layout.title = "Evolução de ρθ";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "ρθ";
                p.layout.yaxis.type = "log";
                p
            )
        )(),
        "plot_rho_theta_ac_hybrid.html",
    )
    PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(;
                    x = 1:length(rhoV_hist),
                    y = rhoV_hist,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρV",
                );
                p.layout.title = "Evolução de ρV";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "ρV";
                p.layout.yaxis.type = "log";
                p
            )
        )(),
        "plot_rho_v_ac_hybrid.html",
    )
    PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(;
                    x = 1:length(rhoP_hist),
                    y = rhoP_hist,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρP",
                );
                p.layout.title = "Evolução de ρP";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "ρP";
                p.layout.yaxis.type = "log";
                p
            )
        )(),
        "plot_rho_p_ac_hybrid.html",
    )
    PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(;
                    x = 1:length(rhoQ_hist),
                    y = rhoQ_hist,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρQ",
                );
                p.layout.title = "Evolução de ρQ";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "ρQ";
                p.layout.yaxis.type = "log";
                p
            )
        )(),
        "plot_rho_q_ac_hybrid.html",
    )

    plot_theta_avg_heatmap(θ_avg_hist, "plot_theta_avg_heatmap_ac_hybrid.html")

    total_real_cost = sum(p.real_generation_cost for p in values(partitions))
    total_expr_cost = sum(p.generation_cost_model_expr for p in values(partitions))

    println("\n[ADMM-AC-hybrid] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println(
        "ρθ_final = $(ρθ[]) | ρV_final = $(ρV[]) | ρP_final = $(ρP[]) | ρQ_final = $(ρQ[])",
    )
    elapsed_s = (time_ns() - t0) / 1e9
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (hybrid/AC): $(_pretty_time(elapsed_s))")

    return total_real_cost, partitions, iter_done, elapsed_s
end