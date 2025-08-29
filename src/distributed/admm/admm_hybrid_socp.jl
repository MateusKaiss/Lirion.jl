function run_admm_socp_hybrid(config::DistributedOPFModel, data::Dict{String, Any})
    @assert config.base_model == SOCP "ADMM-SOCP-hybrid: base_model deve ser SOCP."

    partitions, shared_buses = prepare_admm_partitions(data)
    cut_edges = _find_cut_edges(data)

    # ---------------- duais ----------------
    duals_v2 = Dict{Int, Dict{Int, Float64}}()   # λV2[zone][bus]
    dualsP = Dict{Int, Dict{Int, Float64}}()   # λP [zone][edge]
    dualsQ = Dict{Int, Dict{Int, Float64}}()   # λQ [zone][edge]
    dualsR = Dict{Int, Dict{Int, Float64}}()   # λR [zone][edge]  (wr = Re(W))
    dualsT = Dict{Int, Dict{Int, Float64}}()   # λT [zone][edge]  (wi = Im(W))
    for z in keys(partitions)
        duals_v2[z] = Dict{Int, Float64}()
        dualsP[z] = Dict{Int, Float64}()
        dualsQ[z] = Dict{Int, Float64}()
        dualsR[z] = Dict{Int, Float64}()
        dualsT[z] = Dict{Int, Float64}()
    end

    # --------------- consensos --------------
    v2_avg = Dict{Int, Float64}()
    v2_prev = Dict{Int, Float64}()

    zP = Dict{Int, Float64}()
    zP_prev = Dict{Int, Float64}()
    zQ = Dict{Int, Float64}()
    zQ_prev = Dict{Int, Float64}()
    zR = Dict{Int, Float64}()
    zR_prev = Dict{Int, Float64}()
    zT = Dict{Int, Float64}()
    zT_prev = Dict{Int, Float64}()

    for e in cut_edges
        zP[e] = 0.0
        zQ[e] = 0.0
        zR[e] = 0.0
        zT[e] = 0.0
    end

    # --------------- penalidades ------------
    ρV2 = config.rho_v        # para V^2 (bus)
    ρP = config.rho_flow_p   # para P (branch)
    ρQ = config.rho_flow_q   # para Q (branch)
    ρR = config.rho_flow_p   # para R=Re(W) (branch)
    ρT = config.rho_flow_q   # para T=Im(W) (branch)

    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual

    # --------------- histórico --------------
    res_pri_hist = Float64[]
    res_du_hist = Float64[]
    eps_pri_hist = Float64[]
    eps_du_hist = Float64[]
    rhoV2_hist = Float64[]
    rhoP_hist = Float64[]
    rhoQ_hist = Float64[]
    rhoR_hist = Float64[]
    rhoT_hist = Float64[]

    μ_incr, μ_decr = 2.0, 2.0
    τ_incr, τ_decr = 1.05, 0.98

    @info "[ADMM-SOCP-hybrid] ρV2=$(ρV2), ρP=$(ρP), ρQ=$(ρQ), ρR=$(ρR), ρT=$(ρT) ; eps_pri=$(ϵ_abs) eps_du=$(ϵ_rel)"

    for k in 1:config.max_iter
        println("\n[ADMM-SOCP-hybrid] Iteração $k")

        # === 1) subproblemas locais (SOCP) ===
        for (zone_id, subdata) in partitions
            local_model = build_model(subdata, SOCP; optimizer = SCS.Optimizer)
            set_optimizer_attribute(local_model, "eps_abs", 1e-5)
            set_optimizer_attribute(local_model, "eps_rel", 1e-5)
            set_optimizer_attribute(local_model, "max_iters", 150_000)
            set_optimizer_attribute(local_model, "acceleration_lookback", 10)
            set_optimizer_attribute(local_model, "acceleration_interval", 10)
            set_optimizer_attribute(local_model, "alpha", 1.8)
            set_optimizer_attribute(local_model, "normalize", true)
            set_optimizer_attribute(local_model, "scale", 1.0)

            @expression(local_model, admm_penalty, 0.0)

            # --- BUS (V^2) ---
            bus_ids = subdata["bus"]["bus_i"]
            for b in shared_buses
                if b in bus_ids
                    V2loc = local_model[:v_magnitude_sq][b, 1]
                    V2bar = get(v2_avg, b, 1.0)
                    λV2 = get(duals_v2[zone_id], b, 0.0)
                    admm_penalty += λV2 * (V2loc - V2bar) + (ρV2 / 2) * (V2loc - V2bar)^2
                end
            end

            # --- BRANCH (P, Q, R, T) nas arestas de corte ---
            zone_edges = _zone_cut_edges(subdata, cut_edges)
            for e in zone_edges
                l_local = findfirst(==(e), subdata["branch_local_indices"])
                @assert l_local !== nothing "Aresta global $e não encontrada na partição $zone_id"

                P_loc, Q_loc = _socp_flow_from_exprs_existing(local_model, l_local)
                R_loc, T_loc = _rt_from_exprs_socp(local_model, l_local)

                λP = get(dualsP[zone_id], e, 0.0)
                λQ = get(dualsQ[zone_id], e, 0.0)
                λR = get(dualsR[zone_id], e, 0.0)
                λT = get(dualsT[zone_id], e, 0.0)

                admm_penalty += λP * (P_loc - zP[e]) + (ρP / 2) * (P_loc - zP[e])^2
                admm_penalty += λQ * (Q_loc - zQ[e]) + (ρQ / 2) * (Q_loc - zQ[e])^2
                admm_penalty += λR * (R_loc - zR[e]) + (ρR / 2) * (R_loc - zR[e])^2
                admm_penalty += λT * (T_loc - zT[e]) + (ρT / 2) * (T_loc - zT[e])^2
            end

            @objective(local_model, Min, local_model[:generation_cost] + admm_penalty)
            optimize!(local_model)

            subdata["model"] = local_model
            subdata["generation_cost_model_expr"] = value(local_model[:generation_cost])
            subdata["real_generation_cost"] =
                compute_real_generation_cost(local_model, subdata["objective_params"])
        end

        # === 2) atualizar consensos ===
        v2_prev = isempty(v2_avg) ? v2_avg : deepcopy(v2_avg)
        zP_prev = deepcopy(zP)
        zQ_prev = deepcopy(zQ)
        zR_prev = deepcopy(zR)
        zT_prev = deepcopy(zT)

        # V^2 por bus
        for b in shared_buses
            acc = 0.0
            n = 0
            for (_, sub) in partitions
                if b in sub["bus"]["bus_i"]
                    acc += value(sub["model"][:v_magnitude_sq][b, 1])
                    n += 1
                end
            end
            if n > 0
                v2_avg[b] = acc / n
            end
        end

        # P,Q,R,T por aresta
        for e in cut_edges
            accP = 0.0
            accQ = 0.0
            accR = 0.0
            accT = 0.0
            n = 0
            for (_, sub) in partitions
                l_local = findfirst(==(e), sub["branch_local_indices"])
                if l_local !== nothing
                    P_loc, Q_loc = _socp_flow_from_exprs_existing(sub["model"], l_local)
                    R_loc, T_loc = _rt_from_exprs_socp(sub["model"], l_local)
                    accP += value(P_loc)
                    accQ += value(Q_loc)
                    accR += value(R_loc)
                    accT += value(T_loc)
                    n += 1
                end
            end
            if n > 0
                zP[e] = accP / n
                zQ[e] = accQ / n
                zR[e] = accR / n
                zT[e] = accT / n
            end
        end

        # === 3) atualizar duais ===
        for (zone_id, sub) in partitions
            m = sub["model"]

            # BUS (V^2)
            for b in shared_buses
                if b in sub["bus"]["bus_i"]
                    V2loc = value(m[:v_magnitude_sq][b, 1])
                    duals_v2[zone_id][b] =
                        get(duals_v2[zone_id], b, 0.0) + ρV2 * (V2loc - v2_avg[b])
                end
            end

            # BRANCH (P,Q,R,T)
            for e in _zone_cut_edges(sub, cut_edges)
                l_local = findfirst(==(e), sub["branch_local_indices"])
                P_loc, Q_loc = _socp_flow_from_exprs_existing(m, l_local)
                R_loc, T_loc = _rt_from_exprs_socp(m, l_local)
                dualsP[zone_id][e] =
                    get(dualsP[zone_id], e, 0.0) + ρP * (value(P_loc) - zP[e])
                dualsQ[zone_id][e] =
                    get(dualsQ[zone_id], e, 0.0) + ρQ * (value(Q_loc) - zQ[e])
                dualsR[zone_id][e] =
                    get(dualsR[zone_id], e, 0.0) + ρR * (value(R_loc) - zR[e])
                dualsT[zone_id][e] =
                    get(dualsT[zone_id], e, 0.0) + ρT * (value(T_loc) - zT[e])
            end
        end

        # === 4) resíduos ===
        r_bus_sq = 0.0
        r_branch_sq = 0.0
        for (_, sub) in partitions
            m = sub["model"]
            # BUS
            for b in shared_buses
                if b in sub["bus"]["bus_i"]
                    r_bus_sq += (value(m[:v_magnitude_sq][b, 1]) - v2_avg[b])^2
                end
            end
            # BRANCH
            for e in _zone_cut_edges(sub, cut_edges)
                l_local = findfirst(==(e), sub["branch_local_indices"])
                P_loc, Q_loc = _socp_flow_from_exprs_existing(m, l_local)
                R_loc, T_loc = _rt_from_exprs_socp(m, l_local)
                r_branch_sq +=
                    (value(P_loc) - zP[e])^2 + (value(Q_loc) - zQ[e])^2 +
                    (value(R_loc) - zR[e])^2 + (value(T_loc) - zT[e])^2
            end
        end
        r_pri = sqrt(r_bus_sq + r_branch_sq)

        s_bus = sqrt(sum((v2_avg[b] - get(v2_prev, b, 0.0))^2 for b in shared_buses))
        s_branch = sqrt(
            sum(
                (zP[e] - get(zP_prev, e, 0.0))^2 + (zQ[e] - get(zQ_prev, e, 0.0))^2 +
                (zR[e] - get(zR_prev, e, 0.0))^2 + (zT[e] - get(zT_prev, e, 0.0))^2
                for e in cut_edges
            ),
        )
        s_tot = sqrt(s_bus^2 + s_branch^2)

        push!(res_pri_hist, r_pri)
        push!(res_du_hist, s_tot)

        # dimensões (para ε_pri e ε_du)
        dim_pri = 0
        for (_, sub) in partitions
            dim_pri += sum(b in sub["bus"]["bus_i"] for b in shared_buses)   # V^2
            dim_pri += 4 * length(_zone_cut_edges(sub, cut_edges))           # P,Q,R,T
        end
        dim_du = length(shared_buses) + 4 * length(cut_edges)

        ε_pri =
            sqrt(max(1, dim_pri)) * ϵ_abs +
            ϵ_rel * max(
                sqrt(
                    sum(v2_avg[b]^2 for b in shared_buses) +
                    sum(zP[e]^2 + zQ[e]^2 + zR[e]^2 + zT[e]^2 for e in cut_edges),
                ),
                1e-12,
            )

        ε_du =
            sqrt(max(1, dim_du)) * ϵ_abs +
            ϵ_rel * sqrt(
                sum((v2_avg[b] - get(v2_prev, b, 0.0))^2 for b in shared_buses) +
                sum(
                    (zP[e] - get(zP_prev, e, 0.0))^2 + (zQ[e] - get(zQ_prev, e, 0.0))^2 +
                    (zR[e] - get(zR_prev, e, 0.0))^2 + (zT[e] - get(zT_prev, e, 0.0))^2
                    for e in cut_edges
                ),
            )

        push!(eps_pri_hist, ε_pri)
        push!(eps_du_hist, ε_du)

        println("‖r‖ = $r_pri | ε_pri = $ε_pri")
        println("‖s‖ = $s_tot | ε_du  = $ε_du")

        push!(rhoV2_hist, ρV2)
        push!(rhoP_hist, ρP)
        push!(rhoQ_hist, ρQ)
        push!(rhoR_hist, ρR)
        push!(rhoT_hist, ρT)

        # parada
        if r_pri < ε_pri && s_tot < ε_du
            println("Stopping criterion met.")
            break
        end

        # === 5) atualização adaptativa de ρ (por bloco) ===
        # BUS (V^2)
        if r_bus_sq > (μ_incr^2) * (s_bus^2) && s_bus > 0.0
            ρV2 *= τ_incr
        elseif s_bus^2 > (μ_decr^2) * r_bus_sq && r_bus_sq > 0.0
            ρV2 *= τ_decr
        end
        # BRANCH (P,Q,R,T)
        if r_branch_sq > (μ_incr^2) * (s_branch^2) && s_branch > 0.0
            ρP *= τ_incr
            ρQ *= τ_incr
            ρR *= τ_incr
            ρT *= τ_incr
        elseif s_branch^2 > (μ_decr^2) * r_branch_sq && r_branch_sq > 0.0
            ρP *= τ_decr
            ρQ *= τ_decr
            ρR *= τ_decr
            ρT *= τ_decr
        end
    end

    # === 6) plots e resumo ===
    plot_admm_convergence(res_pri_hist, res_du_hist, eps_pri_hist, eps_du_hist,
        "plot_convergencia_admm_socp_hybrid.html")

    # rhos
    PlotlyLight.save(
        (
            () -> begin
                p = Plot()
                p(;
                    x = 1:length(rhoV2_hist),
                    y = rhoV2_hist,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρ(V²)",
                )
                p.layout.title = "Evolução de ρ(V²)"
                p.layout.xaxis.title = "Iterações"
                p.layout.yaxis.title = "ρ(V²)"
                p.layout.yaxis.type = "log"
                p
            end
        )(), "plot_rho_v2_socp_hybrid.html")

    PlotlyLight.save(
        (
            () -> begin
                p = Plot()
                p(;
                    x = 1:length(rhoP_hist),
                    y = rhoP_hist,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρP",
                )
                p.layout.title = "Evolução de ρP"
                p.layout.xaxis.title = "Iterações"
                p.layout.yaxis.title = "ρP"
                p.layout.yaxis.type = "log"
                p
            end
        )(), "plot_rhoP_socp_hybrid.html")
    PlotlyLight.save(
        (
            () -> begin
                p = Plot()
                p(;
                    x = 1:length(rhoQ_hist),
                    y = rhoQ_hist,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρQ",
                )
                p.layout.title = "Evolução de ρQ"
                p.layout.xaxis.title = "Iterações"
                p.layout.yaxis.title = "ρQ"
                p.layout.yaxis.type = "log"
                p
            end
        )(), "plot_rhoQ_socp_hybrid.html")
    PlotlyLight.save(
        (
            () -> begin
                p = Plot()
                p(;
                    x = 1:length(rhoR_hist),
                    y = rhoR_hist,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρR",
                )
                p.layout.title = "Evolução de ρR"
                p.layout.xaxis.title = "Iterações"
                p.layout.yaxis.title = "ρR"
                p.layout.yaxis.type = "log"
                p
            end
        )(), "plot_rhoR_socp_hybrid.html")
    PlotlyLight.save(
        (
            () -> begin
                p = Plot()
                p(;
                    x = 1:length(rhoT_hist),
                    y = rhoT_hist,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρT",
                )
                p.layout.title = "Evolução de ρT"
                p.layout.xaxis.title = "Iterações"
                p.layout.yaxis.title = "ρT"
                p.layout.yaxis.type = "log"
                p
            end
        )(), "plot_rhoT_socp_hybrid.html")

    total_real_cost = sum(p["real_generation_cost"] for p in values(partitions))
    total_expr_cost = sum(p["generation_cost_model_expr"] for p in values(partitions))

    println("\n[ADMM-SOCP-hybrid] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println(
        "ρ(V²)_final = $ρV2 | ρP_final = $ρP | ρQ_final = $ρQ | ρR_final = $ρR | ρT_final = $ρT",
    )

    return total_real_cost, partitions
end
