@inline function _socp_flow_from_exprs_existing(
    model::Model, l_local::Int; t::Int = 1,
)
    return model[:p_flow_branches][l_local, t], model[:q_flow_branches][l_local, t]
end

function run_admm_socp_branch(config::DistributedOPFModel, data::Dict{String, Any})
    @assert config.base_model == SOCP "ADMM-SOCP-branch: base_model deve ser SOCP."

    partitions, _ = prepare_admm_partitions(data)
    cut_edges = _find_cut_edges(data)

    # duais por zona/aresta
    dualsP = Dict{Int, Dict{Int, Float64}}()
    dualsQ = Dict{Int, Dict{Int, Float64}}()
    for z in keys(partitions)
        dualsP[z] = Dict{Int, Float64}()
        dualsQ[z] = Dict{Int, Float64}()
    end

    # variáveis de consenso por aresta de corte
    zP = Dict{Int, Float64}()  # ativo
    zQ = Dict{Int, Float64}()  # reativo
    for e in cut_edges
        zP[e] = 0.0
        zQ[e] = 0.0
    end
    zP_prev = Dict{Int, Float64}()
    zQ_prev = Dict{Int, Float64}()

    ρP = config.rho_flow_p
    ρQ = config.rho_flow_q
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual

    res_primal_vec = Float64[]
    res_dual_vec = Float64[]
    eps_primal_vec = Float64[]
    eps_dual_vec = Float64[]
    rhoP_vec = Float64[]
    rhoQ_vec = Float64[]

    μ_incr, μ_decr = 10.0, 10.0
    τ_incr, τ_decr = 2.0, 0.5

    for k in 1:config.max_iter
        println("\n[ADMM-SOCP-branch] Iteração $k")

        # 1) resolve subproblemas locais (SOCP)
        for (zone_id, subdata) in partitions
            local_model = build_model(subdata, SOCP; optimizer = SCS.Optimizer)
            # tolerâncias do SCS
            set_optimizer_attribute(local_model, "eps_abs", 1e-5)
            set_optimizer_attribute(local_model, "eps_rel", 1e-5)
            set_optimizer_attribute(local_model, "max_iters", 150_000)

            # aceleradores do SCS:
            set_optimizer_attribute(local_model, "acceleration_lookback", 10)
            set_optimizer_attribute(local_model, "acceleration_interval", 10)
            set_optimizer_attribute(local_model, "alpha", 1.8)
            set_optimizer_attribute(local_model, "normalize", true)
            set_optimizer_attribute(local_model, "scale", 1.0)

            zone_edges = _zone_cut_edges(subdata, cut_edges)

            @expression(local_model, admm_penalty, 0.0)
            for e in zone_edges
                l_local = findfirst(==(e), subdata["branch_local_indices"])
                @assert l_local !== nothing "aresta global $e não encontrada na partição $zone_id"
                P_loc, Q_loc = _socp_flow_from_exprs_existing(local_model, l_local)

                λP = get(dualsP[zone_id], e, 0.0)
                λQ = get(dualsQ[zone_id], e, 0.0)

                admm_penalty += λP * (P_loc - zP[e]) + (ρP / 2) * (P_loc - zP[e])^2
                admm_penalty += λQ * (Q_loc - zQ[e]) + (ρQ / 2) * (Q_loc - zQ[e])^2
            end

            @objective(local_model, Min, local_model[:generation_cost] + admm_penalty)
            optimize!(local_model)

            subdata["model"] = local_model
            subdata["generation_cost_model_expr"] = value(local_model[:generation_cost])
            subdata["real_generation_cost"] =
                compute_real_generation_cost(local_model, subdata["objective_params"])
        end

        # 2) atualiza consensos zP/zQ (média das cópias locais)
        zP_prev = deepcopy(zP)
        zQ_prev = deepcopy(zQ)
        for e in cut_edges
            accP = 0.0
            accQ = 0.0
            num = 0
            for (_, subdata) in partitions
                l_local = findfirst(==(e), subdata["branch_local_indices"])
                if l_local !== nothing
                    m = subdata["model"]
                    P_loc, Q_loc = _socp_flow_from_exprs_existing(m, l_local)
                    accP += value(P_loc)
                    accQ += value(Q_loc)
                    num += 1
                end
            end
            if num > 0
                zP[e] = accP / num
                zQ[e] = accQ / num
            end
        end

        # 3) atualiza duais
        for (zone_id, subdata) in partitions
            for e in _zone_cut_edges(subdata, cut_edges)
                l_local = findfirst(==(e), subdata["branch_local_indices"])
                m = subdata["model"]
                P_loc, Q_loc = _socp_flow_from_exprs_existing(m, l_local)
                dualsP[zone_id][e] =
                    get(dualsP[zone_id], e, 0.0) + ρP * (value(P_loc) - zP[e])
                dualsQ[zone_id][e] =
                    get(dualsQ[zone_id], e, 0.0) + ρQ * (value(Q_loc) - zQ[e])
            end
        end

        # 4) resíduos
        res_primal_sq = 0.0
        for (_, subdata) in partitions
            for e in _zone_cut_edges(subdata, cut_edges)
                l_local = findfirst(==(e), subdata["branch_local_indices"])
                m = subdata["model"]
                P_loc, Q_loc = _socp_flow_from_exprs_existing(m, l_local)
                res_primal_sq += (value(P_loc) - zP[e])^2 + (value(Q_loc) - zQ[e])^2
            end
        end
        r_pri = sqrt(res_primal_sq)

        s_du = sqrt(
            sum(
                (zP[e] - get(zP_prev, e, 0.0))^2 + (zQ[e] - get(zQ_prev, e, 0.0))^2
                for e in cut_edges
            ),
        )

        push!(res_primal_vec, r_pri)
        push!(res_dual_vec, s_du)

        # dimensões para tolerâncias
        dim_pri = 0
        for (_, sub) in partitions
            dim_pri += 2 * length(_zone_cut_edges(sub, cut_edges))
        end
        dim_du = 2 * length(cut_edges)

        ε_pri =
            sqrt(max(1, dim_pri)) * ϵ_abs +
            ϵ_rel * max(sqrt(sum(zP[e]^2 + zQ[e]^2 for e in cut_edges)), 1e-12)

        ε_du =
            sqrt(max(1, dim_du)) * ϵ_abs +
            ϵ_rel * sqrt(
                sum(
                    (zP[e] - get(zP_prev, e, 0.0))^2 + (zQ[e] - get(zQ_prev, e, 0.0))^2
                    for e in cut_edges
                ),
            )

        push!(eps_primal_vec, ε_pri)
        push!(eps_dual_vec, ε_du)

        println("‖r‖ = $r_pri | ε_pri = $ε_pri")
        println("‖s‖ = $s_du  | ε_du  = $ε_du")

        push!(rhoP_vec, ρP)
        push!(rhoQ_vec, ρQ)

        # parada
        if r_pri < ε_pri && s_du < ε_du
            println("Stopping criterion met.")
            break
        end

        # 5) ajuste adaptativo
        if r_pri > μ_incr * s_du
            ρP *= τ_incr
            ρQ *= τ_incr
        elseif s_du > μ_decr * r_pri
            ρP *= τ_decr
            ρQ *= τ_decr
        end
    end

    # 6) plots e resumo
    plot_admm_convergence(res_primal_vec, res_dual_vec, eps_primal_vec, eps_dual_vec,
        "plot_convergencia_admm_socp_branch.html")

    PlotlyLight.save(
        (
            () -> begin
                p = Plot()
                p(;
                    x = 1:length(rhoP_vec),
                    y = rhoP_vec,
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
        )(), "plot_rhoP_socp_branch.html")

    PlotlyLight.save(
        (
            () -> begin
                p = Plot()
                p(;
                    x = 1:length(rhoQ_vec),
                    y = rhoQ_vec,
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
        )(), "plot_rhoQ_socp_branch.html")

    total_real_cost = sum(p["real_generation_cost"] for p in values(partitions))
    total_expr_cost = sum(p["generation_cost_model_expr"] for p in values(partitions))

    println("\n[ADMM-SOCP-branch] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("ρP_final = $ρP | ρQ_final = $ρQ")

    return total_real_cost, partitions
end
