@inline function _ac_flow_from_exprs_existing(
    model::Model,
    subdata::Dict{String, Any},
    l_local::Int;
    t::Int = 1,
)
    return model[:p_flow_branches][l_local, t], model[:q_flow_branches][l_local, t]
end

# ------------------------------------------------------------------------------
# ADMM-AC por branch (consenso em P_flow e Q_flow)
#   - penaliza P_from e Q_from (lado "from")
#   - consenso: média aritmética entre as cópias locais
# ------------------------------------------------------------------------------
function run_admm_ac_branch(config::DistributedOPFModel, data::Dict{String, Any})
    @assert config.base_model == AC "ADMM-AC-branch: base_model deve ser AC."

    partitions, _ = prepare_admm_partitions(data)
    cut_edges = _find_cut_edges(data)

    dualsP = Dict{Int, Dict{Int, Float64}}()
    dualsQ = Dict{Int, Dict{Int, Float64}}()
    for z in keys(partitions)
        dualsP[z] = Dict{Int, Float64}()
        dualsQ[z] = Dict{Int, Float64}()
    end

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

    res_primal_vec = Float64[]
    res_dual_vec = Float64[]
    eps_primal_vec = Float64[]
    eps_dual_vec = Float64[]
    rhoP_vec = Float64[]
    rhoQ_vec = Float64[]

    μ_incr, μ_decr = 2.0, 2.0
    τ_incr, τ_decr = 1.05, 0.98
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual

    for k in 1:config.max_iter
        println("\n[ADMM-AC-branch] Iteração $k")

        for (zone_id, subdata) in partitions
            local_model = build_model(subdata, AC; optimizer = Ipopt.Optimizer)

            zone_edges = _zone_cut_edges(subdata, cut_edges)

            @expression(local_model, admm_penalty, 0.0)
            for e in zone_edges
                l_local = findfirst(==(e), subdata["branch_local_indices"])
                @assert l_local !== nothing "aresta global $e não encontrada na partição $zone_id"

                P_loc, Q_loc = _ac_flow_from_exprs_existing(local_model, subdata, l_local)

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
                    P_loc, Q_loc = _ac_flow_from_exprs_existing(m, subdata, l_local)
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

        for (zone_id, subdata) in partitions
            for e in _zone_cut_edges(subdata, cut_edges)
                l_local = findfirst(==(e), subdata["branch_local_indices"])
                m = subdata["model"]
                P_loc, Q_loc = _ac_flow_from_exprs_existing(m, subdata, l_local)
                dualsP[zone_id][e] =
                    get(dualsP[zone_id], e, 0.0) + ρP * (value(P_loc) - zP[e])
                dualsQ[zone_id][e] =
                    get(dualsQ[zone_id], e, 0.0) + ρQ * (value(Q_loc) - zQ[e])
            end
        end

        res_primal_sq = 0.0
        for (_, subdata) in partitions
            for e in _zone_cut_edges(subdata, cut_edges)
                l_local = findfirst(==(e), subdata["branch_local_indices"])
                m = subdata["model"]
                P_loc, Q_loc = _ac_flow_from_exprs_existing(m, subdata, l_local)
                res_primal_sq += (value(P_loc) - zP[e])^2 + (value(Q_loc) - zQ[e])^2
            end
        end
        res_primal = sqrt(res_primal_sq)

        res_dual = sqrt(
            sum(
                (zP[e] - get(zP_prev, e, 0.0))^2 + (zQ[e] - get(zQ_prev, e, 0.0))^2 for
                e in cut_edges
            ),
        )

        push!(res_primal_vec, res_primal)
        push!(res_dual_vec, res_dual)

        dim_pri = 0
        for (_, subdata) in partitions
            dim_pri += 2 * length(_zone_cut_edges(subdata, cut_edges))
        end
        dim_du = 2 * length(cut_edges)

        ϵ_primal =
            sqrt(max(1, dim_pri)) * ϵ_abs +
            ϵ_rel * max(sqrt(sum(zP[e]^2 + zQ[e]^2 for e in cut_edges)), 1e-12)

        ϵ_dual =
            sqrt(max(1, dim_du)) * ϵ_abs +
            ϵ_rel * sqrt(
                sum(
                    (zP[e] - get(zP_prev, e, 0.0))^2 +
                    (zQ[e] - get(zQ_prev, e, 0.0))^2 for e in cut_edges
                ),
            )

        push!(eps_primal_vec, ϵ_primal)
        push!(eps_dual_vec, ϵ_dual)

        println("‖r‖ = $res_primal | ε_pri = $ϵ_primal")
        println("‖s‖ = $res_dual   | ε_dual = $ϵ_dual")

        push!(rhoP_vec, ρP)
        push!(rhoQ_vec, ρQ)

        if res_primal < ϵ_primal && res_dual < ϵ_dual
            println("Stopping criterion met.")
            break
        end

        if res_primal > μ_incr * res_dual
            ρP *= τ_incr
            ρQ *= τ_incr
        elseif res_dual > μ_decr * res_primal
            ρP *= τ_decr
            ρQ *= τ_decr
        end
    end

    plot_admm_convergence(res_primal_vec, res_dual_vec, eps_primal_vec, eps_dual_vec,
        "plot_convergencia_admm_ac_branch.html")

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
        )(), "plot_rhoP_ac_branch.html")

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
        )(), "plot_rhoQ_ac_branch.html")

    total_real_cost = sum(p["real_generation_cost"] for p in values(partitions))
    total_expr_cost = sum(p["generation_cost_model_expr"] for p in values(partitions))

    println("\n[ADMM-AC-branch] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("ρP_final = $ρP | ρQ_final = $ρQ")

    return total_real_cost, partitions
end
