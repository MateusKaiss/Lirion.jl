function run_admm_branch(config::DistributedOPFModel, data::Dict{String, Any})
    partitions, _shared_buses_unused = prepare_admm_partitions(data)
    base_model = config.base_model
    base_model == DC || error("ADMM-branch está implementado aqui apenas para DC.")
    cut_edges = _find_cut_edges(data)

    duals = Dict{Int, Dict{Int, Float64}}()
    z_flow = Dict{Int, Float64}()
    z_flow_prev = Dict{Int, Float64}()

    current_rho = config.rho_flow

    for zone_id in keys(partitions)
        duals[zone_id] = Dict{Int, Float64}()
    end
    for e in cut_edges
        z_flow[e] = 0.0
    end

    res_primal_vec = Float64[]
    res_dual_vec = Float64[]
    eps_primal_vec = Float64[]
    eps_dual_vec = Float64[]
    rho_vec = Float64[]

    μ_incr, μ_decr = 2.0, 2.0
    τ_incr, τ_decr = 1.05, 0.98
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual

    for k in 1:config.max_iter
        println("\n[ADMM-branch] Iteração $k")
        for (zone_id, subdata) in partitions
            local_model = build_model(subdata, base_model; optimizer = Ipopt.Optimizer)

            zone_edges = _zone_cut_edges(subdata, cut_edges)
            @expression(local_model, admm_penalty, 0.0)
            for e in zone_edges
                l_local = findfirst(==(e), subdata["branch_local_indices"])
                p_loc = _dc_flow_expr(local_model, subdata, l_local)
                λ = get(duals[zone_id], e, 0.0)
                z = z_flow[e]
                admm_penalty += λ * (p_loc - z) + (current_rho / 2) * (p_loc - z)^2
            end

            @objective(local_model, Min, local_model[:generation_cost] + admm_penalty)

            optimize!(local_model)
            subdata["model"] = local_model
            subdata["generation_cost_model_expr"] = value(local_model[:generation_cost])
            subdata["real_generation_cost"] =
                compute_real_generation_cost(local_model, subdata["objective_params"])
        end

        z_flow_prev = deepcopy(z_flow)

        for e in cut_edges
            num = 0
            acc = 0.0
            for (_, subdata) in partitions
                local_idx = findfirst(==(e), subdata["branch_local_indices"])
                if local_idx !== nothing
                    m = subdata["model"]
                    p_val = value(_dc_flow_expr(m, subdata, local_idx))
                    acc += p_val
                    num += 1
                end
            end
            z_flow[e] = (num > 0) ? acc / num : 0.0
        end

        for (zone_id, subdata) in partitions
            for e in _zone_cut_edges(subdata, cut_edges)
                l_local = findfirst(==(e), subdata["branch_local_indices"])
                m = subdata["model"]
                p_val = value(_dc_flow_expr(m, subdata, l_local))
                λ = get(duals[zone_id], e, 0.0)
                duals[zone_id][e] = λ + current_rho * (p_val - z_flow[e])
            end
        end

        res_primal = 0.0
        for (_, subdata) in partitions
            for e in _zone_cut_edges(subdata, cut_edges)
                l_local = findfirst(==(e), subdata["branch_local_indices"])
                m = subdata["model"]
                p_val = value(_dc_flow_expr(m, subdata, l_local))
                res_primal += (p_val - z_flow[e])^2
            end
        end
        res_primal = sqrt(res_primal)

        res_dual = sqrt(sum((z_flow[e] - get(z_flow_prev, e, 0.0))^2 for e in cut_edges))

        push!(res_primal_vec, res_primal)
        push!(res_dual_vec, res_dual)
        dim_pri = sum(length(_zone_cut_edges(sub, cut_edges)) for (_, sub) in partitions)
        dim_du = length(cut_edges)

        ϵ_primal =
            sqrt(max(1, dim_pri)) * ϵ_abs +
            ϵ_rel * max(sqrt(sum(z_flow[e]^2 for e in cut_edges)), 1e-12)
        ϵ_dual =
            sqrt(max(1, dim_du)) * ϵ_abs +
            ϵ_rel * sqrt(sum((z_flow[e] - get(z_flow_prev, e, 0.0))^2 for e in cut_edges))

        push!(eps_primal_vec, ϵ_primal)
        push!(eps_dual_vec, ϵ_dual)

        println("‖r‖ = $res_primal | ε_pri = $ϵ_primal")
        println("‖s‖ = $res_dual   | ε_dual = $ϵ_dual")

        push!(rho_vec, current_rho)

        if res_primal < ϵ_primal && res_dual < ϵ_dual
            println("Stopping criterion met.")
            break
        end

        if res_primal > μ_incr * res_dual
            current_rho *= τ_incr
        elseif res_dual > μ_decr * res_primal
            current_rho *= τ_decr
        end
    end

    plot_admm_convergence(res_primal_vec, res_dual_vec, eps_primal_vec, eps_dual_vec,
        "plot_convergencia_admm_branch.html")
    plot_rho(rho_vec, "plot_rho_admm_branch.html")

    total_real_cost = sum(p["real_generation_cost"] for p in values(partitions))
    total_expr_cost = sum(p["generation_cost_model_expr"] for p in values(partitions))

    println("\n[ADMM-branch] Resumo:")
    println("→ Total (value(:generation_cost)) = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("ρe inicial: $(config.rho_flow) | ρe final: $current_rho")

    return total_real_cost, partitions
end
