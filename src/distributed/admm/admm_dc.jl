function run_admm(config::DistributedOPFModel, data)
    partitions, shared_buses = prepare_admm_partitions(data)
    base_model = config.base_model

    base_model == DC || error("run_admm (bus-consensus) está implementado para DC.")

    current_rho = config.rho_theta

    duals = Dict{Int, Dict{Int, Float64}}()
    θ_avg = Dict{Int, Float64}()
    θ_avg_prev = Dict{Int, Float64}()
    res_primal_vec = Float64[]
    res_dual_vec = Float64[]
    eps_primal_vec = Float64[]
    eps_dual_vec = Float64[]
    rho_vec = Float64[]
    θ_avg_history = Dict{Int, Vector{Float64}}()
    θ_local_history = Dict{Int, Dict{Int, Vector{Float64}}}()

    for zone_id in keys(partitions)
        duals[zone_id] = Dict{Int, Float64}()
    end

    μ_incr, μ_decr = 2.0, 2.0
    τ_incr, τ_decr = 1.05, 0.98
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual

    for k in 1:config.max_iter
        println("\nIteração $k")

        for (zone_id, subdata) in partitions
            local_model = build_model(subdata, base_model; optimizer = Ipopt.Optimizer)

            bus_ids = subdata["bus"]["bus_i"]
            @expression(local_model, admm_penalty, 0.0)
            for b in shared_buses
                if b in bus_ids
                    θ = local_model[:theta][b, 1]
                    θ̄ = get(θ_avg, b, 0.0)
                    λ = get(duals[zone_id], b, 0.0)
                    admm_penalty += λ * (θ - θ̄) + (current_rho / 2) * (θ - θ̄)^2
                end
            end

            @objective(local_model, Min, local_model[:generation_cost] + admm_penalty)

            optimize!(local_model)
            # print(local_model)
            subdata["model"] = local_model
            ### TODO: review whether local_model[:generation_cost]:
            ### The local ADMM subproblem's generation cost is not the same as the real, true generation cost.
            ### Needs reviewing whether [:generation_cost] is the correct variable to use, or re-computing the 
            ### generation cost from the local model at the end of the ADMM execution is more appropriate.
            ### (I used to do the latter and it works as expected.)
            ### This conclusion comes after reducing the ϵ_abs and ϵ_rel to 1e-6, 1e-6, and having much smaller 
            ### ADMM convergence residuals of 498e-12 and 107e-12 for ||r_k|| and ||s_k||, respectively, but
            ### resulting in worse "Total objective value" output.

            # Cost calculated directly from the model (with possible distortions)
            subdata["generation_cost_model_expr"] = value(local_model[:generation_cost])

            # Actual cost recalculated from variable values
            subdata["real_generation_cost"] =
                compute_real_generation_cost(local_model, subdata["objective_params"])
            println("Zone generation cost $zone_id: $(subdata["real_generation_cost"])")

            θ_vals = value.(local_model[:theta])
            for b in shared_buses
                if b in bus_ids
                    θ_b = θ_vals[b, 1]
                    push!(
                        get!(
                            get!(θ_local_history, zone_id, Dict{Int, Vector{Float64}}()),
                            b,
                            Float64[],
                        ),
                        θ_b,
                    )
                end
            end
        end

        # Updates the shared variables (e.g., θ_avg) and records their history
        θ_avg_prev = deepcopy(θ_avg)
        update_shared_variables_dc!(shared_buses, partitions, duals, θ_avg, current_rho)
        for b in shared_buses
            push!(get!(θ_avg_history, b, Float64[]), θ_avg[b])
        end

        res_primal = sqrt(
            sum(
                sum(
                    (value(partitions[z]["model"][:theta][b, 1]) - θ_avg[b])^2
                    for b in shared_buses if b in partitions[z]["bus"]["bus_i"]
                ) for z in keys(partitions)
            ),
        )

        res_dual = sqrt(sum((θ_avg[b] - get(θ_avg_prev, b, 0.0))^2 for b in shared_buses))
        push!(res_primal_vec, res_primal)
        push!(res_dual_vec, res_dual)

        N = length(shared_buses)
        ϵ_primal =
            sqrt(N) * ϵ_abs +
            ϵ_rel * max(
                sqrt(
                    sum(
                        (value(partitions[z]["model"][:theta][b, 1]))^2
                        for z in keys(partitions), b in shared_buses
                        if b in partitions[z]["bus"]["bus_i"]
                    ),
                ),
                sqrt(sum(θ_avg[b]^2 for b in shared_buses)),
            )
        ϵ_dual =
            sqrt(N) * ϵ_abs +
            ϵ_rel * sqrt(sum((θ_avg[b] - get(θ_avg_prev, b, 0.0))^2 for b in shared_buses))

        push!(eps_primal_vec, ϵ_primal)
        push!(eps_dual_vec, ϵ_dual)
        println("‖r‖ = $res_primal | ε_pri = $ϵ_primal")
        println("‖s‖ = $res_dual   | ε_dual = $ϵ_dual")

        push!(rho_vec, current_rho)

        if res_primal < ϵ_primal && res_dual < ϵ_dual
            println("Stopping criterion met.")
            break
        end

        # Adaptive rho update: modify current_rho according to [...]
        ### TODO: add publication reference here

        #=
        Adaptive update of the ADMM penalty parameter ρ based on the relative size of
        primal and dual residuals.

        This strategy, known as *residual balancing*, adjusts ρ to improve convergence:
        - If the primal residual is much larger than the dual, ρ is increased to penalize consensus violation.
        - If the dual residual dominates, ρ is reduced to allow better primal progress.

        This heuristic is commonly adopted in distributed OPF methods to accelerate convergence
        and improve robustness under ill-conditioned settings.

        @article{mhanna2019adaptive,
            title  = {Adaptive ADMM for Distributed AC Optimal Power Flow},
            author = {Mhanna, Sleiman and Verbi{\v{c}}, Gregor and Chapman, Archie C},
            journal = {IEEE Transactions on Power Systems},
            volume = {34},
            number = {3},
            pages  = {2025--2035},
            year   = {2019},
            publisher = {IEEE},
            doi    = {10.1109/TPWRS.2018.2886344}
        }

        @article{gebbran2021fair,
            title  = {Fair coordination of distributed energy resources with Volt-Var control and PV curtailment},
            author = {Gebbran, Daniel and Mhanna, Sleiman and Ma, Yiju and Chapman, Archie C and Verbi{\v{c}}, Gregor},
            journal = {Applied Energy},
            volume  = {286},
            pages   = {116546},
            year    = {2021},
            publisher = {Elsevier},
            doi     = {10.1016/j.apenergy.2021.116546}
        }
        =#

        if res_primal > μ_incr * res_dual
            current_rho *= τ_incr
        elseif res_dual > μ_decr * res_primal
            current_rho *= τ_decr
        end
    end

    plot_admm_convergence(res_primal_vec, res_dual_vec, eps_primal_vec, eps_dual_vec,
        "plot_convergencia_admm.html")
    plot_rho(rho_vec, "plot_rho_admm.html")
    plot_theta_avg_heatmap(θ_avg_history, "plot_theta_avg_heatmap.html")
    plot_theta_local_heatmaps(θ_local_history, "plot_theta_local_heatmap")

    total_real_cost = sum(subdata["real_generation_cost"] for subdata in values(partitions))
    total_expr_cost =
        sum(subdata["generation_cost_model_expr"] for subdata in values(partitions))

    println("\nResumo geral:")
    println("→ Total (value(:generation_cost)) = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("Valor inicial de ρθ: $(config.rho_theta)")
    println("Valor final de ρθ: $current_rho")

    return total_real_cost, partitions
end

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

function run_admm_hybrid(config::DistributedOPFModel, data::Dict{String, Any})
    @assert config.base_model == DC "Hybrid ADMM currently implemented for DC only."

    partitions, shared_buses = prepare_admm_partitions(data)
    cut_edges = _find_cut_edges(data)

    θ_avg = Dict{Int, Float64}()
    θ_avg_prev = Dict{Int, Float64}()
    z_flow = Dict{Int, Float64}()
    z_flow_prev = Dict{Int, Float64}()

    duals_theta = Dict{Int, Dict{Int, Float64}}()  # λθ[zone][bus]
    duals_flow = Dict{Int, Dict{Int, Float64}}()  # λe[zone][edge]
    for z in keys(partitions)
        duals_theta[z] = Dict{Int, Float64}()
        duals_flow[z] = Dict{Int, Float64}()
    end
    for e in cut_edges
        z_flow[e] = 0.0
    end

    ρθ = config.rho_theta
    ρe = config.rho_flow

    res_primal_vec = Float64[]   # híbrido (agregado)
    res_dual_vec = Float64[]   # híbrido (agregado)
    eps_primal_vec = Float64[]
    eps_dual_vec = Float64[]
    rhoθ_vec = Float64[]
    rhoe_vec = Float64[]

    μ_incr, μ_decr = 2.0, 2.0
    τ_incr, τ_decr = 1.05, 0.98
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual
    θ_avg_history = Dict{Int, Vector{Float64}}()

    for k in 1:config.max_iter
        println("\n[ADMM-hybrid] Iteração $k")
        for (zone_id, subdata) in partitions
            local_model = build_model(subdata, DC; optimizer = Ipopt.Optimizer)

            @expression(local_model, admm_penalty, 0.0)
            bus_ids = subdata["bus"]["bus_i"]
            for b in shared_buses
                if b in bus_ids
                    θ = local_model[:theta][b, 1]
                    θ̄ = get(θ_avg, b, 0.0)
                    λθ = get(duals_theta[zone_id], b, 0.0)
                    admm_penalty += λθ * (θ - θ̄) + (ρθ / 2) * (θ - θ̄)^2
                end
            end

            zone_edges = _zone_cut_edges(subdata, cut_edges)
            for e in zone_edges
                l_local = findfirst(==(e), subdata["branch_local_indices"])
                p_loc = _dc_flow_expr(local_model, subdata, l_local)
                λe_loc = get(duals_flow[zone_id], e, 0.0)
                z = z_flow[e]
                admm_penalty += λe_loc * (p_loc - z) + (ρe / 2) * (p_loc - z)^2
            end

            @objective(local_model, Min, local_model[:generation_cost] + admm_penalty)
            optimize!(local_model)

            subdata["model"] = local_model
            subdata["generation_cost_model_expr"] = value(local_model[:generation_cost])
            subdata["real_generation_cost"] =
                compute_real_generation_cost(local_model, subdata["objective_params"])
        end

        θ_avg_prev = isempty(θ_avg) ? θ_avg : deepcopy(θ_avg)
        z_flow_prev = deepcopy(z_flow)

        for b in shared_buses
            acc = 0.0
            num = 0
            for (_, subdata) in partitions
                if b in subdata["bus"]["bus_i"]
                    θb = value(subdata["model"][:theta][b, 1])
                    acc += θb
                    num += 1
                end
            end
            if num > 0
                θ_avg[b] = acc / num
                push!(get!(θ_avg_history, b, Float64[]), θ_avg[b])
            else
                θ_avg[b] = get(θ_avg, b, 0.0)
            end
        end

        for e in cut_edges
            acc = 0.0
            num = 0
            for (_, subdata) in partitions
                local_idx = findfirst(==(e), subdata["branch_local_indices"])
                if local_idx !== nothing
                    p_val = value(_dc_flow_expr(subdata["model"], subdata, local_idx))
                    acc += p_val
                    num += 1
                end
            end
            z_flow[e] = (num > 0) ? acc / num : get(z_flow, e, 0.0)
        end

        for (zone_id, subdata) in partitions
            for b in shared_buses
                if b in subdata["bus"]["bus_i"]
                    θb = value(subdata["model"][:theta][b, 1])
                    duals_theta[zone_id][b] =
                        get(duals_theta[zone_id], b, 0.0) + ρθ * (θb - θ_avg[b])
                end
            end
            for e in _zone_cut_edges(subdata, cut_edges)
                l_local = findfirst(==(e), subdata["branch_local_indices"])
                p_val = value(_dc_flow_expr(subdata["model"], subdata, l_local))
                duals_flow[zone_id][e] =
                    get(duals_flow[zone_id], e, 0.0) + ρe * (p_val - z_flow[e])
            end
        end

        # r_pri^2 = sum_b,sum_z (θ_zb - θ̄_b)^2  +  sum_e,sum_z (p_ze - z_e)^2
        res_primal_sq = 0.0
        for (_, subdata) in partitions
            # ângulos
            for b in shared_buses
                if b in subdata["bus"]["bus_i"]
                    θb = value(subdata["model"][:theta][b, 1])
                    res_primal_sq += (θb - θ_avg[b])^2
                end
            end
            # fluxos
            for e in _zone_cut_edges(subdata, cut_edges)
                l_local = findfirst(==(e), subdata["branch_local_indices"])
                p_val = value(_dc_flow_expr(subdata["model"], subdata, l_local))
                res_primal_sq += (p_val - z_flow[e])^2
            end
        end
        res_primal = sqrt(res_primal_sq)

        # r_dual^2 = sum_b (θ̄_b - θ̄_b_prev)^2 + sum_e (z_e - z_e_prev)^2
        res_dual = sqrt(
            sum((θ_avg[b] - get(θ_avg_prev, b, 0.0))^2 for b in shared_buses) +
            sum((z_flow[e] - get(z_flow_prev, e, 0.0))^2 for e in cut_edges),
        )

        push!(res_primal_vec, res_primal)
        push!(res_dual_vec, res_dual)

        dim_pri = 0
        for (_, subdata) in partitions
            dim_pri += sum(b in subdata["bus"]["bus_i"] for b in shared_buses)  # θ
            dim_pri += length(_zone_cut_edges(subdata, cut_edges))               # fluxos
        end
        dim_du = length(shared_buses) + length(cut_edges)

        ϵ_primal =
            sqrt(max(1, dim_pri)) * ϵ_abs +
            ϵ_rel * max(
                sqrt(
                    sum(θ_avg[b]^2 for b in shared_buses) +
                    sum(z_flow[e]^2 for e in cut_edges),
                ),
                1e-12,
            )

        ϵ_dual =
            sqrt(max(1, dim_du)) * ϵ_abs +
            ϵ_rel * sqrt(
                sum((θ_avg[b] - get(θ_avg_prev, b, 0.0))^2 for b in shared_buses) +
                sum((z_flow[e] - get(z_flow_prev, e, 0.0))^2 for e in cut_edges),
            )

        push!(eps_primal_vec, ϵ_primal)
        push!(eps_dual_vec, ϵ_dual)

        println("‖r‖ = $res_primal | ε_pri = $ϵ_primal")
        println("‖s‖ = $res_dual   | ε_dual = $ϵ_dual")

        push!(rhoθ_vec, ρθ)
        push!(rhoe_vec, ρe)

        if res_primal < ϵ_primal && res_dual < ϵ_dual
            println("Stopping criterion met.")
            break
        end

        if res_primal > μ_incr * res_dual
            ρθ *= τ_incr
            ρe *= τ_incr
        elseif res_dual > μ_decr * res_primal
            ρθ *= τ_decr
            ρe *= τ_decr
        end
    end

    plot_admm_convergence(res_primal_vec, res_dual_vec, eps_primal_vec, eps_dual_vec,
        "plot_convergencia_admm_hybrid.html")

    PlotlyLight.save(
        (
            () -> begin
                p = Plot()
                p(;
                    x = 1:length(rhoθ_vec),
                    y = rhoθ_vec,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρθ",
                )
                p.layout.title = "Evolução de ρθ"
                p.layout.xaxis.title = "Iterações"
                p.layout.yaxis.title = "ρθ"
                p.layout.yaxis.type = "log"
                p
            end
        )(), "plot_rho_theta_hybrid.html")

    PlotlyLight.save(
        (
            () -> begin
                p = Plot()
                p(;
                    x = 1:length(rhoe_vec),
                    y = rhoe_vec,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρe",
                )
                p.layout.title = "Evolução de ρe"
                p.layout.xaxis.title = "Iterações"
                p.layout.yaxis.title = "ρe"
                p.layout.yaxis.type = "log"
                p
            end
        )(), "plot_rho_flow_hybrid.html")

    plot_theta_avg_heatmap(θ_avg_history, "plot_theta_avg_heatmap_hybrid.html")

    total_real_cost = sum(p["real_generation_cost"] for p in values(partitions))
    total_expr_cost = sum(p["generation_cost_model_expr"] for p in values(partitions))

    println("\n[ADMM-hybrid] Resumo:")
    println("→ Total (value(:generation_cost))       = $total_expr_cost")
    println("→ Total (reavaliado via p_generation)   = $total_real_cost")
    println("ρθ_final = $ρθ | ρe_final = $ρe")

    return total_real_cost, partitions
end
