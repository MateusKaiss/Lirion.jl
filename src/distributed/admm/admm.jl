function run_admm(config::DistributedOPFModel, data)
    partitions, interconnection_limits = prepare_admm_partitions(data)
    shared_buses = find_shared_buses(data, partitions)
    base_model = config.base_model

    duals = Dict{Int, Dict{Int, Float64}}()
    θ_avg = Dict{Int, Float64}()
    θ_avg_prev = Dict{Int, Float64}()

    for zone_id in keys(partitions)
        duals[zone_id] = Dict{Int, Float64}()
    end

    for k in 1:config.max_iter
        println("\nIteração $k")

        # --- Step 1: solve the local subproblems ---
        for (zone_id, subdata) in partitions
            local_model = build_model(subdata, base_model; optimizer = Ipopt.Optimizer)
            if base_model == DC
                bus_ids = subdata["bus"]["bus_i"]
                @expression(local_model, admm_penalty, 0.0)

                for b in shared_buses
                    if b in bus_ids
                        θ = local_model[:theta][b, 1]
                        theta_m = get(θ_avg, b, 0.0)
                        λ = get(duals[zone_id], b, 0.0)

                        admm_penalty +=
                            λ * (θ - theta_m) + (config.rho / 2) * (θ - theta_m)^2
                    end
                end

                @objective(local_model, Min, local_model[:generation_cost] + admm_penalty)
            else
                error("Base model $(base_model) not yet implemented in ADMM")
            end

            optimize!(local_model)

            subdata["model"] = local_model
            real_cost = value(local_model[:generation_cost])
            subdata["real_generation_cost"] = real_cost
            println("Zone generation cost $zone_id: $real_cost")
        end

        if base_model == DC
            θ_avg_prev = deepcopy(θ_avg)
            update_shared_variables_dc!(shared_buses, partitions, duals, θ_avg, config.rho)
        end

        res_primal = sqrt(
            sum(
                sum(
                    (value(partitions[z]["model"][:theta][b, 1]) - θ_avg[b])^2
                    for b in shared_buses if b in partitions[z]["bus"]["bus_i"]
                    ; init = 0.0)
                for z in keys(partitions)
                ; init = 0.0),
        )

        res_dual = sqrt(
            sum((θ_avg[b] - get(θ_avg_prev, b, 0.0))^2 for b in shared_buses; init = 0.0),
        )

        println("Res Primal: $res_primal")
        println("Res Dual:   $res_dual")

        if res_primal < 1e-4 && res_dual < 1e-4
            println("Stopping criterion met.")
            break
        end
    end

    # --- Final cost calculation ---
    total_objective = sum(subdata["real_generation_cost"] for subdata in values(partitions))
    println("Total value of the objective function: $total_objective")
    return total_objective, partitions
end
