function run_admm(config::DistributedOPFModel, data)
    partitions = prepare_admm_partitions(data)
    shared_buses = find_shared_buses(data, partitions)
    base_model = config.base_model

    duals = Dict{Int, Dict{Int, Float64}}()
    θ_avg = Dict{Int, Float64}()
    θ_avg_prev = Dict{Int, Float64}()

    for zone_id in keys(partitions)
        duals[zone_id] = Dict{Int, Float64}()
    end

    for k in 1:config.max_iter
        println("Iteração $k")

        # --- Passo 1: resolve os subproblemas locais ---
        for (zone_id, subdata) in partitions
            local_model = build_model(subdata, base_model; optimizer = Ipopt.Optimizer)
            if base_model == DC
                bus_ids = subdata["bus"]["bus_i"]
                @expression(local_model, admm_penalty, 0.0)

                for b in shared_buses
                    if b in bus_ids
                        θ = local_model[:theta][b,1]
                        theta_m= get(θ_avg, b, 0.0)
    

                        λ = get(duals[zone_id], b, 0.0)
    
                        admm_penalty += λ * (θ - theta_m) + (config.rho / 2) * (θ - theta_m)^2
                    end
                end

                @objective(local_model, Min, local_model[:generation_cost] + admm_penalty)

            else
                error("Base model $(base_model) ainda não implementado no ADMM")
            end

            optimize!(local_model)

            # Salvar modelo resolvido
            subdata["model"] = local_model

            # Novo: salvar custo real (sem penalidade)
            real_cost = value(local_model[:generation_cost])
            subdata["real_generation_cost"] = real_cost
            println("→ Custo de geração da zona $zone_id: $real_cost")


            subdata["model"] = local_model
        end

        if base_model == DC
            update_shared_variables_dc!(shared_buses, partitions, duals, θ_avg, config.rho)
        end
    end

    total_objective = 0.0
    for (_, subdata) in partitions
        local_model = subdata["model"]
        real_cost = subdata["real_generation_cost"]
        total_objective += real_cost

    end


    println("→ Valor total da função objetivo: $total_objective")
    return total_objective, partitions

end
