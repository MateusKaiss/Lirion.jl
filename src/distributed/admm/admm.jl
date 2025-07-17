function run_admm(config::DistributedOPFModel, data)
    partitions = prepare_admm_partitions(data)
    shared_buses = find_shared_buses(data, partitions)
    base_model = config.base_model

    duals = Dict{Int, Dict{Int, Float64}}()
    θ_avg = Dict{Int, Float64}()
    for zone_id in keys(partitions)
        duals[zone_id] = Dict{Int, Float64}()
    end

    for k = 1:config.max_iter
        println("Iteração $k")

        # --- Passo 1: resolve os subproblemas locais ---
        for (zone_id, subdata) in partitions
            local_model = build_model(subdata, base_model; optimizer = Ipopt.Optimizer)

            if base_model isa DCModel
                bus_ids = subdata["bus"]["bus_i"]
                @expression(local_model, admm_penalty, 0.0)

                for b in shared_buses
                    if b in bus_ids
                        idx = findfirst(isequal(b), bus_ids)
                        θ = local_model[:θ][idx]
                        θ̄ = get(θ_avg, b, 0.0)
                        λ = get(duals[zone_id], b, 0.0)
                        admm_penalty += λ * (θ - θ̄) + (config.rho / 2) * (θ - θ̄)^2
                    end
                end

                existing_objective = objective_function(local_model)
                @objective(local_model, Min, existing_objective + admm_penalty)
            else
                error("Base model $(base_model) ainda não implementado no ADMM")
            end

            optimize!(local_model)
            subdata[:model] = local_model
        end

        # --- Passo 2: atualiza variáveis compartilhadas ---
        if base_model isa DCModel
            update_shared_variables_dc!(
                shared_buses, partitions, duals, θ_avg, config.rho
            )
        end
    end

    return partitions
end