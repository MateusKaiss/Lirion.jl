function run_admm(config::DistributedOPFModel, data)
    partitions, shared_buses = prepare_admm_partitions(data)
    base_model = config.base_model

    duals = Dict{Int, Dict{Int, Float64}}()
    θ_avg = Dict{Int, Float64}()
    θ_avg_prev = Dict{Int, Float64}()
    res_primal_vec = Float64[]
    res_dual_vec = Float64[]

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
        push!(res_primal_vec, res_primal)
        push!(res_dual_vec, res_dual)


        println("Res Primal: $res_primal")
        println("Res Dual:   $res_dual")

        if res_primal < 1e-4 && res_dual < 1e-4
            println("Stopping criterion met.")
            break
        end
    end
    plot_admm_residuals(res_primal_vec, res_dual_vec, "residuos_admm.html")

    # --- Final cost calculation ---
    total_objective = sum(subdata["real_generation_cost"] for subdata in values(partitions))
    println("Total value of the objective function: $total_objective")
    return total_objective, partitions
end

function plot_admm_residuals(primal::Vector{Float64}, dual::Vector{Float64}, path::String)
    x_vals = 1:length(primal)
    p = Plot()
    p(; x=x_vals, y=primal, type=:scatter, mode="lines+markers", name="Primal")
    p(; x=x_vals, y=dual,   type=:scatter, mode="lines+markers", name="Dual")

    p.layout.title = "Convergência do ADMM"
    p.layout.xaxis.title = "Iteração"
    p.layout.yaxis.title = "Resíduo"
    p.layout.yaxis.type = "log"

    PlotlyLight.save(p, path)
    println(" Gráfico salvo em $path")
end
