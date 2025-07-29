function run_admm(config::DistributedOPFModel, data)
    partitions, shared_buses = prepare_admm_partitions(data)
    base_model = config.base_model

    duals = Dict{Int, Dict{Int, Float64}}()
    θ_avg = Dict{Int, Float64}()
    θ_avg_prev = Dict{Int, Float64}()
    res_primal_vec = Float64[]
    res_dual_vec = Float64[]
    eps_primal_vec = Float64[]
    eps_dual_vec = Float64[]

    for zone_id in keys(partitions)
        duals[zone_id] = Dict{Int, Float64}()
    end

    μ_incr, μ_decr = 1.15, 1.15
    τ_incr, τ_decr = 1.15, 0.9
    ϵ_abs, ϵ_rel = 1e-4, 1e-4

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
                        θ̄ = get(θ_avg, b, 0.0)
                        λ = get(duals[zone_id], b, 0.0)
                        admm_penalty += λ * (θ - θ̄) + (config.rho / 2) * (θ - θ̄)^2
                    end
                end
                @objective(local_model, Min, local_model[:generation_cost] + admm_penalty)
            else
                error("Base model $(base_model) not yet implemented in ADMM")
            end

            optimize!(local_model)
            subdata["model"] = local_model
            subdata["real_generation_cost"] = value(local_model[:generation_cost])
            println("Zone generation cost $zone_id: $(subdata["real_generation_cost"])")
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
                )
                for z in keys(partitions)
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
                        for z in keys(partitions), b in shared_buses if
                        b in partitions[z]["bus"]["bus_i"]
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

        if res_primal < 1e-4 && res_dual < 1e-4
            println("Stopping criterion met.")
            break
        end

        # --- ρ adaptive heuristic ---
        if res_primal > μ_incr * res_dual
            config.rho *= τ_incr
        elseif res_dual > μ_decr * res_primal
            config.rho *= τ_decr
        end
    end

    plot_admm_convergence(
        res_primal_vec, res_dual_vec,
        eps_primal_vec, eps_dual_vec,
        "convergencia_admm.html",
    )

    total_cost = sum(subdata["real_generation_cost"] for subdata in values(partitions))
    println("Total objective value: $total_cost")
    return total_cost, partitions
end

function plot_admm_convergence(
    res_primal::Vector{Float64},
    res_dual::Vector{Float64},
    eps_primal::Vector{Float64},
    eps_dual::Vector{Float64},
    path::String,
)
    x_vals = 1:length(res_primal)
    p = Plot()
    p(;
        x = x_vals,
        y = res_primal,
        type = :scatter,
        mode = "lines",
        name = "‖rᵏ‖",
        line = (; color = "blue"),
    )
    p(;
        x = x_vals,
        y = res_dual,
        type = :scatter,
        mode = "lines",
        name = "‖sᵏ‖",
        line = (; color = "red"),
    )
    p(;
        x = x_vals,
        y = eps_primal,
        type = :scatter,
        mode = "lines",
        name = "εᵖʳⁱ",
        line = (; color = "purple", dash = "dash"),
    )
    p(;
        x = x_vals,
        y = eps_dual,
        type = :scatter,
        mode = "lines",
        name = "εᵈᵘᵃˡ",
        line = (; color = "orange", dash = "dash"),
    )

    p.layout.title = "Convergência do ADMM"
    p.layout.xaxis.title = "Iterações k"
    p.layout.yaxis.title = "Resíduo"
    p.layout.yaxis.type = "log"

    PlotlyLight.save(p, path)
    return println("Graph saved in $path")
end
