function run_admm(config::DistributedOPFModel, data)
    partitions, shared_buses = prepare_admm_partitions(data)
    base_model = config.base_model

    # Initialize local rho variable from config (don't mutate config.rho)
    current_rho = config.rho

    duals = Dict{Int, Dict{Int, Float64}}()
    θ_avg = Dict{Int, Float64}()
    θ_avg_prev = Dict{Int, Float64}()
    res_primal_vec = Float64[]
    res_dual_vec = Float64[]
    eps_primal_vec = Float64[]
    eps_dual_vec = Float64[]
    rho_vec = Float64[]
    θ_avg_history = Dict{Int, Vector{Float64}}()
    θ_local_history = Dict{Int, Dict{Int, Vector{Float64}}}()  # zona → barra → histórico

    for zone_id in keys(partitions)
        duals[zone_id] = Dict{Int, Float64}()
    end

    μ_incr, μ_decr = 1.15, 1.15
    τ_incr, τ_decr = 1.15, 0.9
    ϵ_abs, ϵ_rel = 1e-4, 1e-4

    for k in 1:config.max_iter
        println("\nIteração $k")

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
                        admm_penalty += λ * (θ - θ̄) + (current_rho / 2) * (θ - θ̄)^2
                    end
                end
                @objective(local_model, Min,
                    local_model[:generation_cost] + local_model[:virtual_penalty] + admm_penalty
                )

            else
                error("Base model $(base_model) not yet implemented in ADMM")
            end

            optimize!(local_model)
            subdata["model"] = local_model
            ### TODO: review whether local_model[:generation_cost]:
            ### The local ADMM subproblem's generation cost is not the same as the real, true generation cost.
            ### Needs reviewing whether [:generation_cost] is the correct variable to use, or re-computing the 
            ### generation cost from the local model at the end of the ADMM execution is more appropriate.
            ### (I used to do the latter and it works as expected.)
            ### This conclusion comes after reducing the ϵ_abs and ϵ_rel to 1e-6, 1e-6, and having much smaller 
            ### ADMM convergence residuals of 498e-12 and 107e-12 for ||r_k|| and ||s_k||, respectively, but
            ### resulting in worse "Total objective value" output.
            subdata["real_generation_cost"] = value(local_model[:generation_cost])
            println("Zone generation cost $zone_id: $(subdata["real_generation_cost"])")

            if base_model == DC
                model = subdata["model"]
                θ_vals = value.(model[:theta])
                bus_ids = subdata["bus"]["bus_i"]

                ### θ_local_history is only recorded for shared buses, not for local buses
                ### TODO: add local buses to θ_local_history and ensure consistency with input file of 57 buses
                for b in shared_buses
                    if b in bus_ids
                        θ_b = θ_vals[b, 1]
                        push!(get!(get!(θ_local_history, zone_id, Dict{Int, Vector{Float64}}()), b, Float64[]), θ_b)
                    end
                end
            end
        end

        if base_model == DC
            θ_avg_prev = deepcopy(θ_avg)
            update_shared_variables_dc!(shared_buses, partitions, duals, θ_avg, current_rho)
            
            ### TODO: ensure consistency with input file of 57 buses and actual shared buses. 
            ### There seems to be quite a few shared buses that are not reflected in the output html,
            ### which indicated that θ_avg_history is not consistent with the input file.
            ### Perhaps a processing issue?
            for b in shared_buses
                push!(get!(θ_avg_history, b, Float64[]), θ_avg[b])
            end
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

        # Track current_rho instead of config.rho
        push!(rho_vec, current_rho)

        if res_primal < ϵ_primal && res_dual < ϵ_dual
            println("Stopping criterion met.")
            break
        end

        # Adaptive rho update: modify current_rho according to [...]
        ### TODO: add publication reference here
        if res_primal > μ_incr * res_dual
            current_rho *= τ_incr
        elseif res_dual > μ_decr * res_primal
            current_rho *= τ_decr
        end
    end

    plot_admm_convergence(
        res_primal_vec, res_dual_vec,
        eps_primal_vec, eps_dual_vec,
        "plot_convergencia_admm.html",
    )

    plot_rho(rho_vec, "plot_rho_admm.html")
    plot_theta_avg_heatmap(θ_avg_history, "plot_theta_avg_heatmap.html")
    plot_theta_local_heatmaps(θ_local_history, "plot_theta_local_heatmap")

    total_cost = sum(subdata["real_generation_cost"] for subdata in values(partitions))
    println("Total objective value: $total_cost")
    println("Valor inicial de rho: $(config.rho)")
    println("Valor final de rho: $current_rho")
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
    return println("Gráfico de convergência salvo em $path")
end

function plot_rho(rho_vec::Vector{Float64}, path::String)
    x_vals = 1:length(rho_vec)
    p = Plot()
    p(;
        x = x_vals,
        y = rho_vec,
        type = :scatter,
        mode = "lines+markers",
        name = "ρ",
        line = (; color = "green"),
    )

    p.layout.title = "Evolução de ρ (rho) no ADMM"
    p.layout.xaxis.title = "Iterações k"
    p.layout.yaxis.title = "ρ"
    p.layout.yaxis.type = "log"

    PlotlyLight.save(p, path)
    return println("Gráfico de rho salvo em $path")
end

function plot_theta_avg_heatmap(θ_avg_history::Dict{Int, Vector{Float64}}, path::String)
    buses = sort(collect(keys(θ_avg_history)))
    z_matrix = [θ_avg_history[b] for b in buses]

    p = Plot()
    p(;
        z = z_matrix,
        x = 1:length(z_matrix[1]),
        y = buses,
        type = :heatmap,
        colorscale = "Viridis",
        colorbar = (; title = "θ_avg"),
    )

    p.layout.title = "Evolução de θ_avg por Barra Compartilhada"
    p.layout.xaxis.title = "Iterações"
    p.layout.yaxis.title = "Ângulo de fase (rad) por Barra Compartilhada"

    PlotlyLight.save(p, path)
    return println("Heatmap salvo em $path")
end

function plot_theta_local_heatmaps(θ_local_history::Dict{Int, Dict{Int, Vector{Float64}}}, base_path::String)
    for (zone_id, zone_data) in θ_local_history
        buses = sort(collect(keys(zone_data)))
        z_matrix = [zone_data[b] for b in buses]

        p = Plot()
        p(
            z = z_matrix,
            x = 1:length(z_matrix[1]),
            y = buses,
            type = :heatmap,
            colorscale = "Viridis",
            colorbar = (; title = "θ (local)"),
        )

        p.layout.title = "Zona $zone_id – θ local por Barra Local"
        p.layout.xaxis.title = "Iterações"
        p.layout.yaxis.title = "Ângulo de fase (rad) por Barra Local"

        path = "$(base_path)_zone_$zone_id.html"
        PlotlyLight.save(p, path)
        println(" Heatmap local salvo em $path")
    end
end
