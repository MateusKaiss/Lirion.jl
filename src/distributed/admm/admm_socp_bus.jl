
function run_admm_socp(config::DistributedOPFModel, data::Dict{String, Any})
    @assert config.base_model == SOCP "ADMM-SOCP-bus: base_model deve ser SOCP."

    partitions, shared_buses = prepare_admm_partitions(data)

    duals_theta = Dict{Int, Dict{Int, Float64}}()
    duals_vmag = Dict{Int, Dict{Int, Float64}}()
    for z in keys(partitions)
        duals_theta[z] = Dict{Int, Float64}()
        duals_vmag[z] = Dict{Int, Float64}()
    end

    θ̄ = Dict{Int, Float64}()    # média circular
    V̄ = Dict{Int, Float64}()    # média aritmética
    θ̄_prev = Dict{Int, Float64}()
    V̄_prev = Dict{Int, Float64}()

    ρθ = config.rho_theta
    ρV = config.rho_v

    res_primal_hist = Float64[]
    res_dual_hist = Float64[]
    eps_primal_hist = Float64[]
    eps_dual_hist = Float64[]
    rhoθ_hist = Float64[]
    rhoV_hist = Float64[]
    θ̄_hist = Dict{Int, Vector{Float64}}()

    μ_incr, μ_decr = 2.0, 2.0
    τ_incr, τ_decr = 1.05, 0.98

    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual

    @info "[ADMM-SOCP-bus] ρθ=$(ρθ) ρV=$(ρV) eps_pri=$(ϵ_abs) eps_du=$(ϵ_rel)"

    for k in 1:config.max_iter
        println("\n[ADMM-SOCP-bus] Iteração $k")

        for (zone_id, subdata) in partitions
            local_model = build_model(subdata, SOCP; optimizer = Gurobi.Optimizer)

            @expression(local_model, admm_penalty, 0.0)
            bus_ids = subdata["bus"]["bus_i"]

            for b in shared_buses
                if b in bus_ids
                    θloc = local_model[:theta][b, 1]
                    Vloc = local_model[:v_magnitude][b, 1]

                    θbar = get(θ̄, b, 0.0)
                    Vbar = get(V̄, b, 1.0)

                    λθ = get(duals_theta[zone_id], b, 0.0)
                    λV = get(duals_vmag[zone_id], b, 0.0)

                    admm_penalty += λθ * (θloc - θbar) + (ρθ / 2) * (θloc - θbar)^2
                    admm_penalty += λV * (Vloc - Vbar) + (ρV / 2) * (Vloc - Vbar)^2
                end
            end

            @objective(local_model, Min, local_model[:generation_cost] + admm_penalty)
            optimize!(local_model)

            # print(local_model)

            subdata["model"] = local_model
            subdata["generation_cost_model_expr"] = value(local_model[:generation_cost])
            subdata["real_generation_cost"] =
                compute_real_generation_cost(local_model, subdata["objective_params"])
        end

        θ̄_prev = isempty(θ̄) ? θ̄ : deepcopy(θ̄)
        V̄_prev = isempty(V̄) ? V̄ : deepcopy(V̄)

        for b in shared_buses
            θs = Float64[]
            Vacc = 0.0
            n = 0
            for (_, subdata) in partitions
                if b in subdata["bus"]["bus_i"]
                    m = subdata["model"]
                    push!(θs, value(m[:theta][b, 1]))
                    Vacc += value(m[:v_magnitude][b, 1])
                    n += 1
                end
            end
            if n > 0
                θ̄[b] = _circular_mean(θs)
                V̄[b] = Vacc / n
                push!(get!(θ̄_hist, b, Float64[]), θ̄[b])
            end
        end

        for (zone_id, subdata) in partitions
            m = subdata["model"]
            bus_ids = subdata["bus"]["bus_i"]
            for b in shared_buses
                if b in bus_ids
                    θloc = value(m[:theta][b, 1])
                    Vloc = value(m[:v_magnitude][b, 1])
                    dθ = _wrap_angle_diff(θloc - θ̄[b])
                    dV = Vloc - V̄[b]
                    duals_theta[zone_id][b] = get(duals_theta[zone_id], b, 0.0) + ρθ * dθ
                    duals_vmag[zone_id][b] = get(duals_vmag[zone_id], b, 0.0) + ρV * dV
                end
            end
        end

        r_pri_sq = 0.0
        for (_, subdata) in partitions
            m = subdata["model"]
            bus_ids = subdata["bus"]["bus_i"]
            for b in shared_buses
                if b in bus_ids
                    dθ = _wrap_angle_diff(value(m[:theta][b, 1]) - θ̄[b])
                    dV = value(m[:v_magnitude][b, 1]) - V̄[b]
                    r_pri_sq += dθ^2 + dV^2
                end
            end
        end
        r_pri = sqrt(r_pri_sq)

        s_du = sqrt(
            sum(_wrap_angle_diff(θ̄[b] - get(θ̄_prev, b, 0.0))^2 for b in shared_buses) +
            sum((V̄[b] - get(V̄_prev, b, 0.0))^2 for b in shared_buses),
        )

        push!(res_primal_hist, r_pri)
        push!(res_dual_hist, s_du)

        dim_pri = 0
        for (_, sub) in partitions
            dim_pri += 2 * sum(b in sub["bus"]["bus_i"] for b in shared_buses)
        end
        dim_du = 2 * length(shared_buses)

        ε_pri =
            sqrt(max(1, dim_pri)) * ϵ_abs +
            ϵ_rel * max(
                sqrt(
                    sum(θ̄[b]^2 for b in shared_buses) + sum(V̄[b]^2 for b in shared_buses),
                ),
                1e-12,
            )
        ε_du =
            sqrt(max(1, dim_du)) * ϵ_abs +
            ϵ_rel * sqrt(
                sum(
                    _wrap_angle_diff(θ̄[b] - get(θ̄_prev, b, 0.0))^2 for b in shared_buses
                ) +
                sum((V̄[b] - get(V̄_prev, b, 0.0))^2 for b in shared_buses),
            )

        push!(eps_primal_hist, ε_pri)
        push!(eps_dual_hist, ε_du)

        println("‖r‖ = $r_pri | ε_pri = $ε_pri")
        println("‖s‖ = $s_du  | ε_du  = $ε_du")

        push!(rhoθ_hist, ρθ)
        push!(rhoV_hist, ρV)

        if r_pri < ε_pri && s_du < ε_du
            println("Stopping criterion met.")
            break
        end

        if r_pri > μ_incr * s_du
            ρθ *= τ_incr
            ρV *= τ_incr
        elseif s_du > μ_decr * r_pri
            ρθ *= τ_decr
            ρV *= τ_decr
        end
    end

    plot_admm_convergence(res_primal_hist, res_dual_hist, eps_primal_hist, eps_dual_hist,
        "plot_convergencia_admm_socp.html")

    PlotlyLight.save(
        (
            () -> begin
                p = Plot()
                p(;
                    x = 1:length(rhoθ_hist),
                    y = rhoθ_hist,
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
        )(), "plot_rho_theta_socp.html")

    PlotlyLight.save(
        (
            () -> begin
                p = Plot()
                p(;
                    x = 1:length(rhoV_hist),
                    y = rhoV_hist,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρV",
                )
                p.layout.title = "Evolução de ρV"
                p.layout.xaxis.title = "Iterações"
                p.layout.yaxis.title = "ρV"
                p.layout.yaxis.type = "log"
                p
            end
        )(), "plot_rho_v_socp.html")

    plot_theta_avg_heatmap(θ̄_hist, "plot_theta_avg_heatmap_socp.html")

    total_real_cost = sum(p["real_generation_cost"] for p in values(partitions))
    total_expr_cost = sum(p["generation_cost_model_expr"] for p in values(partitions))

    println("\n[ADMM-SOCP-bus] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("ρθ_final = $ρθ | ρV_final = $ρV")

    return total_real_cost, partitions
end
