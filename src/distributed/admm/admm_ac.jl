# diferença angular "wrap" em (-π, π]
@inline function _wrap_angle_diff(Δ::Float64)
    return atan(sin(Δ), cos(Δ))
end

# média circular de ângulos θ₁, θ₂, ...
# retorna atan(Σsin θ, Σcos θ)
@inline function _circular_mean(θs::AbstractVector{<:Real})
    sum_sin = 0.0
    sum_cos = 0.0
    @inbounds @simd for θ in θs
        sum_sin += sin(θ)
        sum_cos += cos(θ)
    end
    return atan(sum_sin, sum_cos)
end

# ------------------------------------------------------------------------------
# ADMM-AC por barra com consenso em θ e |V|
#   - consenso: θ (média circular) e V (média aritmética)
#   - resíduos e atualização de duais: usam wrap para o θ
# ------------------------------------------------------------------------------
function run_admm_ac(config::DistributedOPFModel, data::Dict{String, Any})
    @assert config.base_model == AC "ADMM-AC: base_model deve ser AC."
    partitions, shared_buses = prepare_admm_partitions(data)
    duals_theta = Dict{Int, Dict{Int, Float64}}()
    duals_vmag = Dict{Int, Dict{Int, Float64}}()
    θ_avg = Dict{Int, Float64}()
    v_avg = Dict{Int, Float64}()
    θ_avg_prev = Dict{Int, Float64}()
    v_avg_prev = Dict{Int, Float64}()

    for z in keys(partitions)
        duals_theta[z] = Dict{Int, Float64}()
        duals_vmag[z] = Dict{Int, Float64}()
    end

    ρθ = config.rho_theta
    ρV = config.rho_v
    @info "[ADMM-AC] rho_theta=$(ρθ), rho_v=$(ρV), eps_pri=$(config.ε_primal), eps_dual=$(config.ε_dual)"

    res_primal_vec = Float64[]
    res_dual_vec = Float64[]
    eps_primal_vec = Float64[]
    eps_dual_vec = Float64[]
    rhoθ_vec = Float64[]
    rhoV_vec = Float64[]
    θ_avg_history = Dict{Int, Vector{Float64}}()

    μ_incr, μ_decr = 2.0, 2.0
    τ_incr, τ_decr = 1.05, 0.98

    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual

    for k in 1:config.max_iter
        println("\n[ADMM-AC] Iteração $k")

        for (zone_id, subdata) in partitions
            local_model = build_model(subdata, AC; optimizer = Ipopt.Optimizer)

            @expression(local_model, admm_penalty, 0.0)
            bus_ids = subdata["bus"]["bus_i"]

            for b in shared_buses
                if b in bus_ids
                    # Ângulo
                    θ = local_model[:theta][b, 1]
                    θ̄ = get(θ_avg, b, 0.0)
                    λθ = get(duals_theta[zone_id], b, 0.0)
                    admm_penalty += λθ * (θ - θ̄) + (ρθ / 2) * (θ - θ̄)^2

                    # Magnitude de tensão
                    V = local_model[:v_magnitude][b, 1]
                    V̄ = get(v_avg, b, 1.0)
                    λV = get(duals_vmag[zone_id], b, 0.0)
                    admm_penalty += λV * (V - V̄) + (ρV / 2) * (V - V̄)^2
                end
            end

            @objective(local_model, Min, local_model[:generation_cost] + admm_penalty)
            optimize!(local_model)

            subdata["model"] = local_model
            subdata["generation_cost_model_expr"] = value(local_model[:generation_cost])
            subdata["real_generation_cost"] =
                compute_real_generation_cost(local_model, subdata["objective_params"])
        end

        θ_avg_prev = isempty(θ_avg) ? θ_avg : deepcopy(θ_avg)
        v_avg_prev = isempty(v_avg) ? v_avg : deepcopy(v_avg)

        for b in shared_buses
            θ_list = Float64[]
            V_acc = 0.0
            n = 0
            for (_, subdata) in partitions
                if b in subdata["bus"]["bus_i"]
                    m = subdata["model"]
                    push!(θ_list, value(m[:theta][b, 1]))
                    V_acc += value(m[:v_magnitude][b, 1])
                    n += 1
                end
            end
            if n > 0
                θ_avg[b] = _circular_mean(θ_list)
                v_avg[b] = V_acc / n
                push!(get!(θ_avg_history, b, Float64[]), θ_avg[b])
            end
        end

        for (zone_id, subdata) in partitions
            m = get(subdata, "model", nothing)
            if m === nothing
                continue
            end
            bus_ids = subdata["bus"]["bus_i"]
            for b in shared_buses
                if b in bus_ids
                    θloc = value(m[:theta][b, 1])
                    Vloc = value(m[:v_magnitude][b, 1])
                    dθ = _wrap_angle_diff(θloc - θ_avg[b])
                    dV = Vloc - v_avg[b]
                    duals_theta[zone_id][b] = get(duals_theta[zone_id], b, 0.0) + ρθ * dθ
                    duals_vmag[zone_id][b] = get(duals_vmag[zone_id], b, 0.0) + ρV * dV
                end
            end
        end

        res_primal_sq = 0.0
        for (_, subdata) in partitions
            m = get(subdata, "model", nothing)
            if m === nothing
                continue
            end
            bus_ids = subdata["bus"]["bus_i"]
            for b in shared_buses
                if b in bus_ids
                    dθ = _wrap_angle_diff(value(m[:theta][b, 1]) - θ_avg[b])
                    dV = value(m[:v_magnitude][b, 1]) - v_avg[b]
                    res_primal_sq += dθ^2 + dV^2
                end
            end
        end
        res_primal = sqrt(res_primal_sq)

        res_dual = sqrt(
            sum(
                _wrap_angle_diff(θ_avg[b] - get(θ_avg_prev, b, 0.0))^2 for
                b in shared_buses
            ) +
            sum((v_avg[b] - get(v_avg_prev, b, 0.0))^2 for b in shared_buses),
        )

        push!(res_primal_vec, res_primal)
        push!(res_dual_vec, res_dual)

        dim_pri =
            2 * sum(
                b in partitions[z]["bus"]["bus_i"] for z in keys(partitions),
                b in shared_buses
            )
        dim_du = 2 * length(shared_buses)

        ϵ_primal =
            sqrt(max(1, dim_pri)) * ϵ_abs +
            ϵ_rel * max(
                sqrt(
                    sum(θ_avg[b]^2 for b in shared_buses) +
                    sum(v_avg[b]^2 for b in shared_buses),
                ),
                1e-12,
            )
        ϵ_dual =
            sqrt(max(1, dim_du)) * ϵ_abs +
            ϵ_rel * sqrt(
                sum(
                    _wrap_angle_diff(θ_avg[b] - get(θ_avg_prev, b, 0.0))^2 for
                    b in shared_buses
                ) +
                sum((v_avg[b] - get(v_avg_prev, b, 0.0))^2 for b in shared_buses),
            )

        push!(eps_primal_vec, ϵ_primal)
        push!(eps_dual_vec, ϵ_dual)

        println("‖r‖ = $res_primal | ε_pri = $ϵ_primal")
        println("‖s‖ = $res_dual   | ε_dual = $ϵ_dual")

        push!(rhoθ_vec, ρθ)
        push!(rhoV_vec, ρV)

        if res_primal < ϵ_primal && res_dual < ϵ_dual
            println("Stopping criterion met.")
            break
        end

        if res_primal > μ_incr * res_dual
            ρθ *= τ_incr
            ρV *= τ_incr
        elseif res_dual > μ_decr * res_primal
            ρθ *= τ_decr
            ρV *= τ_decr
        end
    end

    plot_admm_convergence(res_primal_vec, res_dual_vec, eps_primal_vec, eps_dual_vec,
        "plot_convergencia_admm_ac.html")

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
        )(), "plot_rho_theta_ac.html")

    PlotlyLight.save(
        (
            () -> begin
                p = Plot()
                p(;
                    x = 1:length(rhoV_vec),
                    y = rhoV_vec,
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
        )(), "plot_rho_vmag_ac.html")

    plot_theta_avg_heatmap(θ_avg_history, "plot_theta_avg_heatmap_ac.html")

    total_real_cost = sum(p["real_generation_cost"] for p in values(partitions))
    total_expr_cost = sum(p["generation_cost_model_expr"] for p in values(partitions))
    println("\n[ADMM-AC] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("ρθ_final = $ρθ | ρV_final = $ρV")

    return total_real_cost, partitions
end
