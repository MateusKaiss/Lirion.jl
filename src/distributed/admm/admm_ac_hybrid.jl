# --------------------------------------------------------------------------
# ADMM-AC híbrido:
#   - Consenso por BARRA: θ (média circular) e |V| (média aritmética)
#   - Consenso por RAMO (lado "from"): P_flow e Q_flow (médias aritméticas)
#   - Atualização adaptativa de ρ em 2 blocos:
#        • (ρθ, ρV) pelo balanço de resíduos BUS
#        • (ρP, ρQ) pelo balanço de resíduos BRANCH
# --------------------------------------------------------------------------
function run_admm_ac_hybrid(config::DistributedOPFModel, data::Dict{String, Any})
    @assert config.base_model == AC "Hybrid AC: base_model deve ser AC."

    partitions, shared_buses = prepare_admm_partitions(data)
    cut_edges = _find_cut_edges(data)

    dualsθ = Dict{Int, Dict{Int, Float64}}()
    dualsV = Dict{Int, Dict{Int, Float64}}()
    dualsP = Dict{Int, Dict{Int, Float64}}()
    dualsQ = Dict{Int, Dict{Int, Float64}}()
    for z in keys(partitions)
        dualsθ[z] = Dict{Int, Float64}()
        dualsV[z] = Dict{Int, Float64}()
        dualsP[z] = Dict{Int, Float64}()
        dualsQ[z] = Dict{Int, Float64}()
    end

    θ̄ = Dict{Int, Float64}()
    V̄ = Dict{Int, Float64}()
    zP = Dict{Int, Float64}()
    zQ = Dict{Int, Float64}()

    θ̄_prev = Dict{Int, Float64}()
    V̄_prev = Dict{Int, Float64}()
    zP_prev = Dict{Int, Float64}()
    zQ_prev = Dict{Int, Float64}()

    for e in cut_edges
        zP[e] = 0.0
        zQ[e] = 0.0
    end

    ρθ = config.rho_theta
    ρV = config.rho_v
    ρP = config.rho_flow_p
    ρQ = config.rho_flow_q
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual

    res_pri_hist = Float64[]
    res_du_hist = Float64[]
    eps_pri_hist = Float64[]
    eps_du_hist = Float64[]
    rhoθ_hist = Float64[]
    rhoV_hist = Float64[]
    rhoP_hist = Float64[]
    rhoQ_hist = Float64[]
    θ̄_hist = Dict{Int, Vector{Float64}}()

    μ_incr, μ_decr = 2.0, 2.0
    τ_incr, τ_decr = 1.05, 0.98

    @info "[ADMM-AC-hybrid] ρθ=$(ρθ), ρV=$(ρV), ρP=$(ρP), ρQ=$(ρQ), eps_pri=$(ϵ_abs), eps_dual=$(ϵ_rel)"

    for k in 1:config.max_iter
        println("\n[ADMM-AC-hybrid] Iteração $k")

        for (zone_id, subdata) in partitions
            local_model = build_model(subdata, AC; optimizer = Ipopt.Optimizer)

            @expression(local_model, admm_penalty, 0.0)

            # --- BUS (θ, V) ---
            bus_ids = subdata["bus"]["bus_i"]
            for b in shared_buses
                if b in bus_ids
                    θloc = local_model[:theta][b, 1]
                    Vloc = local_model[:v_magnitude][b, 1]
                    θbar = get(θ̄, b, 0.0)
                    Vbar = get(V̄, b, 1.0)
                    λθ = get(dualsθ[zone_id], b, 0.0)
                    λV = get(dualsV[zone_id], b, 0.0)
                    admm_penalty += λθ * (θloc - θbar) + (ρθ / 2) * (θloc - θbar)^2
                    admm_penalty += λV * (Vloc - Vbar) + (ρV / 2) * (Vloc - Vbar)^2
                end
            end

            # --- BRANCH (P,Q) ---
            zone_edges = _zone_cut_edges(subdata, cut_edges)
            for e in zone_edges
                l_local = findfirst(==(e), subdata["branch_local_indices"])
                @assert l_local !== nothing "aresta global $e não está na partição $zone_id"
                Ploc, Qloc = _ac_flow_from_exprs_existing(local_model, subdata, l_local)
                λP = get(dualsP[zone_id], e, 0.0)
                λQ = get(dualsQ[zone_id], e, 0.0)
                admm_penalty += λP * (Ploc - zP[e]) + (ρP / 2) * (Ploc - zP[e])^2
                admm_penalty += λQ * (Qloc - zQ[e]) + (ρQ / 2) * (Qloc - zQ[e])^2
            end

            @objective(local_model, Min, local_model[:generation_cost] + admm_penalty)
            optimize!(local_model)

            subdata["model"] = local_model
            subdata["generation_cost_model_expr"] = value(local_model[:generation_cost])
            subdata["real_generation_cost"] =
                compute_real_generation_cost(local_model, subdata["objective_params"])
        end

        θ̄_prev = isempty(θ̄) ? θ̄ : deepcopy(θ̄)
        V̄_prev = isempty(V̄) ? V̄ : deepcopy(V̄)
        zP_prev = deepcopy(zP)
        zQ_prev = deepcopy(zQ)

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

        for e in cut_edges
            accP = 0.0
            accQ = 0.0
            num = 0
            for (_, subdata) in partitions
                l_local = findfirst(==(e), subdata["branch_local_indices"])
                if l_local !== nothing
                    Ploc, Qloc =
                        _ac_flow_from_exprs_existing(subdata["model"], subdata, l_local)
                    accP += value(Ploc)
                    accQ += value(Qloc)
                    num += 1
                end
            end
            if num > 0
                zP[e] = accP / num
                zQ[e] = accQ / num
            end
        end

        for (zone_id, subdata) in partitions
            m = subdata["model"]
            bus_ids = subdata["bus"]["bus_i"]

            # BUS duais
            for b in shared_buses
                if b in bus_ids
                    θloc = value(m[:theta][b, 1])
                    Vloc = value(m[:v_magnitude][b, 1])
                    dθ = _wrap_angle_diff(θloc - θ̄[b])
                    dV = Vloc - V̄[b]
                    dualsθ[zone_id][b] = get(dualsθ[zone_id], b, 0.0) + ρθ * dθ
                    dualsV[zone_id][b] = get(dualsV[zone_id], b, 0.0) + ρV * dV
                end
            end

            # BRANCH duais
            for e in _zone_cut_edges(subdata, cut_edges)
                l_local = findfirst(==(e), subdata["branch_local_indices"])
                Ploc, Qloc = _ac_flow_from_exprs_existing(m, subdata, l_local)
                dualsP[zone_id][e] =
                    get(dualsP[zone_id], e, 0.0) + ρP * (value(Ploc) - zP[e])
                dualsQ[zone_id][e] =
                    get(dualsQ[zone_id], e, 0.0) + ρQ * (value(Qloc) - zQ[e])
            end
        end

        r_bus_sq = 0.0
        r_branch_sq = 0.0
        for (_, subdata) in partitions
            m = subdata["model"]
            bus_ids = subdata["bus"]["bus_i"]
            # BUS
            for b in shared_buses
                if b in bus_ids
                    dθ = _wrap_angle_diff(value(m[:theta][b, 1]) - θ̄[b])
                    dV = value(m[:v_magnitude][b, 1]) - V̄[b]
                    r_bus_sq += dθ^2 + dV^2
                end
            end
            # BRANCH
            for e in _zone_cut_edges(subdata, cut_edges)
                l_local = findfirst(==(e), subdata["branch_local_indices"])
                Ploc, Qloc = _ac_flow_from_exprs_existing(m, subdata, l_local)
                r_branch_sq += (value(Ploc) - zP[e])^2 + (value(Qloc) - zQ[e])^2
            end
        end
        r_pri = sqrt(r_bus_sq + r_branch_sq)

        s_bus = sqrt(
            sum(_wrap_angle_diff(θ̄[b] - get(θ̄_prev, b, 0.0))^2 for b in shared_buses) +
            sum((V̄[b] - get(V̄_prev, b, 0.0))^2 for b in shared_buses),
        )
        s_branch = sqrt(
            sum((zP[e] - get(zP_prev, e, 0.0))^2 for e in cut_edges) +
            sum((zQ[e] - get(zQ_prev, e, 0.0))^2 for e in cut_edges),
        )
        s_tot = sqrt(s_bus^2 + s_branch^2)

        push!(res_pri_hist, r_pri)
        push!(res_du_hist, s_tot)

        dim_pri = 0
        for (_, sub) in partitions
            dim_pri += 2 * sum(b in sub["bus"]["bus_i"] for b in shared_buses)     # θ,V
            dim_pri += 2 * length(_zone_cut_edges(sub, cut_edges))                 # P,Q
        end
        dim_du = 2 * length(shared_buses) + 2 * length(cut_edges)                  # θ̄,V̄,P̄,Q̄

        ε_pri =
            sqrt(max(1, dim_pri)) * ϵ_abs +
            ϵ_rel * max(
                sqrt(
                    sum(θ̄[b]^2 for b in shared_buses) +
                    sum(V̄[b]^2 for b in shared_buses) +
                    sum(zP[e]^2 for e in cut_edges) + sum(zQ[e]^2 for e in cut_edges),
                ),
                1e-12,
            )
        ε_du =
            sqrt(max(1, dim_du)) * ϵ_abs +
            ϵ_rel * sqrt(
                sum(
                    _wrap_angle_diff(θ̄[b] - get(θ̄_prev, b, 0.0))^2 for b in shared_buses
                ) +
                sum((V̄[b] - get(V̄_prev, b, 0.0))^2 for b in shared_buses) +
                sum((zP[e] - get(zP_prev, e, 0.0))^2 for e in cut_edges) +
                sum((zQ[e] - get(zQ_prev, e, 0.0))^2 for e in cut_edges),
            )

        push!(eps_pri_hist, ε_pri)
        push!(eps_du_hist, ε_du)

        println("‖r‖ = $r_pri | ε_pri = $ε_pri")
        println("‖s‖ = $s_tot | ε_du  = $ε_du")

        push!(rhoθ_hist, ρθ)
        push!(rhoV_hist, ρV)
        push!(rhoP_hist, ρP)
        push!(rhoQ_hist, ρQ)

        # Parada
        if r_pri < ε_pri && s_tot < ε_du
            println("Stopping criterion met.")
            break
        end

        # BUS: (θ, V)
        if r_bus_sq > (μ_incr^2) * (s_bus^2) && s_bus > 0.0
            ρθ *= τ_incr
            ρV *= τ_incr
        elseif s_bus^2 > (μ_decr^2) * r_bus_sq && r_bus_sq > 0.0
            ρθ *= τ_decr
            ρV *= τ_decr
        end
        # BRANCH: (P, Q)
        if r_branch_sq > (μ_incr^2) * (s_branch^2) && s_branch > 0.0
            ρP *= τ_incr
            ρQ *= τ_incr
        elseif s_branch^2 > (μ_decr^2) * r_branch_sq && r_branch_sq > 0.0
            ρP *= τ_decr
            ρQ *= τ_decr
        end
    end

    plot_admm_convergence(res_pri_hist, res_du_hist, eps_pri_hist, eps_du_hist,
        "plot_convergencia_admm_ac_hybrid.html")

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
        )(), "plot_rho_theta_ac_hybrid.html")

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
        )(), "plot_rho_v_ac_hybrid.html")

    PlotlyLight.save(
        (
            () -> begin
                p = Plot()
                p(;
                    x = 1:length(rhoP_hist),
                    y = rhoP_hist,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρP",
                )
                p.layout.title = "Evolução de ρP"
                p.layout.xaxis.title = "Iterações"
                p.layout.yaxis.title = "ρP"
                p.layout.yaxis.type = "log"
                p
            end
        )(), "plot_rho_p_ac_hybrid.html")

    PlotlyLight.save(
        (
            () -> begin
                p = Plot()
                p(;
                    x = 1:length(rhoQ_hist),
                    y = rhoQ_hist,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρQ",
                )
                p.layout.title = "Evolução de ρQ"
                p.layout.xaxis.title = "Iterações"
                p.layout.yaxis.title = "ρQ"
                p.layout.yaxis.type = "log"
                p
            end
        )(), "plot_rho_q_ac_hybrid.html")

    plot_theta_avg_heatmap(θ̄_hist, "plot_theta_avg_heatmap_ac_hybrid.html")

    total_real_cost = sum(p["real_generation_cost"] for p in values(partitions))
    total_expr_cost = sum(p["generation_cost_model_expr"] for p in values(partitions))

    println("\n[ADMM-AC-hybrid] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("ρθ_final = $ρθ | ρV_final = $ρV | ρP_final = $ρP | ρQ_final = $ρQ")

    return total_real_cost, partitions
end
