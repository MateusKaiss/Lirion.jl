@inline function _rt_from_exprs_socp(model::Model, l_local::Int; t::Int = 1)
    return model[:R][l_local, t], model[:T_var][l_local, t]
end

function run_admm_socp(config::DistributedOPFModel, data::Dict{String,Any})
    @assert config.base_model == SOCP "ADMM-SOCP: base_model deve ser SOCP."

    partitions, shared_buses = prepare_admm_partitions(data)
    cut_edges = _find_cut_edges(data)

    # ------------------- duais -------------------
    duals_v2 = Dict{Int, Dict{Int, Float64}}()   # λV2[zone][bus]
    duals_R  = Dict{Int, Dict{Int, Float64}}()   # λR[zone][edge]
    duals_T  = Dict{Int, Dict{Int, Float64}}()   # λT[zone][edge]
    for z in keys(partitions)
        duals_v2[z] = Dict{Int, Float64}()
        duals_R[z]  = Dict{Int, Float64}()
        duals_T[z]  = Dict{Int, Float64}()
    end

    # ---------------- consenso -------------------
    v2_avg     = Dict{Int, Float64}()  # V^2 médio por barra compartilhada
    v2_prev    = Dict{Int, Float64}()
    zR         = Dict{Int, Float64}()  # consenso R por aresta de corte
    zT         = Dict{Int, Float64}()  # consenso T por aresta de corte
    zR_prev    = Dict{Int, Float64}()
    zT_prev    = Dict{Int, Float64}()
    for e in cut_edges
        zR[e] = 0.0
        zT[e] = 0.0
    end

    # ---------------- penalidades ----------------
    ρV2 = config.rho_v       
    ρR  = config.rho_flow_p  
    ρT  = config.rho_flow_q 
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual

    res_primal_vec = Float64[]
    res_dual_vec   = Float64[]
    eps_primal_vec = Float64[]
    eps_dual_vec   = Float64[]
    rhoV2_vec      = Float64[]
    rhoR_vec       = Float64[]
    rhoT_vec       = Float64[]
    μ_incr, μ_decr = 2.0, 2.0
    τ_incr, τ_decr = 1.05, 0.98

    for k in 1:config.max_iter
        println("\n[ADMM-SOCP(no-θ)] Iteração $k")

        # 1) resolve subproblemas locais (SOCP)
        for (zone_id, subdata) in partitions
            local_model = build_model(subdata, SOCP; optimizer = SCS.Optimizer)
            set_optimizer_attribute(local_model, "eps_abs", 1e-6)
            set_optimizer_attribute(local_model, "eps_rel", 1e-6)
            set_optimizer_attribute(local_model, "max_iters", 200_000)

            @expression(local_model, admm_penalty, 0.0)

            # --- penalização por BUS (V^2) ---
            bus_ids = subdata["bus"]["bus_i"]
            for b in shared_buses
                if b in bus_ids
                    V2  = local_model[:v_magnitude_sq][b, 1]
                    V2̄ = get(v2_avg, b, 1.0)
                    λV2 = get(duals_v2[zone_id], b, 0.0)
                    admm_penalty += λV2 * (V2 - V2̄) + (ρV2 / 2) * (V2 - V2̄)^2
                end
            end

            # --- penalização por BRANCH (R,T) nas arestas de corte ---
            zone_edges = _zone_cut_edges(subdata, cut_edges)
            for e in zone_edges
                l_local = findfirst(==(e), subdata["branch_local_indices"])
                @assert l_local !== nothing "Aresta global $e não está na partição $zone_id"
                Rloc, Tloc = _rt_from_exprs_socp(local_model, l_local)
                λR = get(duals_R[zone_id], e, 0.0)
                λT = get(duals_T[zone_id], e, 0.0)
                admm_penalty += λR * (Rloc - zR[e]) + (ρR / 2) * (Rloc - zR[e])^2
                admm_penalty += λT * (Tloc - zT[e]) + (ρT / 2) * (Tloc - zT[e])^2
            end

            @objective(local_model, Min, local_model[:generation_cost] + admm_penalty)
            optimize!(local_model)

            subdata["model"] = local_model
            subdata["generation_cost_model_expr"] = value(local_model[:generation_cost])
            subdata["real_generation_cost"] =
                compute_real_generation_cost(local_model, subdata["objective_params"])
        end

        # 2) atualiza consensos
        v2_prev = isempty(v2_avg) ? v2_avg : deepcopy(v2_avg)
        zR_prev = deepcopy(zR)
        zT_prev = deepcopy(zT)

        # V^2 (bus)
        for b in shared_buses
            acc = 0.0; n = 0
            for (_, subdata) in partitions
                if b in subdata["bus"]["bus_i"]
                    acc += value(subdata["model"][:v_magnitude_sq][b, 1])
                    n += 1
                end
            end
            if n > 0
                v2_avg[b] = acc / n
            end
        end

        # R,T (edges)
        for e in cut_edges
            accR = 0.0; accT = 0.0; n = 0
            for (_, subdata) in partitions
                l_local = findfirst(==(e), subdata["branch_local_indices"])
                if l_local !== nothing
                    Rloc, Tloc = _rt_from_exprs_socp(subdata["model"], l_local)
                    accR += value(Rloc); accT += value(Tloc); n += 1
                end
            end
            if n > 0
                zR[e] = accR / n
                zT[e] = accT / n
            end
        end

        # 3) atualiza duais
        for (zone_id, subdata) in partitions
            m = subdata["model"]
            # BUS
            for b in shared_buses
                if b in subdata["bus"]["bus_i"]
                    V2loc = value(m[:v_magnitude_sq][b, 1])
                    duals_v2[zone_id][b] =
                        get(duals_v2[zone_id], b, 0.0) + ρV2 * (V2loc - v2_avg[b])
                end
            end
            # BRANCH
            for e in _zone_cut_edges(subdata, cut_edges)
                l_local = findfirst(==(e), subdata["branch_local_indices"])
                Rloc, Tloc = _rt_from_exprs_socp(m, l_local)
                duals_R[zone_id][e] =
                    get(duals_R[zone_id], e, 0.0) + ρR * (value(Rloc) - zR[e])
                duals_T[zone_id][e] =
                    get(duals_T[zone_id], e, 0.0) + ρT * (value(Tloc) - zT[e])
            end
        end

        # 4) resíduos
        r_pri_sq = 0.0
        for (_, subdata) in partitions
            m = subdata["model"]
            # V^2
            for b in shared_buses
                if b in subdata["bus"]["bus_i"]
                    r_pri_sq += (value(m[:v_magnitude_sq][b, 1]) - v2_avg[b])^2
                end
            end
            # R,T
            for e in _zone_cut_edges(subdata, cut_edges)
                l_local = findfirst(==(e), subdata["branch_local_indices"])
                Rloc, Tloc = _rt_from_exprs_socp(m, l_local)
                r_pri_sq += (value(Rloc) - zR[e])^2 + (value(Tloc) - zT[e])^2
            end
        end
        r_pri = sqrt(r_pri_sq)

        s_du = sqrt(
            sum((v2_avg[b] - get(v2_prev, b, 0.0))^2 for b in shared_buses) +
            sum((zR[e] - get(zR_prev, e, 0.0))^2 for e in cut_edges) +
            sum((zT[e] - get(zT_prev, e, 0.0))^2 for e in cut_edges)
        )

        push!(res_primal_vec, r_pri)
        push!(res_dual_vec,   s_du)

        # dimensões para as tolerâncias
        dim_pri = 0
        for (_, sub) in partitions
            dim_pri += sum(b in sub["bus"]["bus_i"] for b in shared_buses)   # V^2
            dim_pri += 2 * length(_zone_cut_edges(sub, cut_edges))           # R e T
        end
        dim_du = length(shared_buses) + 2 * length(cut_edges)

        ε_pri =
            sqrt(max(1, dim_pri)) * ϵ_abs +
            ϵ_rel * max(
                sqrt(
                    sum(v2_avg[b]^2 for b in shared_buses) +
                    sum(zR[e]^2 for e in cut_edges) + sum(zT[e]^2 for e in cut_edges)
                ),
                1e-12,
            )
        ε_du =
            sqrt(max(1, dim_du)) * ϵ_abs +
            ϵ_rel * sqrt(
                sum((v2_avg[b] - get(v2_prev, b, 0.0))^2 for b in shared_buses) +
                sum((zR[e] - get(zR_prev, e, 0.0))^2 for e in cut_edges) +
                sum((zT[e] - get(zT_prev, e, 0.0))^2 for e in cut_edges)
            )

        push!(eps_primal_vec, ε_pri)
        push!(eps_dual_vec,   ε_du)

        println("‖r‖ = $r_pri | ε_pri = $ε_pri")
        println("‖s‖ = $s_du  | ε_du  = $ε_du")

        push!(rhoV2_vec, ρV2)
        push!(rhoR_vec,  ρR)
        push!(rhoT_vec,  ρT)

        if r_pri < ε_pri && s_du < ε_du
            println("Stopping criterion met.")
            break
        end

        # 5) ajuste adaptativo (residual balancing)
        if r_pri > μ_incr * s_du
            ρV2 *= τ_incr;  ρR *= τ_incr;  ρT *= τ_incr
        elseif s_du > μ_decr * r_pri
            ρV2 *= τ_decr;  ρR *= τ_decr;  ρT *= τ_decr
        end
    end

    # 6) plots e resumo
    plot_admm_convergence(res_primal_vec, res_dual_vec, eps_primal_vec, eps_dual_vec,
        "plot_convergencia_admm_socp_no_theta.html")

    PlotlyLight.save((()->begin
        p = Plot(); p(; x=1:length(rhoV2_vec), y=rhoV2_vec, type=:scatter, mode="lines+markers", name="ρ(V²)")
        p.layout.title="Evolução de ρ(V²)"; p.layout.xaxis.title="Iterações"; p.layout.yaxis.title="ρ(V²)"; p.layout.yaxis.type="log"; p
    end)(), "plot_rho_v2_socp_no_theta.html")

    PlotlyLight.save((()->begin
        p = Plot(); p(; x=1:length(rhoR_vec), y=rhoR_vec, type=:scatter, mode="lines+markers", name="ρR")
        p.layout.title="Evolução de ρR"; p.layout.xaxis.title="Iterações"; p.layout.yaxis.title="ρR"; p.layout.yaxis.type="log"; p
    end)(), "plot_rho_R_socp_no_theta.html")

    PlotlyLight.save((()->begin
        p = Plot(); p(; x=1:length(rhoT_vec), y=rhoT_vec, type=:scatter, mode="lines+markers", name="ρT")
        p.layout.title="Evolução de ρT"; p.layout.xaxis.title="Iterações"; p.layout.yaxis.title="ρT"; p.layout.yaxis.type="log"; p
    end)(), "plot_rho_T_socp_no_theta.html")

    total_real_cost = sum(p["real_generation_cost"] for p in values(partitions))
    total_expr_cost = sum(p["generation_cost_model_expr"] for p in values(partitions))
    println("\n[ADMM-SOCP(no-θ)] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("ρ(V²)_final = $ρV2 | ρR_final = $ρR | ρT_final = $ρT")

    return total_real_cost, partitions
end
