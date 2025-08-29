# ================= AC Utils =================
# ---------- Solver ----------
function resolve_optimizer(config)
    if hasproperty(config, :optimizer) && getproperty(config, :optimizer) !== nothing
        return getproperty(config, :optimizer)
    else
        return Ipopt.Optimizer
    end
end

function add_penalty_bus_ac!(m, subdata, shared_buses::Vector{Int},
                             θ_avg::Dict{Int,Float64}, v_avg::Dict{Int,Float64},
                             dualsθ_zone::Dict{Int,Float64}, dualsV_zone::Dict{Int,Float64},
                             ρθ::Float64, ρV::Float64)
    bus_ids = subdata["bus"]["bus_i"]
    @expression(m, admm_penalty, 0.0)
    for b in shared_buses
        if b in bus_ids
            θ  = m[:theta][b,1]; θ̄ = get(θ_avg, b, 0.0); λθ = get(dualsθ_zone, b, 0.0)
            V  = m[:v_magnitude][b,1]; V̄ = get(v_avg, b, 1.0); λV = get(dualsV_zone, b, 0.0)
            admm_penalty += λθ*(θ - θ̄) + (ρθ/2)*(θ - θ̄)^2
            admm_penalty += λV*(V - V̄) + (ρV/2)*(V - V̄)^2
        end
    end
    return admm_penalty
end

function add_penalty_branch_ac!(m, subdata, zone_edges::Vector{Int},
                                zP::Dict{Int,Float64}, zQ::Dict{Int,Float64},
                                dualsP_zone::Dict{Int,Float64}, dualsQ_zone::Dict{Int,Float64},
                                ρP::Float64, ρQ::Float64)
    @expression(m, admm_penalty, 0.0)
    for e in zone_edges
        l_local = findfirst(==(e), subdata["branch_local_indices"])
        @assert l_local !== nothing "aresta global $e não encontrada na partição"
        P_loc, Q_loc = _ac_flow_from_exprs_existing(m, subdata, l_local)
        λP = get(dualsP_zone, e, 0.0); λQ = get(dualsQ_zone, e, 0.0)
        admm_penalty += λP*(P_loc - zP[e]) + (ρP/2)*(P_loc - zP[e])^2
        admm_penalty += λQ*(Q_loc - zQ[e]) + (ρQ/2)*(Q_loc - zQ[e])^2
    end
    return admm_penalty
end

# ---------- Consensos ----------
function compute_ac_bus_consensus!(θ_avg::Dict{Int,Float64}, v_avg::Dict{Int,Float64},
                                   partitions, shared_buses::Vector{Int})
    for b in shared_buses
        thetas = Float64[]; Vacc = 0.0; n = 0
        for (_, sub) in partitions
            if b in sub["bus"]["bus_i"]
                m = sub["model"]
                push!(thetas, value(m[:theta][b,1]))
                Vacc += value(m[:v_magnitude][b,1])
                n += 1
            end
        end
        if n > 0
            θ_avg[b] = _circular_mean(thetas)
            v_avg[b] = Vacc / n
        end
    end
end

function compute_ac_branch_consensus!(zP::Dict{Int,Float64}, zQ::Dict{Int,Float64},
                                      partitions, cut_edges::Vector{Int})
    for e in cut_edges
        accP=0.0; accQ=0.0; n=0
        for (_, sub) in partitions
            l_local = findfirst(==(e), sub["branch_local_indices"])
            if l_local !== nothing
                P_loc, Q_loc = _ac_flow_from_exprs_existing(sub["model"], sub, l_local)
                accP += value(P_loc); accQ += value(Q_loc); n += 1
            end
        end
        if n>0
            zP[e] = accP/n; zQ[e] = accQ/n
        end
    end
end

# ---------- Atualização de duais ----------
"λθ_zb ← λθ_zb + ρθ * wrap(θ_zb - θ̄_b); λV_zb ← λV_zb + ρV (V_zb - V̄_b)."
function update_duals_bus_ac!(dualsθ, dualsV, partitions, θ_avg, v_avg, shared_buses, ρθ, ρV)
    for (z, sub) in partitions
        m = sub["model"]; bus_ids = sub["bus"]["bus_i"]
        for b in shared_buses
            if b in bus_ids
                θloc = value(m[:theta][b,1]); Vloc = value(m[:v_magnitude][b,1])
                dθ = _wrap_angle_diff(θloc - θ_avg[b]); dV = Vloc - v_avg[b]
                dualsθ[z][b] = get(dualsθ[z], b, 0.0) + ρθ*dθ
                dualsV[z][b] = get(dualsV[z], b, 0.0) + ρV*dV
            end
        end
    end
end

"λP_ze ← λP_ze + ρP (P_ze - zP_e); λQ_ze ← λQ_ze + ρQ (Q_ze - zQ_e)."
function update_duals_branch_ac!(dualsP, dualsQ, partitions, zP, zQ, cut_edges, ρP, ρQ)
    for (z, sub) in partitions
        for e in _zone_cut_edges(sub, cut_edges)
            l_local = findfirst(==(e), sub["branch_local_indices"])
            P_loc, Q_loc = _ac_flow_from_exprs_existing(sub["model"], sub, l_local)
            dualsP[z][e] = get(dualsP[z], e, 0.0) + ρP*(value(P_loc) - zP[e])
            dualsQ[z][e] = get(dualsQ[z], e, 0.0) + ρQ*(value(Q_loc) - zQ[e])
        end
    end
end

# ---------- Resíduos ----------
"‖r‖ (bus/AC) = sqrt(Σ_z Σ_b [wrap(θ_zb-θ̄_b)]² + (V_zb-V̄_b)²)."
function primal_residual_bus_ac(partitions, θ_avg, v_avg, shared_buses)
    acc = 0.0
    for (_, sub) in partitions
        m = sub["model"]; bus_ids = sub["bus"]["bus_i"]
        for b in shared_buses
            if b in bus_ids
                dθ = _wrap_angle_diff(value(m[:theta][b,1]) - θ_avg[b])
                dV = value(m[:v_magnitude][b,1]) - v_avg[b]
                acc += dθ^2 + dV^2
            end
        end
    end
    return sqrt(acc)
end

"‖s‖ (bus/AC) = sqrt(Σ_b wrap(θ̄_b-θ̄_b_prev)² + (V̄_b - V̄_b_prev)²)."
function dual_residual_bus_ac(θ_avg, v_avg, θ_prev, v_prev, shared_buses)
    sθ = sum(_wrap_angle_diff(θ_avg[b] - get(θ_prev, b, 0.0))^2 for b in shared_buses)
    sV = sum((v_avg[b] - get(v_prev, b, 0.0))^2 for b in shared_buses)
    return sqrt(sθ + sV)
end

function eps_primal_bus_ac(partitions, θ_avg, v_avg, shared_buses, ϵ_abs, ϵ_rel)
    dim = 2 * sum(b in partitions[z]["bus"]["bus_i"] for z in keys(partitions), b in shared_buses)
    z_norm = sqrt(sum(θ_avg[b]^2 for b in shared_buses) + sum(v_avg[b]^2 for b in shared_buses))
    return sqrt(max(1, dim))*ϵ_abs + ϵ_rel*max(z_norm, 1e-12)
end

function eps_dual_bus_ac(θ_avg, v_avg, θ_prev, v_prev, shared_buses, ϵ_abs, ϵ_rel)
    dim = 2 * length(shared_buses)
    sθ = sum(_wrap_angle_diff(θ_avg[b] - get(θ_prev, b, 0.0))^2 for b in shared_buses)
    sV = sum((v_avg[b] - get(v_prev, b, 0.0))^2 for b in shared_buses)
    return sqrt(max(1, dim))*ϵ_abs + ϵ_rel*sqrt(sθ + sV)
end

"‖r‖ (branch/AC) = sqrt(Σ_z Σ_e [(P_ze - zP_e)² + (Q_ze - zQ_e)²])."
function primal_residual_branch_ac(partitions, zP, zQ, cut_edges)
    acc = 0.0
    for (_, sub) in partitions
        for e in _zone_cut_edges(sub, cut_edges)
            l_local = findfirst(==(e), sub["branch_local_indices"])
            P_loc, Q_loc = _ac_flow_from_exprs_existing(sub["model"], sub, l_local)
            acc += (value(P_loc) - zP[e])^2 + (value(Q_loc) - zQ[e])^2
        end
    end
    return sqrt(acc)
end

function adaptive_rho!(ρ::Base.RefValue{Float64}, r::Float64, s::Float64;
                       μ_incr=2.0, μ_decr=2.0, τ_incr=1.05, τ_decr=0.98)
    if r > μ_incr*s
        ρ[] *= τ_incr
    elseif s > μ_decr*r
        ρ[] *= τ_decr
    end
end

"‖s‖ (branch/AC) = sqrt(Σ_e [(zP - zP_prev)² + (zQ - zQ_prev)²])."
function dual_residual_branch_ac(zP, zQ, zP_prev, zQ_prev, cut_edges)
    sP = sum((zP[e] - get(zP_prev, e, 0.0))^2 for e in cut_edges)
    sQ = sum((zQ[e] - get(zQ_prev, e, 0.0))^2 for e in cut_edges)
    return sqrt(sP + sQ)
end

function eps_primal_branch_ac(partitions, zP, zQ, cut_edges, ϵ_abs, ϵ_rel)
    dim_pri = 2 * sum(length(_zone_cut_edges(sub, cut_edges)) for (_, sub) in partitions)
    z_norm = sqrt(sum(zP[e]^2 + zQ[e]^2 for e in cut_edges))
    return sqrt(max(1, dim_pri))*ϵ_abs + ϵ_rel*max(z_norm, 1e-12)
end

function eps_dual_branch_ac(zP, zQ, zP_prev, zQ_prev, cut_edges, ϵ_abs, ϵ_rel)
    dim_du = 2 * length(cut_edges)
    s = dual_residual_branch_ac(zP, zQ, zP_prev, zQ_prev, cut_edges)
    return sqrt(max(1, dim_du))*ϵ_abs + ϵ_rel*s
end


# ---------- Plots (wrappers) ----------
function save_ac_bus_plots(res_pri, res_du, eps_pri, eps_du, rhoθ_hist, rhoV_hist, θ_avg_hist;
                           stem::String="admm_ac")
    plot_admm_convergence(res_pri, res_du, eps_pri, eps_du, "plot_convergencia_$(stem).html")
    PlotlyLight.save( (()->(p=Plot(); p(;x=1:length(rhoθ_hist), y=rhoθ_hist, type=:scatter, mode="lines+markers", name="ρθ");
        p.layout.title="Evolução de ρθ"; p.layout.xaxis.title="Iterações"; p.layout.yaxis.title="ρθ"; p.layout.yaxis.type="log"; p))(),
        "plot_rho_theta_$(stem).html")
    PlotlyLight.save( (()->(p=Plot(); p(;x=1:length(rhoV_hist), y=rhoV_hist, type=:scatter, mode="lines+markers", name="ρV");
        p.layout.title="Evolução de ρV"; p.layout.xaxis.title="Iterações"; p.layout.yaxis.title="ρV"; p.layout.yaxis.type="log"; p))(),
        "plot_rho_vmag_$(stem).html")
    plot_theta_avg_heatmap(θ_avg_hist, "plot_theta_avg_heatmap_$(stem).html")
end

function save_ac_branch_plots(res_pri, res_du, eps_pri, eps_du, rhoP_hist, rhoQ_hist;
                              stem::String="admm_ac_branch")
    plot_admm_convergence(res_pri, res_du, eps_pri, eps_du, "plot_convergencia_$(stem).html")
    PlotlyLight.save( (()->(p=Plot(); p(;x=1:length(rhoP_hist), y=rhoP_hist, type=:scatter, mode="lines+markers", name="ρP");
        p.layout.title="Evolução de ρP"; p.layout.xaxis.title="Iterações"; p.layout.yaxis.title="ρP"; p.layout.yaxis.type="log"; p))(),
        "plot_rhoP_$(stem).html")
    PlotlyLight.save( (()->(p=Plot(); p(;x=1:length(rhoQ_hist), y=rhoQ_hist, type=:scatter, mode="lines+markers", name="ρQ");
        p.layout.title="Evolução de ρQ"; p.layout.xaxis.title="Iterações"; p.layout.yaxis.title="ρQ"; p.layout.yaxis.type="log"; p))(),
        "plot_rhoQ_$(stem).html")
end

function run_admm_ac(config::DistributedOPFModel, data::Dict{String, Any})
    @assert config.base_model == AC "ADMM-AC: base_model deve ser AC."
    partitions, shared_buses = prepare_admm_partitions(data)

    # estado
    ρθ = Ref(config.rho_theta)
    ρV = Ref(config.rho_v)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual

    dualsθ = Dict{Int, Dict{Int, Float64}}(z => Dict{Int,Float64}() for z in keys(partitions))
    dualsV = Dict{Int, Dict{Int, Float64}}(z => Dict{Int,Float64}() for z in keys(partitions))
    θ_avg = Dict{Int, Float64}(); v_avg = Dict{Int, Float64}()
    θ_prev = Dict{Int, Float64}(); v_prev = Dict{Int, Float64}()

    # histórico
    res_pri_hist = Float64[]; res_du_hist = Float64[]
    eps_pri_hist = Float64[]; eps_du_hist = Float64[]
    rhoθ_hist = Float64[]; rhoV_hist = Float64[]
    θ_avg_hist = Dict{Int, Vector{Float64}}()

    optimizer_factory = resolve_optimizer(config)

    for k in 1:config.max_iter
        println("\n[ADMM-AC] Iteração $k")

        # 1) resolve locais com penalidade de BUS
        for (zone_id, sub) in partitions
            m = build_model(sub, AC; optimizer = optimizer_factory)
            admm_pen = add_penalty_bus_ac!(m, sub, shared_buses, θ_avg, v_avg, dualsθ[zone_id], dualsV[zone_id], ρθ[], ρV[])
            @objective(m, Min, m[:generation_cost] + admm_pen)
            optimize!(m)

            sub["model"] = m
            sub["generation_cost_model_expr"] = value(m[:generation_cost])
            sub["real_generation_cost"] = compute_real_generation_cost(m, sub["objective_params"])
        end

        # 2) consensos
        θ_prev = isempty(θ_avg) ? θ_avg : deepcopy(θ_avg)
        v_prev = isempty(v_avg) ? v_avg : deepcopy(v_avg)
        compute_ac_bus_consensus!(θ_avg, v_avg, partitions, shared_buses)
        for b in shared_buses
            push!(get!(θ_avg_hist, b, Float64[]), θ_avg[b])
        end

        # 3) duais
        update_duals_bus_ac!(dualsθ, dualsV, partitions, θ_avg, v_avg, shared_buses, ρθ[], ρV[])

        # 4) resíduos + tolerâncias
        r = primal_residual_bus_ac(partitions, θ_avg, v_avg, shared_buses)
        s = dual_residual_bus_ac(θ_avg, v_avg, θ_prev, v_prev, shared_buses)
        push!(res_pri_hist, r); push!(res_du_hist, s)

        ε_pri = eps_primal_bus_ac(partitions, θ_avg, v_avg, shared_buses, ϵ_abs, ϵ_rel)
        ε_du  = eps_dual_bus_ac(θ_avg, v_avg, θ_prev, v_prev, shared_buses, ϵ_abs, ϵ_rel)
        push!(eps_pri_hist, ε_pri); push!(eps_du_hist, ε_du)

        println("‖r‖ = $r | ε_pri = $ε_pri")
        println("‖s‖ = $s | ε_du  = $ε_du")

        push!(rhoθ_hist, ρθ[]); push!(rhoV_hist, ρV[])

        if r < ε_pri && s < ε_du
            println("Stopping criterion met.")
            break
        end

        adaptive_rho!(ρθ, r, s)
        adaptive_rho!(ρV, r, s)
    end

    save_ac_bus_plots(res_pri_hist, res_du_hist, eps_pri_hist, eps_du_hist, rhoθ_hist, rhoV_hist, θ_avg_hist; stem="admm_ac")

    total_real_cost = sum(p["real_generation_cost"] for p in values(partitions))
    total_expr_cost = sum(p["generation_cost_model_expr"] for p in values(partitions))
    println("\n[ADMM-AC] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("ρθ_final = $(ρθ[]) | ρV_final = $(ρV[])")

    return total_real_cost, partitions
end

function run_admm_ac_branch(config::DistributedOPFModel, data::Dict{String, Any})
    @assert config.base_model == AC "ADMM-AC-branch: base_model deve ser AC."
    partitions, _ = prepare_admm_partitions(data)
    cut_edges = _find_cut_edges(data)

    # estado
    ρP = Ref(config.rho_flow_p); ρQ = Ref(config.rho_flow_q)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual

    dualsP = Dict{Int, Dict{Int, Float64}}(z => Dict{Int,Float64}() for z in keys(partitions))
    dualsQ = Dict{Int, Dict{Int, Float64}}(z => Dict{Int,Float64}() for z in keys(partitions))
    zP = Dict{Int,Float64}(e => 0.0 for e in cut_edges)
    zQ = Dict{Int,Float64}(e => 0.0 for e in cut_edges)
    zP_prev = Dict{Int,Float64}(); zQ_prev = Dict{Int,Float64}()

    # histórico
    res_pri_hist = Float64[]; res_du_hist = Float64[]
    eps_pri_hist = Float64[]; eps_du_hist = Float64[]
    rhoP_hist = Float64[]; rhoQ_hist = Float64[]

    optimizer_factory = resolve_optimizer(config)

    for k in 1:config.max_iter
        println("\n[ADMM-AC-branch] Iteração $k")

        # 1) resolve locais com penalidade de BRANCH (P,Q)
        for (zone_id, sub) in partitions
            m = build_model(sub, AC; optimizer = optimizer_factory)
            zone_edges = _zone_cut_edges(sub, cut_edges)
            admm_pen = add_penalty_branch_ac!(m, sub, zone_edges, zP, zQ, dualsP[zone_id], dualsQ[zone_id], ρP[], ρQ[])
            @objective(m, Min, m[:generation_cost] + admm_pen)
            optimize!(m)

            sub["model"] = m
            sub["generation_cost_model_expr"] = value(m[:generation_cost])
            sub["real_generation_cost"] = compute_real_generation_cost(m, sub["objective_params"])
        end

        # 2) consensos (P,Q)
        zP_prev = isempty(zP) ? zP : deepcopy(zP)
        zQ_prev = isempty(zQ) ? zQ : deepcopy(zQ)
        compute_ac_branch_consensus!(zP, zQ, partitions, cut_edges)

        # 3) duais
        update_duals_branch_ac!(dualsP, dualsQ, partitions, zP, zQ, cut_edges, ρP[], ρQ[])

        # 4) resíduos + tolerâncias
        r = primal_residual_branch_ac(partitions, zP, zQ, cut_edges)
        s = dual_residual_branch_ac(zP, zQ, zP_prev, zQ_prev, cut_edges)
        push!(res_pri_hist, r); push!(res_du_hist, s)

        ε_pri = eps_primal_branch_ac(partitions, zP, zQ, cut_edges, ϵ_abs, ϵ_rel)
        ε_du  = eps_dual_branch_ac(zP, zQ, zP_prev, zQ_prev, cut_edges, ϵ_abs, ϵ_rel)
        push!(eps_pri_hist, ε_pri); push!(eps_du_hist, ε_du)

        println("‖r‖ = $r | ε_pri = $ε_pri")
        println("‖s‖ = $s | ε_du  = $ε_du")

        push!(rhoP_hist, ρP[]); push!(rhoQ_hist, ρQ[])

        if r < ε_pri && s < ε_du
            println("Stopping criterion met.")
            break
        end

        adaptive_rho!(ρP, r, s)
        adaptive_rho!(ρQ, r, s)
    end

    save_ac_branch_plots(res_pri_hist, res_du_hist, eps_pri_hist, eps_du_hist, rhoP_hist, rhoQ_hist; stem="admm_ac_branch")

    total_real_cost = sum(p["real_generation_cost"] for p in values(partitions))
    total_expr_cost = sum(p["generation_cost_model_expr"] for p in values(partitions))
    println("\n[ADMM-AC-branch] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("ρP_final = $(ρP[]) | ρQ_final = $(ρQ[])")

    return total_real_cost, partitions
end

function run_admm_ac_hybrid(config::DistributedOPFModel, data::Dict{String, Any})
    @assert config.base_model == AC "Hybrid AC: base_model deve ser AC."
    partitions, shared_buses = prepare_admm_partitions(data)
    cut_edges = _find_cut_edges(data)

    # estado
    ρθ = Ref(config.rho_theta); ρV = Ref(config.rho_v)
    ρP = Ref(config.rho_flow_p); ρQ = Ref(config.rho_flow_q)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual

    dualsθ = Dict{Int, Dict{Int, Float64}}(z => Dict{Int,Float64}() for z in keys(partitions))
    dualsV = Dict{Int, Dict{Int, Float64}}(z => Dict{Int,Float64}() for z in keys(partitions))
    dualsP = Dict{Int, Dict{Int, Float64}}(z => Dict{Int,Float64}() for z in keys(partitions))
    dualsQ = Dict{Int, Dict{Int, Float64}}(z => Dict{Int,Float64}() for z in keys(partitions))

    θ_avg = Dict{Int,Float64}(); v_avg = Dict{Int,Float64}()
    θ_prev = Dict{Int,Float64}(); v_prev = Dict{Int,Float64}()
    zP = Dict{Int,Float64}(e => 0.0 for e in cut_edges)
    zQ = Dict{Int,Float64}(e => 0.0 for e in cut_edges)
    zP_prev = Dict{Int,Float64}(); zQ_prev = Dict{Int,Float64}()

    # histórico
    res_pri_hist = Float64[]; res_du_hist = Float64[]
    eps_pri_hist = Float64[]; eps_du_hist = Float64[]
    rhoθ_hist = Float64[]; rhoV_hist = Float64[]; rhoP_hist = Float64[]; rhoQ_hist = Float64[]
    θ_avg_hist = Dict{Int, Vector{Float64}}()

    optimizer_factory = resolve_optimizer(config)

    for k in 1:config.max_iter
        println("\n[ADMM-AC-hybrid] Iteração $k")

        # 1) resolve locais com penalidades de BUS e BRANCH
        for (zone_id, sub) in partitions
            m = build_model(sub, AC; optimizer = optimizer_factory)
            @expression(m, admm_penalty, 0.0)

            # BUS (θ, V)
            admm_penalty += add_penalty_bus_ac!(m, sub, shared_buses, θ_avg, v_avg, dualsθ[zone_id], dualsV[zone_id], ρθ[], ρV[])
            # BRANCH (P, Q)
            zone_edges = _zone_cut_edges(sub, cut_edges)
            admm_penalty += add_penalty_branch_ac!(m, sub, zone_edges, zP, zQ, dualsP[zone_id], dualsQ[zone_id], ρP[], ρQ[])

            @objective(m, Min, m[:generation_cost] + admm_penalty)
            optimize!(m)

            sub["model"] = m
            sub["generation_cost_model_expr"] = value(m[:generation_cost])
            sub["real_generation_cost"] = compute_real_generation_cost(m, sub["objective_params"])
        end

        # 2) consensos
        θ_prev = isempty(θ_avg) ? θ_avg : deepcopy(θ_avg)
        v_prev = isempty(v_avg) ? v_avg : deepcopy(v_avg)
        zP_prev = isempty(zP) ? zP : deepcopy(zP)
        zQ_prev = isempty(zQ) ? zQ : deepcopy(zQ)

        compute_ac_bus_consensus!(θ_avg, v_avg, partitions, shared_buses)
        for b in shared_buses
            push!(get!(θ_avg_hist, b, Float64[]), θ_avg[b])
        end
        compute_ac_branch_consensus!(zP, zQ, partitions, cut_edges)

        # 3) duais
        update_duals_bus_ac!(dualsθ, dualsV, partitions, θ_avg, v_avg, shared_buses, ρθ[], ρV[])
        update_duals_branch_ac!(dualsP, dualsQ, partitions, zP, zQ, cut_edges, ρP[], ρQ[])

        # 4) resíduos + tolerâncias (agregados)
        r_bus = primal_residual_bus_ac(partitions, θ_avg, v_avg, shared_buses)
        r_br  = primal_residual_branch_ac(partitions, zP, zQ, cut_edges)
        r_tot = sqrt(r_bus^2 + r_br^2)

        s_bus = dual_residual_bus_ac(θ_avg, v_avg, θ_prev, v_prev, shared_buses)
        s_br  = dual_residual_branch_ac(zP, zQ, zP_prev, zQ_prev, cut_edges)
        s_tot = sqrt(s_bus^2 + s_br^2)

        push!(res_pri_hist, r_tot); push!(res_du_hist, s_tot)

        ε_pri = begin
            εb = eps_primal_bus_ac(partitions, θ_avg, v_avg, shared_buses, ϵ_abs, ϵ_rel)
            εe = eps_primal_branch_ac(partitions, zP, zQ, cut_edges, ϵ_abs, ϵ_rel)
            sqrt(εb^2 + εe^2)
        end
        ε_du = begin
            εb = eps_dual_bus_ac(θ_avg, v_avg, θ_prev, v_prev, shared_buses, ϵ_abs, ϵ_rel)
            εe = eps_dual_branch_ac(zP, zQ, zP_prev, zQ_prev, cut_edges, ϵ_abs, ϵ_rel)
            sqrt(εb^2 + εe^2)
        end
        push!(eps_pri_hist, ε_pri); push!(eps_du_hist, ε_du)

        println("‖r‖ = $r_tot | ε_pri = $ε_pri")
        println("‖s‖ = $s_tot | ε_du  = $ε_du")

        push!(rhoθ_hist, ρθ[]); push!(rhoV_hist, ρV[]); push!(rhoP_hist, ρP[]); push!(rhoQ_hist, ρQ[])

        if r_tot < ε_pri && s_tot < ε_du
            println("Stopping criterion met.")
            break
        end

        # ρ adaptativo em blocos (mantendo tua lógica)
        adaptive_rho!(ρθ, r_bus, s_bus)
        adaptive_rho!(ρV, r_bus, s_bus)
        adaptive_rho!(ρP, r_br,  s_br)
        adaptive_rho!(ρQ, r_br,  s_br)
    end

    plot_admm_convergence(res_pri_hist, res_du_hist, eps_pri_hist, eps_du_hist,
        "plot_convergencia_admm_ac_hybrid.html")

    PlotlyLight.save( (()->(p=Plot(); p(;x=1:length(rhoθ_hist), y=rhoθ_hist, type=:scatter, mode="lines+markers", name="ρθ");
        p.layout.title="Evolução de ρθ"; p.layout.xaxis.title="Iterações"; p.layout.yaxis.title="ρθ"; p.layout.yaxis.type="log"; p))(),
        "plot_rho_theta_ac_hybrid.html")
    PlotlyLight.save( (()->(p=Plot(); p(;x=1:length(rhoV_hist), y=rhoV_hist, type=:scatter, mode="lines+markers", name="ρV");
        p.layout.title="Evolução de ρV"; p.layout.xaxis.title="Iterações"; p.layout.yaxis.title="ρV"; p.layout.yaxis.type="log"; p))(),
        "plot_rho_v_ac_hybrid.html")
    PlotlyLight.save( (()->(p=Plot(); p(;x=1:length(rhoP_hist), y=rhoP_hist, type=:scatter, mode="lines+markers", name="ρP");
        p.layout.title="Evolução de ρP"; p.layout.xaxis.title="Iterações"; p.layout.yaxis.title="ρP"; p.layout.yaxis.type="log"; p))(),
        "plot_rho_p_ac_hybrid.html")
    PlotlyLight.save( (()->(p=Plot(); p(;x=1:length(rhoQ_hist), y=rhoQ_hist, type=:scatter, mode="lines+markers", name="ρQ");
        p.layout.title="Evolução de ρQ"; p.layout.xaxis.title="Iterações"; p.layout.yaxis.title="ρQ"; p.layout.yaxis.type="log"; p))(),
        "plot_rho_q_ac_hybrid.html")

    plot_theta_avg_heatmap(θ_avg_hist, "plot_theta_avg_heatmap_ac_hybrid.html")

    total_real_cost = sum(p["real_generation_cost"] for p in values(partitions))
    total_expr_cost = sum(p["generation_cost_model_expr"] for p in values(partitions))
    println("\n[ADMM-AC-hybrid] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("ρθ_final = $(ρθ[]) | ρV_final = $(ρV[]) | ρP_final = $(ρP[]) | ρQ_final = $(ρQ[])")

    return total_real_cost, partitions
end
