# =============== DC Utils ===============
# ---------- Solver ----------
function resolve_optimizer(config)
    if hasproperty(config, :optimizer) && getproperty(config, :optimizer) !== nothing
        return getproperty(config, :optimizer)
    else
        return Ipopt.Optimizer
    end
end

# ---------- Penalidades (ADMM) ----------
function add_penalty_bus_dc!(m, subdata, shared_buses::Vector{Int},
                             θ_avg::Dict{Int,Float64}, duals_zone::Dict{Int,Float64}, ρ::Float64)
    bus_ids = subdata["bus"]["bus_i"]
    @expression(m, admm_penalty, 0.0)
    for b in shared_buses
        if b in bus_ids
            θ  = m[:theta][b,1]
            θ̄ = get(θ_avg, b, 0.0)
            λ  = get(duals_zone, b, 0.0)
            admm_penalty += λ*(θ - θ̄) + (ρ/2)*(θ - θ̄)^2
        end
    end
    return admm_penalty
end

function add_penalty_branch_dc!(m, subdata, zone_edges::Vector{Int},
                                z_flow::Dict{Int,Float64}, duals_zone::Dict{Int,Float64}, ρ::Float64)
    @expression(m, admm_penalty, 0.0)
    for e in zone_edges
        l_local = findfirst(==(e), subdata["branch_local_indices"])
        p_loc   = _dc_flow_expr(m, subdata, l_local)
        λ       = get(duals_zone, e, 0.0)
        z       = z_flow[e]
        admm_penalty += λ*(p_loc - z) + (ρ/2)*(p_loc - z)^2
    end
    return admm_penalty
end

# ---------- Consensos ----------
function compute_theta_consensus!(θ_avg::Dict{Int,Float64}, partitions, shared_buses::Vector{Int})
    for b in shared_buses
        acc = 0.0; n = 0
        for (_, subdata) in partitions
            if b in subdata["bus"]["bus_i"]
                acc += value(subdata["model"][:theta][b,1]); n += 1
            end
        end
        θ_avg[b] = (n > 0) ? (acc / n) : get(θ_avg, b, 0.0)
    end
end

function compute_z_flow_consensus!(z_flow::Dict{Int,Float64}, partitions, cut_edges::Vector{Int})
    for e in cut_edges
        acc = 0.0; n = 0
        for (_, subdata) in partitions
            local_idx = findfirst(==(e), subdata["branch_local_indices"])
            if local_idx !== nothing
                acc += value(_dc_flow_expr(subdata["model"], subdata, local_idx))
                n += 1
            end
        end
        z_flow[e] = (n > 0) ? acc/n : get(z_flow, e, 0.0)
    end
end

# ---------- Atualização de duais ----------
function update_duals_theta!(duals::Dict{Int,Dict{Int,Float64}}, partitions,
                             θ_avg::Dict{Int,Float64}, shared_buses::Vector{Int}, ρ::Float64)
    for (z, subdata) in partitions
        for b in shared_buses
            if b in subdata["bus"]["bus_i"]
                θb = value(subdata["model"][:theta][b,1])
                duals[z][b] = get(duals[z], b, 0.0) + ρ*(θb - θ_avg[b])
            end
        end
    end
end

"λe_ze ← λe_ze + ρ (p_ze - z_e)"
function update_duals_flow!(duals::Dict{Int,Dict{Int,Float64}}, partitions,
                            z_flow::Dict{Int,Float64}, cut_edges::Vector{Int}, ρ::Float64)
    for (z, subdata) in partitions
        zone_edges = _zone_cut_edges(subdata, cut_edges)
        for e in zone_edges
            l_local = findfirst(==(e), subdata["branch_local_indices"])
            p_val = value(_dc_flow_expr(subdata["model"], subdata, l_local))
            duals[z][e] = get(duals[z], e, 0.0) + ρ*(p_val - z_flow[e])
        end
    end
end

# ---------- Resíduos ----------
"‖r‖ (bus): sqrt(Σ_z Σ_b (θ_zb - θ̄_b)^2)."
function primal_residual_bus_dc(partitions, θ_avg::Dict{Int,Float64}, shared_buses::Vector{Int})
    acc = 0.0
    for (_, subdata) in partitions
        for b in shared_buses
            if b in subdata["bus"]["bus_i"]
                acc += (value(subdata["model"][:theta][b,1]) - θ_avg[b])^2
            end
        end
    end
    return sqrt(acc)
end

"‖s‖ (bus): sqrt(Σ_b (θ̄_b - θ̄_b_prev)^2)."
dual_residual_bus_dc(θ_avg::Dict{Int,Float64}, θ_avg_prev::Dict{Int,Float64}) =
    sqrt(sum((θ_avg[b] - get(θ_avg_prev, b, 0.0))^2 for b in keys(θ_avg)))

"ε_pri (bus) = √N ε_abs + ε_rel max(‖θ_local‖, ‖θ̄‖)."
function eps_primal_bus_dc(partitions, θ_avg::Dict{Int,Float64}, shared_buses::Vector{Int},
                           ϵ_abs::Float64, ϵ_rel::Float64)
    N = length(shared_buses)
    norm_x_sq = 0.0
    for (_, subdata) in partitions
        for b in shared_buses
            if b in subdata["bus"]["bus_i"]
                norm_x_sq += (value(subdata["model"][:theta][b,1]))^2
            end
        end
    end
    norm_z_sq = sum((θ_avg[b])^2 for b in shared_buses)
    return sqrt(N)*ϵ_abs + ϵ_rel*max(sqrt(norm_x_sq), sqrt(norm_z_sq))
end

"ε_du (bus) = √N ε_abs + ε_rel ‖θ̄ - θ̄_prev‖."
function eps_dual_bus_dc(θ_avg::Dict{Int,Float64}, θ_avg_prev::Dict{Int,Float64},
                         ϵ_abs::Float64, ϵ_rel::Float64)
    N = max(1, length(θ_avg))
    return sqrt(N)*ϵ_abs + ϵ_rel*sqrt(sum((θ_avg[b] - get(θ_avg_prev, b, 0.0))^2 for b in keys(θ_avg)))
end

"‖r‖ (branch): sqrt(Σ_z Σ_e (p_ze - z_e)^2)."
function primal_residual_branch_dc(partitions, z_flow::Dict{Int,Float64}, cut_edges::Vector{Int})
    acc = 0.0
    for (_, subdata) in partitions
        for e in _zone_cut_edges(subdata, cut_edges)
            l_local = findfirst(==(e), subdata["branch_local_indices"])
            acc += (value(_dc_flow_expr(subdata["model"], subdata, l_local)) - z_flow[e])^2
        end
    end
    return sqrt(acc)
end

"‖s‖ (branch): sqrt(Σ_e (z_e - z_e_prev)^2)."
dual_residual_branch_dc(z_flow::Dict{Int,Float64}, z_flow_prev::Dict{Int,Float64}) =
    sqrt(sum((z_flow[e] - get(z_flow_prev, e, 0.0))^2 for e in keys(z_flow)))

"ε_pri (branch) = √dim_pri ε_abs + ε_rel max(‖z‖, 1e-12)."
function eps_primal_branch_dc(partitions, z_flow::Dict{Int,Float64}, cut_edges::Vector{Int},
                              ϵ_abs::Float64, ϵ_rel::Float64)
    dim_pri = sum(length(_zone_cut_edges(sub, cut_edges)) for (_, sub) in partitions)
    return sqrt(max(1, dim_pri))*ϵ_abs + ϵ_rel*max(sqrt(sum(z_flow[e]^2 for e in cut_edges)), 1e-12)
end

"ε_du (branch) = √dim_du ε_abs + ε_rel ‖z - z_prev‖."
function eps_dual_branch_dc(z_flow::Dict{Int,Float64}, z_flow_prev::Dict{Int,Float64},
                            ϵ_abs::Float64, ϵ_rel::Float64)
    dim_du = max(1, length(z_flow))
    return sqrt(dim_du)*ϵ_abs + ϵ_rel*sqrt(sum((z_flow[e] - get(z_flow_prev, e, 0.0))^2 for e in keys(z_flow)))
end

# ---------- ρ adaptativo ----------
function adaptive_rho!(ρ::Base.RefValue{Float64}, r::Float64, s::Float64;
                       μ_incr=2.0, μ_decr=2.0, τ_incr=1.05, τ_decr=0.98)
    if r > μ_incr*s
        ρ[] *= τ_incr
    elseif s > μ_decr*r
        ρ[] *= τ_decr
    end
end

# ---------- Plots (wrappers) ----------
function save_dc_bus_plots(res_pri, res_du, eps_pri, eps_du, ρ_hist, θ_avg_history;
                           stem::String="admm")
    plot_admm_convergence(res_pri, res_du, eps_pri, eps_du,
        "plot_convergencia_$(stem)_dc.html")
    plot_rho(ρ_hist, "plot_rho_$(stem)_dc.html")
    plot_theta_avg_heatmap(θ_avg_history, "plot_theta_avg_heatmap_$(stem)_dc.html")
end

function save_dc_branch_plots(res_pri, res_du, eps_pri, eps_du, ρ_hist;
                              stem::String="admm_branch")
    plot_admm_convergence(res_pri, res_du, eps_pri, eps_du,
        "plot_convergencia_$(stem).html")
    plot_rho(ρ_hist, "plot_rho_$(stem).html")
end

function run_admm(config::DistributedOPFModel, data)
    partitions, shared_buses = prepare_admm_partitions(data)
    base_model = config.base_model
    base_model == DC || error("run_admm (bus-consensus) está implementado para DC.")

    # --- estado ADMM
    ρ = Ref(config.rho_theta)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual
    duals = Dict{Int, Dict{Int, Float64}}(z => Dict{Int,Float64}() for z in keys(partitions))
    θ_avg = Dict{Int,Float64}(); θ_avg_prev = Dict{Int,Float64}()

    # --- histórico
    res_primal_vec = Float64[]; res_dual_vec = Float64[]
    eps_primal_vec = Float64[]; eps_dual_vec = Float64[]
    rho_vec = Float64[]; θ_avg_history = Dict{Int, Vector{Float64}}()
    θ_local_history = Dict{Int, Dict{Int, Vector{Float64}}}()

    optimizer_factory = resolve_optimizer(config)

    for k in 1:config.max_iter
        println("\nIteração $k (bus/DC)")

        # 1) resolve subproblemas locais + penalidade
        for (zone_id, subdata) in partitions
            local_model = build_model(subdata, base_model; optimizer = optimizer_factory)

            # penalidade por barra (θ)
            admm_penalty = add_penalty_bus_dc!(local_model, subdata, shared_buses, θ_avg, duals[zone_id], ρ[])
            @objective(local_model, Min, local_model[:generation_cost] + admm_penalty)
            optimize!(local_model)

            subdata["model"] = local_model
            subdata["generation_cost_model_expr"] = value(local_model[:generation_cost])
            subdata["real_generation_cost"] = compute_real_generation_cost(local_model, subdata["objective_params"])

            # histórico local de θ por barra compartilhada
            θ_vals = value.(local_model[:theta])
            for b in shared_buses
                if b in subdata["bus"]["bus_i"]
                    push!(get!(get!(θ_local_history, zone_id, Dict{Int,Vector{Float64}}()), b, Float64[]),
                          θ_vals[b,1])
                end
            end
        end

        # 2) consenso
        θ_avg_prev = isempty(θ_avg) ? θ_avg : deepcopy(θ_avg)
        compute_theta_consensus!(θ_avg, partitions, shared_buses)
        for b in shared_buses
            push!(get!(θ_avg_history, b, Float64[]), θ_avg[b])
        end

        # 3) duais
        update_duals_theta!(duals, partitions, θ_avg, shared_buses, ρ[])

        # 4) resíduos + tolerâncias
        res_primal = primal_residual_bus_dc(partitions, θ_avg, shared_buses)
        res_dual   = dual_residual_bus_dc(θ_avg, θ_avg_prev)
        push!(res_primal_vec, res_primal); push!(res_dual_vec, res_dual)

        ϵ_pri = eps_primal_bus_dc(partitions, θ_avg, shared_buses, ϵ_abs, ϵ_rel)
        ϵ_du  = eps_dual_bus_dc(θ_avg, θ_avg_prev, ϵ_abs, ϵ_rel)
        push!(eps_primal_vec, ϵ_pri); push!(eps_dual_vec, ϵ_du)

        println("‖r‖ = $res_primal | ε_pri = $ϵ_pri")
        println("‖s‖ = $res_dual   | ε_dual = $ϵ_du")

        push!(rho_vec, ρ[])

        # 5) checagem de parada
        if res_primal < ϵ_pri && res_dual < ϵ_du
            println("Stopping criterion met.")
            break
        end

        # 6) ρ adaptativo
        adaptive_rho!(ρ, res_primal, res_dual)
    end

    # plots
    save_dc_bus_plots(res_primal_vec, res_dual_vec, eps_primal_vec, eps_dual_vec, rho_vec, θ_avg_history; stem="admm")

    total_real_cost = sum(sub["real_generation_cost"] for sub in values(partitions))
    total_expr_cost = sum(sub["generation_cost_model_expr"] for sub in values(partitions))

    println("\nResumo (bus/DC):")
    println("→ Total (value(:generation_cost)) = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("ρθ inicial: $(config.rho_theta) | ρθ final: $(ρ[])")

    return total_real_cost, partitions
end

function run_admm_branch(config::DistributedOPFModel, data::Dict{String, Any})
    partitions, _ = prepare_admm_partitions(data)
    base_model = config.base_model
    base_model == DC || error("ADMM-branch está implementado aqui apenas para DC.")
    cut_edges = _find_cut_edges(data)

    # estado ADMM
    ρ = Ref(config.rho_flow)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual
    duals = Dict{Int, Dict{Int, Float64}}(z => Dict{Int,Float64}() for z in keys(partitions))
    z_flow = Dict{Int,Float64}(e => 0.0 for e in cut_edges)
    z_flow_prev = Dict{Int,Float64}()

    # histórico
    res_primal_vec = Float64[]; res_dual_vec = Float64[]
    eps_primal_vec = Float64[]; eps_dual_vec = Float64[]
    rho_vec = Float64[]

    optimizer_factory = resolve_optimizer(config)

    for k in 1:config.max_iter
        println("\n[ADMM-branch/DC] Iteração $k")

        # 1) resolve locais + penalidade de fluxo
        for (zone_id, subdata) in partitions
            local_model = build_model(subdata, base_model; optimizer = optimizer_factory)
            zone_edges = _zone_cut_edges(subdata, cut_edges)
            admm_penalty = add_penalty_branch_dc!(local_model, subdata, zone_edges, z_flow, duals[zone_id], ρ[])
            @objective(local_model, Min, local_model[:generation_cost] + admm_penalty)
            optimize!(local_model)

            subdata["model"] = local_model
            subdata["generation_cost_model_expr"] = value(local_model[:generation_cost])
            subdata["real_generation_cost"] = compute_real_generation_cost(local_model, subdata["objective_params"])
        end

        # 2) consenso de fluxo
        z_flow_prev = isempty(z_flow) ? z_flow : deepcopy(z_flow)
        compute_z_flow_consensus!(z_flow, partitions, cut_edges)

        # 3) duais
        update_duals_flow!(duals, partitions, z_flow, cut_edges, ρ[])

        # 4) resíduos + tolerâncias
        res_primal = primal_residual_branch_dc(partitions, z_flow, cut_edges)
        res_dual   = dual_residual_branch_dc(z_flow, z_flow_prev)
        push!(res_primal_vec, res_primal); push!(res_dual_vec, res_dual)

        ϵ_pri = eps_primal_branch_dc(partitions, z_flow, cut_edges, ϵ_abs, ϵ_rel)
        ϵ_du  = eps_dual_branch_dc(z_flow, z_flow_prev, ϵ_abs, ϵ_rel)
        push!(eps_primal_vec, ϵ_pri); push!(eps_dual_vec, ϵ_du)

        println("‖r‖ = $res_primal | ε_pri = $ϵ_pri")
        println("‖s‖ = $res_dual   | ε_dual = $ϵ_du")

        push!(rho_vec, ρ[])

        if res_primal < ϵ_pri && res_dual < ϵ_du
            println("Stopping criterion met.")
            break
        end

        adaptive_rho!(ρ, res_primal, res_dual)
    end

    save_dc_branch_plots(res_primal_vec, res_dual_vec, eps_primal_vec, eps_dual_vec, rho_vec; stem="admm_branch")

    total_real_cost = sum(p["real_generation_cost"] for p in values(partitions))
    total_expr_cost = sum(p["generation_cost_model_expr"] for p in values(partitions))

    println("\n[ADMM-branch/DC] Resumo:")
    println("→ Total (value(:generation_cost)) = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("ρe inicial: $(config.rho_flow) | ρe final: $(ρ[])")

    return total_real_cost, partitions
end

function run_admm_hybrid(config::DistributedOPFModel, data::Dict{String, Any})
    @assert config.base_model == DC "Hybrid ADMM currently implemented for DC only."

    partitions, shared_buses = prepare_admm_partitions(data)
    cut_edges = _find_cut_edges(data)

    # estado ADMM
    ρθ = Ref(config.rho_theta)
    ρe = Ref(config.rho_flow)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual

    θ_avg = Dict{Int,Float64}(); θ_avg_prev = Dict{Int,Float64}()
    z_flow = Dict{Int,Float64}(e => 0.0 for e in cut_edges); z_flow_prev = Dict{Int,Float64}()

    duals_theta = Dict{Int, Dict{Int, Float64}}(z => Dict{Int,Float64}() for z in keys(partitions))
    duals_flow  = Dict{Int, Dict{Int, Float64}}(z => Dict{Int,Float64}() for z in keys(partitions))

    # histórico
    res_primal_vec = Float64[]; res_dual_vec = Float64[]
    eps_primal_vec = Float64[]; eps_dual_vec = Float64[]
    rhoθ_vec = Float64[]; rhoe_vec = Float64[]
    θ_avg_history = Dict{Int, Vector{Float64}}()

    optimizer_factory = resolve_optimizer(config)

    for k in 1:config.max_iter
        println("\n[ADMM-hybrid/DC] Iteração $k")

        # 1) resolve locais + penalidades (θ e fluxo)
        for (zone_id, subdata) in partitions
            local_model = build_model(subdata, DC; optimizer = optimizer_factory)

            @expression(local_model, admm_penalty, 0.0)
            # θ
            admm_penalty += add_penalty_bus_dc!(local_model, subdata, shared_buses, θ_avg, duals_theta[zone_id], ρθ[])
            # fluxo
            zone_edges = _zone_cut_edges(subdata, cut_edges)
            admm_penalty += add_penalty_branch_dc!(local_model, subdata, zone_edges, z_flow, duals_flow[zone_id], ρe[])

            @objective(local_model, Min, local_model[:generation_cost] + admm_penalty)
            optimize!(local_model)

            subdata["model"] = local_model
            subdata["generation_cost_model_expr"] = value(local_model[:generation_cost])
            subdata["real_generation_cost"] = compute_real_generation_cost(local_model, subdata["objective_params"])
        end

        # 2) consensos
        θ_avg_prev = isempty(θ_avg) ? θ_avg : deepcopy(θ_avg)
        z_flow_prev = deepcopy(z_flow)

        compute_theta_consensus!(θ_avg, partitions, shared_buses)
        for b in shared_buses
            push!(get!(θ_avg_history, b, Float64[]), θ_avg[b])
        end
        compute_z_flow_consensus!(z_flow, partitions, cut_edges)

        # 3) duais
        update_duals_theta!(duals_theta, partitions, θ_avg, shared_buses, ρθ[])
        update_duals_flow!(duals_flow, partitions, z_flow, cut_edges, ρe[])

        # 4) resíduos + tolerâncias (agregados)
        res_pri_bus = primal_residual_bus_dc(partitions, θ_avg, shared_buses)
        res_pri_br  = primal_residual_branch_dc(partitions, z_flow, cut_edges)
        res_primal  = sqrt(res_pri_bus^2 + res_pri_br^2)

        res_du_bus = dual_residual_bus_dc(θ_avg, θ_avg_prev)
        res_du_br  = dual_residual_branch_dc(z_flow, z_flow_prev)
        res_dual   = sqrt(res_du_bus^2 + res_du_br^2)

        push!(res_primal_vec, res_primal); push!(res_dual_vec, res_dual)

        # dimensões p/ eps_primal (soma das dimensões de bus e branch)
        ϵ_pri = begin
            ϵb = eps_primal_bus_dc(partitions, θ_avg, shared_buses, ϵ_abs, ϵ_rel)
            ϵe = eps_primal_branch_dc(partitions, z_flow, cut_edges, ϵ_abs, ϵ_rel)
            # heurística simples: combina como norma euclidiana
            sqrt(ϵb^2 + ϵe^2)
        end
        ϵ_du = begin
            ϵb = eps_dual_bus_dc(θ_avg, θ_avg_prev, ϵ_abs, ϵ_rel)
            ϵe = eps_dual_branch_dc(z_flow, z_flow_prev, ϵ_abs, ϵ_rel)
            sqrt(ϵb^2 + ϵe^2)
        end
        push!(eps_primal_vec, ϵ_pri); push!(eps_dual_vec, ϵ_du)

        println("‖r‖ = $res_primal | ε_pri = $ϵ_pri")
        println("‖s‖ = $res_dual   | ε_dual = $ϵ_du")

        push!(rhoθ_vec, ρθ[]); push!(rhoe_vec, ρe[])

        if res_primal < ϵ_pri && res_dual < ϵ_du
            println("Stopping criterion met.")
            break
        end

        # Atualiza AMBOS os ρ para manter escalas parecidas
        adaptive_rho!(ρθ, res_primal, res_dual)
        adaptive_rho!(ρe, res_primal, res_dual)
    end

    # Plots (reusa wrappers + dois gráficos separados de ρ)
    plot_admm_convergence(res_primal_vec, res_dual_vec, eps_primal_vec, eps_dual_vec,
        "plot_convergencia_admm_hybrid.html")

    PlotlyLight.save( (()->(p=Plot(); p(;x=1:length(rhoθ_vec), y=rhoθ_vec, type=:scatter, mode="lines+markers", name="ρθ"); p.layout.title="Evolução de ρθ"; p.layout.xaxis.title="Iterações"; p.layout.yaxis.title="ρθ"; p.layout.yaxis.type="log"; p))(),
        "plot_rho_theta_hybrid.html")
    PlotlyLight.save( (()->(p=Plot(); p(;x=1:length(rhoe_vec), y=rhoe_vec, type=:scatter, mode="lines+markers", name="ρe"); p.layout.title="Evolução de ρe"; p.layout.xaxis.title="Iterações"; p.layout.yaxis.title="ρe"; p.layout.yaxis.type="log"; p))(),
        "plot_rho_flow_hybrid.html")

    plot_theta_avg_heatmap(θ_avg_history, "plot_theta_avg_heatmap_hybrid.html")

    total_real_cost = sum(p["real_generation_cost"] for p in values(partitions))
    total_expr_cost = sum(p["generation_cost_model_expr"] for p in values(partitions))

    println("\n[ADMM-hybrid/DC] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("ρθ_final = $(ρθ[]) | ρe_final = $(ρe[])")

    return total_real_cost, partitions
end