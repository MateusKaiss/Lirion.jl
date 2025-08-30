# ================= AC Utils =================
# ---------- Solver ----------
function resolve_optimizer_ac(config)
    if hasproperty(config, :optimizer) && getproperty(config, :optimizer) !== nothing
        return getproperty(config, :optimizer)
    else
        return Ipopt.Optimizer
    end
end

# ---------- Helpers ----------
@inline function _edge_local(sub::Dict{String,Any}, e::Int)
    edgemap = get(sub, "edge_to_local", nothing)
    return edgemap === nothing ? findfirst(==(e), sub["branch_local_indices"]) :
                                 get(edgemap, e, nothing)
end

@inline function _bus_set(sub::Dict{String,Any})
    get(sub, "bus_set", Set(Int[b for b in sub["bus"]["bus_i"]]))
end

@inline function _assert_nonzero(x::Real, msg::AbstractString)
    @assert abs(x) > 1e-8 msg
end

# ---------- Penalidades (ADMM) ----------

"Penalidade por BARRA no AC: θ (cópia local) vs θ̄ e V (cópia local) vs V̄.
 Retorna AffExpr ANÔNIMO (não registra nome no modelo)."
function add_penalty_bus_ac!(
    m::JuMP.Model,
    subdata::Dict{String,Any},
    shared_buses::AbstractVector{<:Integer},
    θ_avg::AbstractDict{<:Integer,<:Real},
    v_avg::AbstractDict{<:Integer,<:Real},
    dualsθ_zone::AbstractDict{<:Integer,<:Real},
    dualsV_zone::AbstractDict{<:Integer,<:Real},
    ρθ::Real,
    ρV::Real,
)
    pen = @expression(m, 0.0)  # anônima → sem colisão de nome
    busset = _bus_set(subdata)
    @inbounds for b in shared_buses
        if b in busset
            θ   = m[:theta][b, 1]
            θ̄  = Float64(get(θ_avg, b, 0.0))
            λθ  = Float64(get(dualsθ_zone, b, 0.0))
            pen += λθ * (θ - θ̄) + (float(ρθ) / 2) * (θ - θ̄)^2

            V   = m[:v_magnitude][b, 1]
            V̄  = Float64(get(v_avg, b, 1.0))
            λV  = Float64(get(dualsV_zone, b, 0.0))
            pen += λV * (V - V̄) + (float(ρV) / 2) * (V - V̄)^2
        end
    end
    return pen
end

"Penalidade por RAMO no AC: P_from/Q_from vs zP/zQ.
 Retorna AffExpr ANÔNIMO."
function add_penalty_branch_ac!(
    m::JuMP.Model,
    subdata::Dict{String,Any},
    zone_edges::AbstractVector{<:Integer},
    zP::AbstractDict{<:Integer,<:Real},
    zQ::AbstractDict{<:Integer,<:Real},
    dualsP_zone::AbstractDict{<:Integer,<:Real},
    dualsQ_zone::AbstractDict{<:Integer,<:Real},
    ρP::Real,
    ρQ::Real,
)
    pen = @expression(m, 0.0)  # anônima
    @inbounds for e in zone_edges
        l_local = _edge_local(subdata, e)
        @assert l_local !== nothing "aresta global $e não pertence a esta partição"
        # usa função existente que retorna as EXPRESSÕES JuMP de P_from/Q_from
        P_loc, Q_loc = _ac_flow_from_exprs_existing(m, subdata, l_local)
        λP = Float64(get(dualsP_zone, e, 0.0))
        λQ = Float64(get(dualsQ_zone, e, 0.0))
        zPe = Float64(get(zP, e, 0.0))
        zQe = Float64(get(zQ, e, 0.0))
        pen += λP * (P_loc - zPe) + (float(ρP) / 2) * (P_loc - zPe)^2
        pen += λQ * (Q_loc - zQe) + (float(ρQ) / 2) * (Q_loc - zQe)^2
    end
    return pen
end

# ---------- Consensos ----------
function compute_ac_bus_consensus!(θ_avg::Dict{Int,Float64}, v_avg::Dict{Int,Float64},
                                   partitions, shared_buses::Vector{Int})
    # cacheia value uma vez por modelo
    thetas_map = Dict{Int,Vector{Float64}}()  # acumula θ por bus (para média circular)
    Vacc = Dict{Int,Float64}()               # soma de V por bus
    counts = Dict{Int,Int}()
    for (_, sub) in partitions
        m = sub["model"]
        θvals = value.(m[:theta])
        Vvals = value.(m[:v_magnitude])
        busset = _bus_set(sub)
        @inbounds for b in shared_buses
            if b in busset
                push!(get!(thetas_map, b, Float64[]), θvals[b,1])
                Vacc[b] = get(Vacc,b,0.0) + Vvals[b,1]
                counts[b] = get(counts,b,0) + 1
            end
        end
    end
    @inbounds for b in shared_buses
        n = get(counts,b,0)
        if n > 0
            θ_avg[b] = _circular_mean(get(thetas_map,b,Float64[]))
            v_avg[b] = Vacc[b] / n
        end
    end
end

function compute_ac_branch_consensus!(zP::Dict{Int,Float64}, zQ::Dict{Int,Float64},
                                      partitions, cut_edges::Vector{Int})
    sumsP = Dict{Int,Float64}(); sumsQ = Dict{Int,Float64}(); counts = Dict{Int,Int}()
    for (_, sub) in partitions
        m = sub["model"]
        @inbounds for e in cut_edges
            l_local = _edge_local(sub, e)
            if l_local !== nothing
                P_loc, Q_loc = _ac_flow_from_exprs_existing(m, sub, l_local)
                # value apenas aqui (1x)
                p = value(P_loc); q = value(Q_loc)
                sumsP[e] = get(sumsP,e,0.0) + p
                sumsQ[e] = get(sumsQ,e,0.0) + q
                counts[e] = get(counts,e,0) + 1
            end
        end
    end
    @inbounds for e in cut_edges
        n = get(counts,e,0)
        if n > 0
            zP[e] = sumsP[e] / n
            zQ[e] = sumsQ[e] / n
        end
    end
end

# ---------- Atualização de duais ----------
"λθ_zb ← λθ_zb + ρθ * wrap(θ_zb - θ̄_b); λV_zb ← λV_zb + ρV (V_zb - V̄_b)."
function update_duals_bus_ac!(dualsθ, dualsV, partitions, θ_avg, v_avg, shared_buses, ρθ, ρV)
    for (z, sub) in partitions
        m = sub["model"]
        θvals = value.(m[:theta])
        Vvals = value.(m[:v_magnitude])
        busset = _bus_set(sub)
        @inbounds for b in shared_buses
            if b in busset
                dθ = _wrap_angle_diff(θvals[b,1] - get(θ_avg,b,0.0))
                dV = Vvals[b,1] - get(v_avg,b,0.0)
                dualsθ[z][b] = get(dualsθ[z], b, 0.0) + ρθ*dθ
                dualsV[z][b] = get(dualsV[z], b, 0.0) + ρV*dV
            end
        end
    end
end

"λP_ze ← λP_ze + ρP (P_ze - zP_e); λQ_ze ← λQ_ze + ρQ (Q_ze - zQ_e)."
function update_duals_branch_ac!(dualsP, dualsQ, partitions, zP, zQ, cut_edges, ρP, ρQ)
    for (z, sub) in partitions
        m = sub["model"]
        @inbounds for e in cut_edges
            l_local = _edge_local(sub, e)
            if l_local !== nothing
                P_loc, Q_loc = _ac_flow_from_exprs_existing(m, sub, l_local)
                dualsP[z][e] = get(dualsP[z], e, 0.0) + ρP*(value(P_loc) - zP[e])
                dualsQ[z][e] = get(dualsQ[z], e, 0.0) + ρQ*(value(Q_loc) - zQ[e])
            end
        end
    end
end

# ---------- Resíduos ----------
"‖r‖ (bus/AC) = sqrt(Σ_z Σ_b [wrap(θ_zb-θ̄_b)]² + (V_zb-V̄_b)²)."
function primal_residual_bus_ac(partitions, θ_avg, v_avg, shared_buses)
    acc = 0.0
    for (_, sub) in partitions
        m = sub["model"]
        θvals = value.(m[:theta])
        Vvals = value.(m[:v_magnitude])
        busset = _bus_set(sub)
        @inbounds for b in shared_buses
            if b in busset
                dθ = _wrap_angle_diff(θvals[b,1] - get(θ_avg,b,0.0))
                dV = Vvals[b,1] - get(v_avg,b,0.0)
                acc += dθ^2 + dV^2
            end
        end
    end
    return sqrt(acc)
end

"‖s‖ (bus/AC) = sqrt(Σ_b wrap(θ̄_b-θ̄_b_prev)² + (V̄_b - V̄_b_prev)²)."
function dual_residual_bus_ac(θ_avg, v_avg, θ_prev, v_prev, shared_buses)
    sθ = 0.0; sV = 0.0
    @inbounds for b in shared_buses
        sθ += _wrap_angle_diff(get(θ_avg,b,0.0) - get(θ_prev,b,0.0))^2
        sV += (get(v_avg,b,0.0) - get(v_prev,b,0.0))^2
    end
    return sqrt(sθ + sV)
end

function eps_primal_bus_ac(partitions, θ_avg, v_avg, shared_buses, ϵ_abs, ϵ_rel)
    # dim = 2 * número de cópias locais de buses compartilhadas
    dim = 0
    for (_, sub) in partitions
        busset = _bus_set(sub)
        @inbounds for b in shared_buses
            dim += (b in busset) ? 2 : 0
        end
    end
    # norma de z = concat(θ̄, V̄)
    z_norm_sq = 0.0
    @inbounds for b in shared_buses
        z_norm_sq += get(θ_avg,b,0.0)^2 + get(v_avg,b,0.0)^2
    end
    return sqrt(max(1, dim))*ϵ_abs + ϵ_rel*max(sqrt(z_norm_sq), 1e-12)
end

function eps_dual_bus_ac(θ_avg, v_avg, θ_prev, v_prev, shared_buses, ϵ_abs, ϵ_rel)
    dim = 2 * max(1, length(shared_buses))
    sθ = 0.0; sV = 0.0
    @inbounds for b in shared_buses
        sθ += _wrap_angle_diff(get(θ_avg,b,0.0) - get(θ_prev,b,0.0))^2
        sV += (get(v_avg,b,0.0) - get(v_prev,b,0.0))^2
    end
    return sqrt(dim)*ϵ_abs + ϵ_rel*sqrt(sθ + sV)
end

"‖r‖ (branch/AC) = sqrt(Σ_z Σ_e [(P_ze - zP_e)² + (Q_ze - zQ_e)²])."
function primal_residual_branch_ac(partitions, zP, zQ, cut_edges)
    acc = 0.0
    for (_, sub) in partitions
        m = sub["model"]
        @inbounds for e in cut_edges
            l_local = _edge_local(sub, e)
            if l_local !== nothing
                P_loc, Q_loc = _ac_flow_from_exprs_existing(m, sub, l_local)
                acc += (value(P_loc) - zP[e])^2 + (value(Q_loc) - zQ[e])^2
            end
        end
    end
    return sqrt(acc)
end

function dual_residual_branch_ac(zP, zQ, zP_prev, zQ_prev, cut_edges)
    sP = 0.0; sQ = 0.0
    @inbounds for e in cut_edges
        sP += (get(zP,e,0.0) - get(zP_prev,e,0.0))^2
        sQ += (get(zQ,e,0.0) - get(zQ_prev,e,0.0))^2
    end
    return sqrt(sP + sQ)
end

function eps_primal_branch_ac(partitions, zP, zQ, cut_edges, ϵ_abs, ϵ_rel)
    dim_pri = 0
    for (_, sub) in partitions
        dim_pri += 2 * length(_zone_cut_edges(sub, cut_edges))
    end
    z_norm_sq = 0.0
    @inbounds for e in cut_edges
        z_norm_sq += get(zP,e,0.0)^2 + get(zQ,e,0.0)^2
    end
    return sqrt(max(1, dim_pri))*ϵ_abs + ϵ_rel*max(sqrt(z_norm_sq), 1e-12)
end

function eps_dual_branch_ac(zP, zQ, zP_prev, zQ_prev, cut_edges, ϵ_abs, ϵ_rel)
    dim_du = 2 * max(1, length(cut_edges))
    s = dual_residual_branch_ac(zP, zQ, zP_prev, zQ_prev, cut_edges)
    return sqrt(dim_du)*ϵ_abs + ϵ_rel*s
end

# ---------- ρ adaptativo ----------
function adaptive_rho_ac!(ρ::Base.RefValue{Float64}, r::Float64, s::Float64;
                          μ_incr=2.0, μ_decr=2.0, τ_incr=1.05, τ_decr=0.98)
    if r > μ_incr*s
        ρ[] *= τ_incr
    elseif s > μ_decr*r
        ρ[] *= τ_decr
    end
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

    optimizer_factory = resolve_optimizer_ac(config)

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

        adaptive_rho_ac!(ρθ, r, s)
        adaptive_rho_ac!(ρV, r, s)
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

    optimizer_factory = resolve_optimizer_ac(config)

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

        adaptive_rho_ac!(ρP, r, s)
        adaptive_rho_ac!(ρQ, r, s)
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

    optimizer_factory = resolve_optimizer_ac(config)

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
        adaptive_rho_ac!(ρθ, r_bus, s_bus)
        adaptive_rho_ac!(ρV, r_bus, s_bus)
        adaptive_rho_ac!(ρP, r_br,  s_br)
        adaptive_rho_ac!(ρQ, r_br,  s_br)
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
