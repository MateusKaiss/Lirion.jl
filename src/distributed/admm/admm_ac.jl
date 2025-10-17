# ================= AC Utils =================
# ---------- Solver ----------
function resolve_optimizer_ac(config)
    if hasproperty(config, :optimizer) && getproperty(config, :optimizer) !== nothing
        return getproperty(config, :optimizer)
    else
        return Ipopt.Optimizer
    end
end

@inline function _ac_flow_from_exprs_existing(
    model::Model,
    subdata::Dict{String, Any},
    l_local::Int;
    t::Int = 1,
)
    return model[:p_flow_branches][l_local, t], model[:q_flow_branches][l_local, t]
end

@inline function _wrap_angle_diff(Δ::Float64)
    return atan(sin(Δ), cos(Δ))
end

@inline function _circular_mean(θs::AbstractVector{<:Real})
    sum_sin = 0.0
    sum_cos = 0.0
    @inbounds @simd for θ in θs
        sum_sin += sin(θ)
        sum_cos += cos(θ)
    end
    return atan(sum_sin, sum_cos)
end

@inline function _pretty_time_ac(sec::Real)
    if sec < 60
        return "$(round(sec; digits=3)) s"
    elseif sec < 3600
        return "$(round(sec/60; digits=2)) min"
    else
        return "$(round(sec/3600; digits=2)) h"
    end
end

@inline function _val_or_start_ac(m::JuMP.Model, x::JuMP.VariableRef)
    if JuMP.has_values(m)
        return JuMP.value(x)
    else
        sv = JuMP.start_value(x)
        return sv === nothing ? 0.0 : sv
    end
end

@inline function _theta_val_ac(
    m::JuMP.Model,
    θ::JuMP.Containers.DenseAxisArray,
    b::Int,
    t::Int,
)
    return _val_or_start_ac(m, θ[b, t])
end

@inline function _vmag_val(m::JuMP.Model, V::JuMP.Containers.DenseAxisArray, b::Int, t::Int)
    return _val_or_start_ac(m, V[b, t])
end

# ---------- Hotpaths (buffers e mapeamentos densos) ----------
function prepare_hotpaths_ac!(sub::Dict{String, Any}, shared_buses::Vector{Int})
    bus_ids = sub["bus"]["bus_i"]
    bus_local_index = Dict{Int, Int}()
    @inbounds for (i, b) in enumerate(bus_ids)
        bus_local_index[b] = i
    end
    sub["bus_local_index"] = bus_local_index

    busset = Set(bus_ids)
    local_shared = [b for b in shared_buses if b in busset]
    sub["shared_buses_local"] = local_shared

    shared_pos_local = Vector{Int}(undef, length(local_shared))
    @inbounds for (k, b) in enumerate(local_shared)
        shared_pos_local[k] = bus_local_index[b]
    end
    sub["shared_pos_local"] = shared_pos_local

    sub["theta_buf"] = zeros(Float64, length(bus_ids))
    sub["vmag_buf"] = ones(Float64, length(bus_ids))
    sub["theta_shared_buf"] = zeros(Float64, length(shared_pos_local))
    sub["vmag_shared_buf"] = ones(Float64, length(shared_pos_local))
    return nothing
end

@inline function snapshot_theta_ac!(sub::Dict{String, Any})
    m = sub["model"]
    θ = m[:theta]
    buf = sub["theta_buf"]
    bus_ids = sub["bus"]["bus_i"]
    @inbounds for i in eachindex(bus_ids)
        buf[i] = _theta_val_ac(m, θ, bus_ids[i], 1)
    end
    return buf
end

@inline function snapshot_vmag!(sub::Dict{String, Any})
    m = sub["model"]
    V = m[:v_magnitude]
    buf = sub["vmag_buf"]
    bus_ids = sub["bus"]["bus_i"]
    @inbounds for i in eachindex(bus_ids)
        buf[i] = _vmag_val(m, V, bus_ids[i], 1)
    end
    return buf
end

@inline function snapshot_theta_ac_shared!(sub::Dict{String, Any})
    m = sub["model"]
    θ = m[:theta]
    buf = sub["theta_shared_buf"]
    pos = sub["shared_pos_local"]
    bus_ids = sub["bus"]["bus_i"]
    @inbounds for i in eachindex(pos)
        bidx = pos[i]
        buf[i] = _theta_val_ac(m, θ, bus_ids[bidx], 1)
    end
    return buf
end

@inline function snapshot_vmag_shared!(sub::Dict{String, Any})
    m = sub["model"]
    V = m[:v_magnitude]
    buf = sub["vmag_shared_buf"]
    pos = sub["shared_pos_local"]
    bus_ids = sub["bus"]["bus_i"]
    @inbounds for i in eachindex(pos)
        bidx = pos[i]
        buf[i] = _vmag_val(m, V, bus_ids[bidx], 1)
    end
    return buf
end

# ---------- Maps densos globais para shared buses ----------
function prepare_dense_maps_ac(shared_buses::Vector{Int})
    shared_buses_sorted = sort(shared_buses)
    b2pos = Dict{Int, Int}(b => i for (i, b) in enumerate(shared_buses_sorted))
    return shared_buses_sorted, b2pos
end

# ---------- Penalidades (pré-constroem termos quadráticos e só trocam parte linear) ----------
function _build_quad_bus_expr_ac!(sub::Dict{String, Any}, ρθ::Float64, ρV::Float64)
    m = sub["model"]
    θ = m[:theta]
    V = m[:v_magnitude]
    q = zero(JuMP.QuadExpr)
    for b in sub["shared_buses_local"]
        q += 0.5 * ρθ * θ[b, 1]^2 + 0.5 * ρV * V[b, 1]^2
    end
    sub["quad_bus_expr_ac"] = q
    sub["ρ ακ_bus_cached"] = (ρθ, ρV)
    return q
end

function _build_linear_bus_expr_ac!(
    sub::Dict{String, Any},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
    θ̄_vec::Vector{Float64},
    V̄_vec::Vector{Float64},
    λθ_zone::AbstractVector{<:Real},
    λV_zone::AbstractVector{<:Real},
    ρθ::Float64,
    ρV::Float64,
)
    m = sub["model"]
    θ = m[:theta]
    V = m[:v_magnitude]
    lin = zero(JuMP.AffExpr)
    for b in sub["shared_buses_local"]
        j = b2pos[b]
        lin += (λθ_zone[j] - ρθ * θ̄_vec[j]) * θ[b, 1]
        lin += (λV_zone[j] - ρV * V̄_vec[j]) * V[b, 1]
    end
    sub["lin_bus_expr_ac"] = lin
    return lin
end

function _build_quad_branch_expr_ac!(
    sub::Dict{String, Any}, ρP::Float64, ρQ::Float64, zone_edges::Vector{Int},
)
    m = sub["model"]
    q = zero(JuMP.QuadExpr)
    for e in zone_edges
        l_local = get(sub["edge_to_local"], e, nothing)
        l_local === nothing && continue
        P_loc, Q_loc = _ac_flow_from_exprs_existing(m, sub, l_local)
        q += 0.5 * ρP * P_loc^2 + 0.5 * ρQ * Q_loc^2
    end
    sub["quad_branch_expr_ac"] = q
    sub["ρ ακ_branch_cached"] = (ρP, ρQ)
    return q
end

function _build_linear_branch_expr_ac!(
    sub::Dict{String, Any},
    zone_edges::Vector{Int},
    zP::Dict{Int, Float64},
    zQ::Dict{Int, Float64},
    dualsP_zone::Dict{Int, Float64},
    dualsQ_zone::Dict{Int, Float64},
    ρP::Float64,
    ρQ::Float64,
)
    m = sub["model"]
    lin = zero(JuMP.AffExpr)
    for e in zone_edges
        l_local = get(sub["edge_to_local"], e, nothing)
        l_local === nothing && continue
        P_loc, Q_loc = _ac_flow_from_exprs_existing(m, sub, l_local)
        λP = get(dualsP_zone, e, 0.0)
        λQ = get(dualsQ_zone, e, 0.0)
        zPe = get(zP, e, 0.0)
        zQe = get(zQ, e, 0.0)
        lin += λP * (P_loc - zPe) + λQ * (Q_loc - zQe) - ρP * zPe * P_loc - ρQ * zQe * Q_loc
        # o termo 0.5*ρ*(P^2/Q^2) está no quadrático pré-construído,
        # aqui só ficam os termos lineares dependentes de z/λ.
    end
    sub["lin_branch_expr_ac"] = lin
    return lin
end

# ---------- Consensos densos (bus AC: θ por média circular + |V| média aritmética) ----------
function compute_ac_bus_consensus_vec!(
    θ̄_vec::Vector{Float64},
    V̄_vec::Vector{Float64},
    partitions::Dict{Int, Dict{String, Any}},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
)
    fill!(θ̄_vec, 0.0)
    fill!(V̄_vec, 0.0)
    sums_sin = zeros(Float64, length(shared_buses_sorted))
    sums_cos = zeros(Float64, length(shared_buses_sorted))
    counts = zeros(Int, length(shared_buses_sorted))

    for (_, sub) in partitions
        θs = snapshot_theta_ac_shared!(sub)
        Vs = snapshot_vmag_shared!(sub)
        lbs = sub["shared_buses_local"]
        @inbounds for i in eachindex(lbs)
            j = b2pos[lbs[i]]
            s, c = sincos(θs[i])
            sums_sin[j] += s
            sums_cos[j] += c
            V̄_vec[j] += Vs[i]
            counts[j] += 1
        end
    end

    @inbounds for j in eachindex(θ̄_vec)
        if counts[j] != 0
            θ̄_vec[j] = atan(sums_sin[j], sums_cos[j])
            V̄_vec[j] /= counts[j]
        end
    end
end

# ---------- Consenso para BRANCH (AC) ----------
function compute_ac_branch_consensus!(
    zP::Dict{Int, Float64},
    zQ::Dict{Int, Float64},
    partitions::Dict{Int, Dict{String, Any}},
    cut_edges::Vector{Int},
)
    # acumula P/Q locais das zonas e faz média por aresta de corte
    sumsP = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    sumsQ = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    counts = Dict{Int, Int}(e => 0 for e in cut_edges)

    for (_, sub) in partitions
        m = sub["model"]
        # use as arestas desta zona já pré-filtradas
        for e in get(sub, "zone_edges", Int[])
            l_local = get(sub["edge_to_local"], e, nothing)
            l_local === nothing && continue
            P_loc, Q_loc = _ac_flow_from_exprs_existing(m, sub, l_local)
            pval = JuMP.value(P_loc)
            qval = JuMP.value(Q_loc)
            # protege contra NaN/Inf ocasionais em iterações iniciais
            if isfinite(pval)
                sumsP[e] += pval
                counts[e] += 1
            end
            if isfinite(qval)
                sumsQ[e] += qval
            end
        end
    end

    for e in cut_edges
        c = counts[e]
        if c > 0
            zP[e] = sumsP[e] / c
            zQ[e] = sumsQ[e] / c
        end
        # se c==0, mantém zP/zQ como estavam
    end
    return nothing
end

# ---------- Duais densos para BUS (vetorizados por zona) ----------
@inline function update_duals_bus_ac_vec!(
    dualsθ_vec_zone::Vector{Float64},
    dualsV_vec_zone::Vector{Float64},
    sub::Dict{String, Any},
    θ̄_vec::Vector{Float64},
    V̄_vec::Vector{Float64},
    ρθ::Float64,
    ρV::Float64,
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
)
    θs = snapshot_theta_ac_shared!(sub)
    Vs = snapshot_vmag_shared!(sub)
    lbs = sub["shared_buses_local"]
    @inbounds for i in eachindex(lbs)
        j = b2pos[lbs[i]]
        dθ = _wrap_angle_diff(θs[i] - θ̄_vec[j])
        dV = (Vs[i] - V̄_vec[j])
        dualsθ_vec_zone[j] += ρθ * dθ
        dualsV_vec_zone[j] += ρV * dV
    end
    return nothing
end

# ---------- Duais para BRANCH (AC) ----------
function update_duals_branch_ac!(
    dualsP::Dict{Int, Dict{Int, Float64}},
    dualsQ::Dict{Int, Dict{Int, Float64}},
    partitions::Dict{Int, Dict{String, Any}},
    zP::Dict{Int, Float64},
    zQ::Dict{Int, Float64},
    cut_edges::Vector{Int},   # não é estritamente necessário aqui, mas mantém a assinatura
    ρP::Float64,
    ρQ::Float64,
)
    for (z, sub) in partitions
        m = sub["model"]
        λPz = get!(dualsP, z, Dict{Int, Float64}())
        λQz = get!(dualsQ, z, Dict{Int, Float64}())

        for e in get(sub, "zone_edges", Int[])
            l_local = get(sub["edge_to_local"], e, nothing)
            l_local === nothing && continue
            P_loc, Q_loc = _ac_flow_from_exprs_existing(m, sub, l_local)
            pval = JuMP.value(P_loc)
            qval = JuMP.value(Q_loc)

            if isfinite(pval)
                λPz[e] = get(λPz, e, 0.0) + ρP * (pval - get(zP, e, 0.0))
            end
            if isfinite(qval)
                λQz[e] = get(λQz, e, 0.0) + ρQ * (qval - get(zQ, e, 0.0))
            end
        end
    end
    return nothing
end

# ---------- Resíduos densos para BUS ----------
function primal_residual_bus_ac_vec(
    partitions::Dict{Int, Dict{String, Any}},
    θ̄_vec::Vector{Float64},
    V̄_vec::Vector{Float64},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
)
    acc = 0.0
    for (_, sub) in partitions
        θbuf = snapshot_theta_ac!(sub)
        Vbuf = snapshot_vmag!(sub)
        pos = sub["bus_local_index"]
        @inbounds for b in sub["shared_buses_local"]
            j = b2pos[b]
            acc +=
                _wrap_angle_diff(θbuf[pos[b]] - θ̄_vec[j])^2 + (Vbuf[pos[b]] - V̄_vec[j])^2
        end
    end
    return sqrt(acc)
end

function eps_primal_bus_ac_vec(
    partitions::Dict{Int, Dict{String, Any}},
    θ̄_vec::Vector{Float64},
    V̄_vec::Vector{Float64},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
    ϵ_abs::Float64,
    ϵ_rel::Float64,
)
    dim = 0
    for (_, sub) in partitions
        dim += 2 * length(sub["shared_buses_local"])
    end
    z_norm_sq = 0.0
    @inbounds for j in eachindex(shared_buses_sorted)
        z_norm_sq += θ̄_vec[j]^2 + V̄_vec[j]^2
    end
    return sqrt(max(1, dim)) * ϵ_abs + ϵ_rel * max(sqrt(z_norm_sq), 1e-12)
end

function dual_residual_bus_ac_vec(
    θ̄_vec::Vector{Float64}, V̄_vec::Vector{Float64},
    θ̄_prev::Vector{Float64}, V̄_prev::Vector{Float64},
)
    sθ = 0.0
    sV = 0.0
    @inbounds for j in eachindex(θ̄_vec)
        sθ += _wrap_angle_diff(θ̄_vec[j] - θ̄_prev[j])^2
        sV += (V̄_vec[j] - V̄_prev[j])^2
    end
    return sqrt(sθ + sV)
end

function eps_dual_bus_ac_vec(
    θ̄_vec::Vector{Float64}, V̄_vec::Vector{Float64},
    θ̄_prev::Vector{Float64}, V̄_prev::Vector{Float64},
    ϵ_abs::Float64, ϵ_rel::Float64,
)
    dim = 2 * max(1, length(θ̄_vec))
    s = dual_residual_bus_ac_vec(θ̄_vec, V̄_vec, θ̄_prev, V̄_prev)
    return sqrt(dim) * ϵ_abs + ϵ_rel * s
end

# ---------- ρ adaptativo ----------
function adaptive_rho_ac!(ρ::Base.RefValue{Float64}, r::Float64, s::Float64;
    μ_incr = 2.0, μ_decr = 2.0, τ_incr = 1.05, τ_decr = 0.98)
    if r > μ_incr * s
        ρ[] *= τ_incr
    elseif s > μ_decr * r
        ρ[] *= τ_decr
    end
end

# ---------- Resíduos para BRANCH (AC) ----------
function primal_residual_branch_ac(
    partitions::Dict{Int, Dict{String, Any}},
    zP::Dict{Int, Float64},
    zQ::Dict{Int, Float64},
    cut_edges::Vector{Int},
)
    acc = 0.0
    for (_, sub) in partitions
        m = sub["model"]
        for e in get(sub, "zone_edges", Int[])
            l_local = get(sub["edge_to_local"], e, nothing)
            l_local === nothing && continue
            P_loc, Q_loc = _ac_flow_from_exprs_existing(m, sub, l_local)
            p = JuMP.value(P_loc)
            q = JuMP.value(Q_loc)
            if isfinite(p)
                acc += (p - get(zP, e, 0.0))^2
            end
            if isfinite(q)
                acc += (q - get(zQ, e, 0.0))^2
            end
        end
    end
    return sqrt(acc)
end

function dual_residual_branch_ac(
    zP::Dict{Int, Float64},
    zQ::Dict{Int, Float64},
    zP_prev::Dict{Int, Float64},
    zQ_prev::Dict{Int, Float64},
    cut_edges::Vector{Int},
)
    s = 0.0
    for e in cut_edges
        s += (get(zP, e, 0.0) - get(zP_prev, e, 0.0))^2
        s += (get(zQ, e, 0.0) - get(zQ_prev, e, 0.0))^2
    end
    return sqrt(s)
end

function eps_primal_branch_ac(
    partitions::Dict{Int, Dict{String, Any}},
    zP::Dict{Int, Float64},
    zQ::Dict{Int, Float64},
    cut_edges::Vector{Int},
    ϵ_abs::Float64,
    ϵ_rel::Float64,
)
    # dimensão = 2 variáveis (P,Q) por aresta local em cada zona
    dim = 0
    for (_, sub) in partitions
        dim += 2 * length(get(sub, "zone_edges", Int[]))
    end
    # norma de z para termo relativo
    z_norm_sq = 0.0
    for e in cut_edges
        z_norm_sq += get(zP, e, 0.0)^2 + get(zQ, e, 0.0)^2
    end
    return sqrt(max(1, dim)) * ϵ_abs + ϵ_rel * max(sqrt(z_norm_sq), 1e-12)
end

function eps_dual_branch_ac(
    zP::Dict{Int, Float64},
    zQ::Dict{Int, Float64},
    zP_prev::Dict{Int, Float64},
    zQ_prev::Dict{Int, Float64},
    cut_edges::Vector{Int},
    ϵ_abs::Float64,
    ϵ_rel::Float64,
)
    dim = 2 * max(1, length(cut_edges))
    s = dual_residual_branch_ac(zP, zQ, zP_prev, zQ_prev, cut_edges)
    return sqrt(dim) * ϵ_abs + ϵ_rel * s
end

# ---------- Plots (mesma API) ----------
function save_ac_bus_plots(res_pri, res_du, eps_pri, eps_du, rhoθ_hist, rhoV_hist,
    θ_avg_hist;
    stem::String = "admm_ac")
    plot_admm_convergence(
        res_pri,
        res_du,
        eps_pri,
        eps_du,
        "plot_convergencia_$(stem).html",
    )
    PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(;
                    x = 1:length(rhoθ_hist),
                    y = rhoθ_hist,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρθ",
                );
                p.layout.title = "Evolução de ρθ";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "ρθ";
                p.layout.yaxis.type = "log";
                p
            )
        )(), "plot_rho_theta_$(stem).html")
    PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(;
                    x = 1:length(rhoV_hist),
                    y = rhoV_hist,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρV",
                );
                p.layout.title = "Evolução de ρV";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "ρV";
                p.layout.yaxis.type = "log";
                p
            )
        )(), "plot_rho_vmag_$(stem).html")
    return plot_theta_avg_heatmap(θ_avg_hist, "plot_theta_avg_heatmap_$(stem).html")
end

function save_ac_branch_plots(res_pri, res_du, eps_pri, eps_du, rhoP_hist, rhoQ_hist;
    stem::String = "admm_ac_branch")
    plot_admm_convergence(
        res_pri,
        res_du,
        eps_pri,
        eps_du,
        "plot_convergencia_$(stem).html",
    )
    PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(;
                    x = 1:length(rhoP_hist),
                    y = rhoP_hist,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρP",
                );
                p.layout.title = "Evolução de ρP";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "ρP";
                p.layout.yaxis.type = "log";
                p
            )
        )(), "plot_rhoP_$(stem).html")
    return PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(;
                    x = 1:length(rhoQ_hist),
                    y = rhoQ_hist,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρQ",
                );
                p.layout.title = "Evolução de ρQ";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "ρQ";
                p.layout.yaxis.type = "log";
                p
            )
        )(), "plot_rhoQ_$(stem).html")
end

# ============================ ADMM: BUS (AC) ============================
function run_admm_ac(config::DistributedOPFModel, data::Dict{String, Any})
    t0 = time_ns()
    @assert config.base_model == AC "ADMM-AC: base_model deve ser AC."
    partitions, shared_buses = prepare_admm_partitions(data)

    ρθ = Ref(config.rho_theta)
    ρV = Ref(config.rho_v)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual

    shared_buses_sorted, b2pos = prepare_dense_maps_ac(shared_buses)

    # duais densos por zona
    dualsθ_dense = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )
    dualsV_dense = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )

    θ̄_vec = zeros(Float64, length(shared_buses_sorted))
    V̄_vec = ones(Float64, length(shared_buses_sorted))
    θ̄_prev = similar(θ̄_vec)
    V̄_prev = similar(V̄_vec)

    res_pri_hist = Float64[]
    res_du_hist = Float64[]
    eps_pri_hist = Float64[]
    eps_du_hist = Float64[]
    rhoθ_hist = Float64[]
    rhoV_hist = Float64[]
    θ_avg_hist = Dict{Int, Vector{Float64}}()

    optimizer_factory = resolve_optimizer_ac(config)

    # build inicial + hotpaths + pré-quad
    for (zone_id, sub) in partitions
        m = build_model(sub, AC; optimizer = optimizer_factory, build_only = true)
        sub["model"] = m
        sub["obj_base"] = m[:generation_cost]
        sub["theta_var"] = m[:theta]
        sub["vmag_var"] = m[:v_magnitude]
        sub["has_solution"] = false
        _set_model_silent!(m)

        # start values leves
        @inbounds for I in keys(sub["theta_var"])
            set_start_value(sub["theta_var"][I], 0.0)
        end
        @inbounds for I in keys(sub["vmag_var"])
            set_start_value(sub["vmag_var"][I], 1.0)
        end

        prepare_hotpaths_ac!(sub, shared_buses)
        _build_quad_bus_expr_ac!(sub, ρθ[], ρV[])
        @objective(m, Min, sub["obj_base"] + sub["quad_bus_expr_ac"])
    end

    prev_ρθ, prev_ρV = ρθ[], ρV[]
    zone_ids = collect(keys(partitions))
    # buffers para juntar λ por zona (evita alocar toda iteração)
    λθ_zone_buf = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in zone_ids
    )
    λV_zone_buf = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in zone_ids
    )

    for k in 1:config.max_iter
        _FAST_USE_LOGS && println("\n[ADMM-AC] Iteração $k")

        # se rho mudou, reconstrói apenas termo quadrático
        if ρθ[] != prev_ρθ || ρV[] != prev_ρV
            for (_, sub) in partitions
                _build_quad_bus_expr_ac!(sub, ρθ[], ρV[])
            end
            prev_ρθ, prev_ρV = ρθ[], ρV[]
        end

        # alvo do termo linear: θ̄_vec, V̄_vec reais
        # (são atualizados depois da otimização de zona)
        # solve por zona em paralelo
        nz = length(zone_ids)
        results = Vector{Any}(undef, nz)
        Threads.@threads for i in 1:nz
            z = zone_ids[i]
            sub = partitions[z]
            # junta λ densos da zona
            λθ = (λθ_zone_buf[z] .= dualsθ_dense[z])
            λV = (λV_zone_buf[z] .= dualsV_dense[z])
            lin = _build_linear_bus_expr_ac!(sub, shared_buses_sorted, b2pos,
                θ̄_vec, V̄_vec, λθ, λV, ρθ[], ρV[])
            @objective(sub["model"], Min, sub["obj_base"] + sub["quad_bus_expr_ac"] + lin)
            optimize!(sub["model"])
            m = sub["model"]
            results[i] = (expr_cost = value(m[:generation_cost]),
                real_cost = compute_real_generation_cost(m, sub["objective_params"]))
        end

        for (i, z) in enumerate(zone_ids)
            sub = partitions[z]
            out = results[i]
            sub["has_solution"] = true
            sub["generation_cost_model_expr"] = out.expr_cost
            sub["real_generation_cost"] = out.real_cost
        end

        # consenso (densos)
        θ̄_prev .= θ̄_vec
        V̄_prev .= V̄_vec
        compute_ac_bus_consensus_vec!(
            θ̄_vec,
            V̄_vec,
            partitions,
            shared_buses_sorted,
            b2pos,
        )
        @inbounds for b in shared_buses_sorted
            j = b2pos[b]
            push!(get!(θ_avg_hist, b, Float64[]), θ̄_vec[j])
        end

        # duais (densos)
        for z in zone_ids
            update_duals_bus_ac_vec!(dualsθ_dense[z], dualsV_dense[z], partitions[z],
                θ̄_vec, V̄_vec, ρθ[], ρV[], shared_buses_sorted, b2pos)
        end

        # resíduos/tolerâncias
        r = primal_residual_bus_ac_vec(
            partitions,
            θ̄_vec,
            V̄_vec,
            shared_buses_sorted,
            b2pos,
        )
        s = dual_residual_bus_ac_vec(θ̄_vec, V̄_vec, θ̄_prev, V̄_prev)
        εp = eps_primal_bus_ac_vec(
            partitions,
            θ̄_vec,
            V̄_vec,
            shared_buses_sorted,
            b2pos,
            ϵ_abs,
            ϵ_rel,
        )
        εd = eps_dual_bus_ac_vec(θ̄_vec, V̄_vec, θ̄_prev, V̄_prev, ϵ_abs, ϵ_rel)

        push!(res_pri_hist, r)
        push!(res_du_hist, s)
        push!(eps_pri_hist, εp)
        push!(eps_du_hist, εd)
        _FAST_USE_LOGS && println("‖r‖ = $r | ε_pri = $εp")
        _FAST_USE_LOGS && println("‖s‖ = $s | ε_du  = $εd")

        push!(rhoθ_hist, ρθ[])
        push!(rhoV_hist, ρV[])

        if r < εp && s < εd
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        adaptive_rho_ac!(ρθ, r, s)
        adaptive_rho_ac!(ρV, r, s)
    end

    save_ac_bus_plots(
        res_pri_hist,
        res_du_hist,
        eps_pri_hist,
        eps_du_hist,
        rhoθ_hist,
        rhoV_hist,
        θ_avg_hist;
        stem = "admm_ac",
    )

    total_real_cost = sum(p["real_generation_cost"] for p in values(partitions))
    total_expr_cost = sum(p["generation_cost_model_expr"] for p in values(partitions))
    println("\n[ADMM-AC] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("ρθ_final = $(ρθ[]) | ρV_final = $(ρV[])")
    elapsed_s = (time_ns() - t0) / 1e9
    println("Tempo ADMM (bus/AC): $(_pretty_time_ac(elapsed_s))")

    return total_real_cost, partitions
end

# ========================== ADMM: BRANCH (AC) ==========================
function run_admm_ac_branch(config::DistributedOPFModel, data::Dict{String, Any})
    t0 = time_ns()
    @assert config.base_model == AC "ADMM-AC-branch: base_model deve ser AC."
    partitions, _ = prepare_admm_partitions(data)
    cut_edges = _find_cut_edges(data)

    ρP = Ref(config.rho_flow_p)
    ρQ = Ref(config.rho_flow_q)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual

    dualsP =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))
    dualsQ =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))
    zP = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zQ = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zP_prev = Dict{Int, Float64}()
    zQ_prev = Dict{Int, Float64}()

    res_pri_hist = Float64[]
    res_du_hist = Float64[]
    eps_pri_hist = Float64[]
    eps_du_hist = Float64[]
    rhoP_hist = Float64[]
    rhoQ_hist = Float64[]

    optimizer_factory = resolve_optimizer_ac(config)

    for (zone_id, sub) in partitions
        m = build_model(sub, AC; optimizer = optimizer_factory, build_only = true)
        sub["model"] = m
        sub["obj_base"] = m[:generation_cost]
        sub["theta_var"] = m[:theta]
        sub["vmag_var"] = m[:v_magnitude]
        sub["zone_edges"] = _zone_cut_edges(sub, cut_edges)
        sub["has_solution"] = false
        _set_model_silent!(m)

        @inbounds for I in keys(sub["theta_var"])
            set_start_value(sub["theta_var"][I], 0.0)
        end
        @inbounds for I in keys(sub["vmag_var"])
            set_start_value(sub["vmag_var"][I], 1.0)
        end

        # pré-quad branch
        _build_quad_branch_expr_ac!(sub, ρP[], ρQ[], sub["zone_edges"])
        @objective(m, Min, sub["obj_base"] + sub["quad_branch_expr_ac"])
    end

    prev_ρP, prev_ρQ = ρP[], ρQ[]
    zone_ids = collect(keys(partitions))

    for k in 1:config.max_iter
        _FAST_USE_LOGS && println("\n[ADMM-AC-branch] Iteração $k")

        if ρP[] != prev_ρP || ρQ[] != prev_ρQ
            for (_, sub) in partitions
                _build_quad_branch_expr_ac!(sub, ρP[], ρQ[], sub["zone_edges"])
            end
            prev_ρP, prev_ρQ = ρP[], ρQ[]
        end

        nz = length(zone_ids)
        results = Vector{Any}(undef, nz)
        Threads.@threads for i in 1:nz
            z = zone_ids[i]
            sub = partitions[z]
            lin = _build_linear_branch_expr_ac!(
                sub, sub["zone_edges"], zP, zQ, dualsP[z], dualsQ[z], ρP[], ρQ[],
            )
            @objective(
                sub["model"],
                Min,
                sub["obj_base"] + sub["quad_branch_expr_ac"] + lin
            )
            optimize!(sub["model"])
            m = sub["model"]
            results[i] = (expr_cost = value(m[:generation_cost]),
                real_cost = compute_real_generation_cost(m, sub["objective_params"]))
        end

        for (i, z) in enumerate(zone_ids)
            sub = partitions[z]
            out = results[i]
            sub["has_solution"] = true
            sub["generation_cost_model_expr"] = out.expr_cost
            sub["real_generation_cost"] = out.real_cost
        end

        zP_prev = copy(zP)
        zQ_prev = copy(zQ)
        compute_ac_branch_consensus!(zP, zQ, partitions, cut_edges)
        update_duals_branch_ac!(dualsP, dualsQ, partitions, zP, zQ, cut_edges, ρP[], ρQ[])

        r = primal_residual_branch_ac(partitions, zP, zQ, cut_edges)
        s = dual_residual_branch_ac(zP, zQ, zP_prev, zQ_prev, cut_edges)
        εp = eps_primal_branch_ac(partitions, zP, zQ, cut_edges, ϵ_abs, ϵ_rel)
        εd = eps_dual_branch_ac(zP, zQ, zP_prev, zQ_prev, cut_edges, ϵ_abs, ϵ_rel)

        push!(res_pri_hist, r)
        push!(res_du_hist, s)
        push!(eps_pri_hist, εp)
        push!(eps_du_hist, εd)
        _FAST_USE_LOGS && println("‖r‖ = $r | ε_pri = $εp")
        _FAST_USE_LOGS && println("‖s‖ = $s | ε_du  = $εd")

        push!(rhoP_hist, ρP[])
        push!(rhoQ_hist, ρQ[])

        if r < εp && s < εd
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        adaptive_rho_ac!(ρP, r, s)
        adaptive_rho_ac!(ρQ, r, s)
    end

    save_ac_branch_plots(
        res_pri_hist,
        res_du_hist,
        eps_pri_hist,
        eps_du_hist,
        rhoP_hist,
        rhoQ_hist;
        stem = "admm_ac_branch",
    )

    total_real_cost = sum(p["real_generation_cost"] for p in values(partitions))
    total_expr_cost = sum(p["generation_cost_model_expr"] for p in values(partitions))
    println("\n[ADMM-AC-branch] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("ρP_final = $(ρP[]) | ρQ_final = $(ρQ[])")
    elapsed_s = (time_ns() - t0) / 1e9
    println("Tempo ADMM (branch/AC): $(_pretty_time_ac(elapsed_s))")

    return total_real_cost, partitions
end

# ========================== ADMM: HÍBRIDO (AC) ==========================
function run_admm_ac_hybrid(config::DistributedOPFModel, data::Dict{String, Any})
    t0 = time_ns()
    @assert config.base_model == AC "Hybrid AC: base_model deve ser AC."
    partitions, shared_buses = prepare_admm_partitions(data)
    cut_edges = _find_cut_edges(data)

    ρθ = Ref(config.rho_theta)
    ρV = Ref(config.rho_v)
    ρP = Ref(config.rho_flow_p)
    ρQ = Ref(config.rho_flow_q)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual

    shared_buses_sorted, b2pos = prepare_dense_maps_ac(shared_buses)

    dualsθ_dense = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )
    dualsV_dense = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )
    dualsP =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))
    dualsQ =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))

    θ̄_vec = zeros(Float64, length(shared_buses_sorted))
    V̄_vec = ones(Float64, length(shared_buses_sorted))
    θ̄_prev = similar(θ̄_vec)
    V̄_prev = similar(V̄_vec)

    zP = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zQ = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zP_prev = Dict{Int, Float64}()
    zQ_prev = Dict{Int, Float64}()

    res_pri_hist = Float64[]
    res_du_hist = Float64[]
    eps_pri_hist = Float64[]
    eps_du_hist = Float64[]
    rhoθ_hist = Float64[]
    rhoV_hist = Float64[]
    rhoP_hist = Float64[]
    rhoQ_hist = Float64[]
    θ_avg_hist = Dict{Int, Vector{Float64}}()

    optimizer_factory = resolve_optimizer_ac(config)

    # build inicial + hotpaths + termos quadráticos
    for (zone_id, sub) in partitions
        m = build_model(sub, AC; optimizer = optimizer_factory, build_only = true)
        sub["model"] = m
        sub["obj_base"] = m[:generation_cost]
        sub["theta_var"] = m[:theta]
        sub["vmag_var"] = m[:v_magnitude]
        sub["zone_edges"] = _zone_cut_edges(sub, cut_edges)
        sub["has_solution"] = false
        _set_model_silent!(m)

        @inbounds for I in keys(sub["theta_var"])
            set_start_value(sub["theta_var"][I], 0.0)
        end
        @inbounds for I in keys(sub["vmag_var"])
            set_start_value(sub["vmag_var"][I], 1.0)
        end

        prepare_hotpaths_ac!(sub, shared_buses)
        _build_quad_bus_expr_ac!(sub, ρθ[], ρV[])
        _build_quad_branch_expr_ac!(sub, ρP[], ρQ[], sub["zone_edges"])
        @objective(
            m,
            Min,
            sub["obj_base"] + sub["quad_bus_expr_ac"] + sub["quad_branch_expr_ac"]
        )
    end

    prev_ρθ, prev_ρV = ρθ[], ρV[]
    prev_ρP, prev_ρQ = ρP[], ρQ[]
    zone_ids = collect(keys(partitions))
    λθ_zone_buf = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in zone_ids
    )
    λV_zone_buf = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in zone_ids
    )

    for k in 1:config.max_iter
        _FAST_USE_LOGS && println("\n[ADMM-AC-hybrid] Iteração $k")

        # rebuild apenas se algum ρ mudar
        if ρθ[] != prev_ρθ || ρV[] != prev_ρV
            for (_, sub) in partitions
                _build_quad_bus_expr_ac!(sub, ρθ[], ρV[])
            end
            prev_ρθ, prev_ρV = ρθ[], ρV[]
        end
        if ρP[] != prev_ρP || ρQ[] != prev_ρQ
            for (_, sub) in partitions
                _build_quad_branch_expr_ac!(sub, ρP[], ρQ[], sub["zone_edges"])
            end
            prev_ρP, prev_ρQ = ρP[], ρQ[]
        end

        nz = length(zone_ids)
        results = Vector{Any}(undef, nz)
        Threads.@threads for i in 1:nz
            z = zone_ids[i]
            sub = partitions[z]
            λθ = (λθ_zone_buf[z] .= dualsθ_dense[z])
            λV = (λV_zone_buf[z] .= dualsV_dense[z])

            lin_bus = _build_linear_bus_expr_ac!(
                sub, shared_buses_sorted, b2pos, θ̄_vec, V̄_vec, λθ, λV, ρθ[], ρV[],
            )
            lin_br = _build_linear_branch_expr_ac!(
                sub, sub["zone_edges"], zP, zQ, dualsP[z], dualsQ[z], ρP[], ρQ[],
            )

            @objective(sub["model"], Min,
                sub["obj_base"] + sub["quad_bus_expr_ac"] + sub["quad_branch_expr_ac"] +
                (lin_bus + lin_br))
            optimize!(sub["model"])
            m = sub["model"]
            results[i] = (expr_cost = value(m[:generation_cost]),
                real_cost = compute_real_generation_cost(m, sub["objective_params"]))
        end

        for (i, z) in enumerate(zone_ids)
            sub = partitions[z]
            out = results[i]
            sub["has_solution"] = true
            sub["generation_cost_model_expr"] = out.expr_cost
            sub["real_generation_cost"] = out.real_cost
        end

        θ̄_prev .= θ̄_vec
        V̄_prev .= V̄_vec
        zP_prev = copy(zP)
        zQ_prev = copy(zQ)

        compute_ac_bus_consensus_vec!(
            θ̄_vec,
            V̄_vec,
            partitions,
            shared_buses_sorted,
            b2pos,
        )
        @inbounds for b in shared_buses_sorted
            j = b2pos[b]
            push!(get!(θ_avg_hist, b, Float64[]), θ̄_vec[j])
        end
        compute_ac_branch_consensus!(zP, zQ, partitions, cut_edges)

        for z in zone_ids
            update_duals_bus_ac_vec!(dualsθ_dense[z], dualsV_dense[z], partitions[z],
                θ̄_vec, V̄_vec, ρθ[], ρV[], shared_buses_sorted, b2pos)
        end
        update_duals_branch_ac!(dualsP, dualsQ, partitions, zP, zQ, cut_edges, ρP[], ρQ[])

        r_bus = primal_residual_bus_ac_vec(
            partitions,
            θ̄_vec,
            V̄_vec,
            shared_buses_sorted,
            b2pos,
        )
        r_br = primal_residual_branch_ac(partitions, zP, zQ, cut_edges)
        r_tot = hypot(r_bus, r_br)

        s_bus = dual_residual_bus_ac_vec(θ̄_vec, V̄_vec, θ̄_prev, V̄_prev)
        s_br = dual_residual_branch_ac(zP, zQ, zP_prev, zQ_prev, cut_edges)
        s_tot = hypot(s_bus, s_br)

        εb_pri = eps_primal_bus_ac_vec(
            partitions,
            θ̄_vec,
            V̄_vec,
            shared_buses_sorted,
            b2pos,
            ϵ_abs,
            ϵ_rel,
        )
        εe_pri = eps_primal_branch_ac(partitions, zP, zQ, cut_edges, ϵ_abs, ϵ_rel)
        ε_pri = hypot(εb_pri, εe_pri)

        εb_du = eps_dual_bus_ac_vec(θ̄_vec, V̄_vec, θ̄_prev, V̄_prev, ϵ_abs, ϵ_rel)
        εe_du = eps_dual_branch_ac(zP, zQ, zP_prev, zQ_prev, cut_edges, ϵ_abs, ϵ_rel)
        ε_du = hypot(εb_du, εe_du)

        push!(res_pri_hist, r_tot)
        push!(res_du_hist, s_tot)
        push!(eps_pri_hist, ε_pri)
        push!(eps_du_hist, ε_du)
        _FAST_USE_LOGS && println("‖r‖ = $r_tot | ε_pri = $ε_pri")
        _FAST_USE_LOGS && println("‖s‖ = $s_tot | ε_du  = $ε_du")

        push!(rhoθ_hist, ρθ[])
        push!(rhoV_hist, ρV[])
        push!(rhoP_hist, ρP[])
        push!(rhoQ_hist, ρQ[])

        if r_tot < ε_pri && s_tot < ε_du
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        adaptive_rho_ac!(ρθ, r_bus, s_bus)
        adaptive_rho_ac!(ρV, r_bus, s_bus)
        adaptive_rho_ac!(ρP, r_br, s_br)
        adaptive_rho_ac!(ρQ, r_br, s_br)
    end

    plot_admm_convergence(
        res_pri_hist,
        res_du_hist,
        eps_pri_hist,
        eps_du_hist,
        "plot_convergencia_admm_ac_hybrid.html",
    )
    PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(;
                    x = 1:length(rhoθ_hist),
                    y = rhoθ_hist,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρθ",
                );
                p.layout.title = "Evolução de ρθ";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "ρθ";
                p.layout.yaxis.type = "log";
                p
            )
        )(), "plot_rho_theta_ac_hybrid.html")
    PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(;
                    x = 1:length(rhoV_hist),
                    y = rhoV_hist,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρV",
                );
                p.layout.title = "Evolução de ρV";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "ρV";
                p.layout.yaxis.type = "log";
                p
            )
        )(), "plot_rho_v_ac_hybrid.html")
    PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(;
                    x = 1:length(rhoP_hist),
                    y = rhoP_hist,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρP",
                );
                p.layout.title = "Evolução de ρP";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "ρP";
                p.layout.yaxis.type = "log";
                p
            )
        )(), "plot_rho_p_ac_hybrid.html")
    PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(;
                    x = 1:length(rhoQ_hist),
                    y = rhoQ_hist,
                    type = :scatter,
                    mode = "lines+markers",
                    name = "ρQ",
                );
                p.layout.title = "Evolução de ρQ";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "ρQ";
                p.layout.yaxis.type = "log";
                p
            )
        )(), "plot_rho_q_ac_hybrid.html")

    plot_theta_avg_heatmap(θ_avg_hist, "plot_theta_avg_heatmap_ac_hybrid.html")

    total_real_cost = sum(p["real_generation_cost"] for p in values(partitions))
    total_expr_cost = sum(p["generation_cost_model_expr"] for p in values(partitions))
    println("\n[ADMM-AC-hybrid] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println(
        "ρθ_final = $(ρθ[]) | ρV_final = $(ρV[]) | ρP_final = $(ρP[]) | ρQ_final = $(ρQ[])",
    )
    elapsed_s = (time_ns() - t0) / 1e9
    println("Tempo ADMM (hybrid/AC): $(_pretty_time_ac(elapsed_s))")

    return total_real_cost, partitions
end
