# ================= SOCP Utils (revisado para espelhar DC/AC) =================

# ---- Helper: warm-start on/off (espelha AC/DC) ----
# Se já existir no seu core (ex.: Core/args.jl), remova esta definição.
_warm_on(config) =
    (hasproperty(config, :warm_start) && getproperty(config, :warm_start) == true)

"Retorna true se o modelo `m` está usando Gurobi."
is_gurobi(m::JuMP.Model)::Bool = occursin("Gurobi", MOI.get(m, MOI.SolverName()))

"Ridge suave para estabilizar QCP no runner branch (P,Q)."
function _soft_anchor_branch_expr!(sub::Dict{String, Any};
    κV2::Float64 = 1e-6, κPQ::Float64 = 1e-8,
    use_bus::Bool = true, use_pq::Bool = true,
)
    m = sub["model"]
    q = zero(JuMP.QuadExpr)

    if use_bus && haskey(m, :v_magnitude_sq)
        V2 = m[:v_magnitude_sq]
        sub["v2_var"] = V2
        # ancora TODAS as barras da zona (não só compartilhadas)
        for I in keys(V2)
            x = V2[I]
            # start_value se houver, senão 1.0 (flat seguro p/ SOCP)
            xref = (JuMP.start_value(x) === nothing) ? 1.0 : JuMP.start_value(x)
            q += 0.5 * κV2 * (x - xref)^2
        end
    end

    if use_pq && haskey(sub, "zone_edges")
        for e in sub["zone_edges"]
            l = get(sub["edge_to_local"], e, nothing)
            l === nothing && continue
            P, Q = _socp_flow_from_exprs_existing(m, l)   # deve ser AffExpr!
            q += 0.5 * κPQ * (P^2 + Q^2)
        end
    end

    sub["soft_anchor_expr"] = q
    return q
end

# ---- Otimizador: tenta Gurobi -> COSMO -> SCS (com override) ----
function resolve_socp_optimizer(config)
    # 1) prioridade absoluta: se vier um otimizador explícito no config
    if hasproperty(config, :optimizer) && getproperty(config, :optimizer) !== nothing
        return getproperty(config, :optimizer)
    end

    # 2) ordem preferida (opcional): :auto | :gurobi | :cosmo | :scs
    prefer =
        hasproperty(config, :prefer_socp_solver) ?
        getproperty(config, :prefer_socp_solver) : :auto

    # helper para montar a cadeia na ordem desejada
    function _chain(prefer_sym)
        if prefer_sym == :gurobi
            return (:gurobi, :cosmo, :scs)
        elseif prefer_sym == :cosmo
            return (:cosmo, :gurobi, :scs)
        elseif prefer_sym == :scs
            return (:scs, :gurobi, :cosmo)
        else
            return (:gurobi, :cosmo, :scs)  # :auto
        end
    end
    order = _chain(prefer)

    # 3) factory que instancia o primeiro solver que funcionar
    return () -> begin
        for s in order
            try
                if s === :gurobi
                    return Gurobi.Optimizer()
                elseif s === :cosmo
                    return COSMO.Optimizer()
                else
                    return SCS.Optimizer()
                end
            catch
                # falhou criar esse solver? tenta o próximo da cadeia
            end
        end
        # fallback paranoico (improvável chegar aqui)
        return SCS.Optimizer()
    end
end

# (mantido, mas sem forçar nada além do que você já usa)
safe_set_optimizer_attr!(m, key, val) = (
    try
        set_optimizer_attribute(m, key, val)
    catch
        nothing
    end
)

"Configura tolerâncias/accel do SCS (mantido; use se quiser)."
function set_socp_defaults!(m; eps_abs = 1e-6, eps_rel = 1e-6, max_iters = 200_000,
    acceleration_lookback = 10, acceleration_interval = 10,
    alpha = 1.8, normalize = true, scale = 1.0)
    safe_set_optimizer_attr!(m, "eps_abs", eps_abs)
    safe_set_optimizer_attr!(m, "eps_rel", eps_rel)
    safe_set_optimizer_attr!(m, "max_iters", max_iters)
    safe_set_optimizer_attr!(m, "acceleration_lookback", acceleration_lookback)
    safe_set_optimizer_attr!(m, "acceleration_interval", acceleration_interval)
    safe_set_optimizer_attr!(m, "alpha", alpha)
    safe_set_optimizer_attr!(m, "normalize", normalize)
    safe_set_optimizer_attr!(m, "scale", scale)
    return nothing
end

# ----------------- Helpers (paridade com DC/AC) -----------------
mutable struct FISTAStateSOCP
    β::Float64
    α::Float64
    prev_metric::Float64
end
FISTAStateSOCP() = FISTAStateSOCP(1.0, 0.0, Inf)

function fista_update_socp!(st::FISTAStateSOCP; metric::Float64)
    if metric <= st.prev_metric
        βnext = (1.0 + sqrt(1.0 + 4.0 * st.β^2)) / 2.0
        st.α = (st.β - 1.0) / βnext
        st.β = βnext
        st.prev_metric = metric
        return false
    else
        st.β = 1.0
        st.α = 0.0
        st.prev_metric = metric
        return true
    end
end

@inline function _extrapolate_socp!(y_hat::AbstractVector{T}, y::AbstractVector{T},
    y_prev::AbstractVector{T}, γ::T) where {T <: Real}
    @inbounds @simd for i in eachindex(y)
        y_hat[i] = y[i] + γ * (y[i] - y_prev[i])
    end
    return nothing
end

# ----------------- Hot-path / buffers (bus) -----------------
function prepare_socp_hotpaths!(sub::Dict{String, Any}, shared_buses_sorted::Vector{Int})
    # mapas de índice local para acesso rápido
    bus_ids = sub["bus"]["bus_i"]
    bus_local_index = Dict{Int, Int}()
    @inbounds for (i, b) in enumerate(bus_ids)
        bus_local_index[b] = i
    end
    sub["bus_local_index"] = bus_local_index

    # filtra shared da zona
    local_shared = get(sub, "shared_buses_local", Int[])
    if isempty(local_shared)
        busset = _bus_set(sub)
        sub["shared_buses_local"] = [b for b in shared_buses_sorted if b in busset]
        local_shared = sub["shared_buses_local"]
    end

    # posições locais das compartilhadas
    shared_pos_local = Vector{Int}(undef, length(local_shared))
    @inbounds for (k, b) in enumerate(local_shared)
        shared_pos_local[k] = bus_local_index[b]
    end
    sub["shared_pos_local"] = shared_pos_local

    # buffers
    sub["v2_buf"] = zeros(Float64, length(bus_ids))
    sub["v2_shared_buf"] = zeros(Float64, length(shared_pos_local))
    return nothing
end

@inline function snapshot_v2!(sub::Dict{String, Any})
    m = sub["model"]
    V2 = m[:v_magnitude_sq]
    buf = sub["v2_buf"]
    ids = sub["bus"]["bus_i"]
    @inbounds for i in eachindex(ids)
        buf[i] =
            JuMP.has_values(m) ? value(V2[ids[i], 1]) :
            (
                JuMP.start_value(V2[ids[i], 1]) === nothing ? 1.0 :
                JuMP.start_value(V2[ids[i], 1])
            )
    end
    return buf
end

@inline function snapshot_v2_shared!(sub::Dict{String, Any})
    m = sub["model"]
    V2 = m[:v_magnitude_sq]
    buf = sub["v2_shared_buf"]
    pos = sub["shared_pos_local"]
    ids = sub["bus"]["bus_i"]
    @inbounds for i in eachindex(pos)
        bidx = pos[i]
        buf[i] =
            JuMP.has_values(m) ? value(V2[ids[bidx], 1]) :
            (
                JuMP.start_value(V2[ids[bidx], 1]) === nothing ? 1.0 :
                JuMP.start_value(V2[ids[bidx], 1])
            )
    end
    return buf
end

# ----------------- Construção de termos (cache de quadrático) -----------------
function _build_quad_bus_expr_socp!(sub::Dict{String, Any}, ρ::Float64)
    m = sub["model"]
    V2 = m[:v_magnitude_sq]
    q = zero(JuMP.QuadExpr)
    # penaliza APENAS barras compartilhadas (padrão ADMM/DC/AC)
    for b in sub["shared_buses_local"]
        q += 0.5 * ρ * V2[b, 1]^2
    end
    sub["quad_bus_expr"] = q
    sub["ρ_bus_cached"] = ρ
    return sub["quad_bus_expr"]
end

function _build_quad_branch_rt_expr!(sub::Dict{String, Any}, ρR::Float64, ρT::Float64,
    zone_edges::Vector{Int})
    m = sub["model"]
    q = zero(JuMP.QuadExpr)
    @inbounds for e in zone_edges
        l_local = get(sub["edge_to_local"], e, nothing)
        l_local === nothing && continue
        Rloc, Tloc = _rt_from_exprs_socp(m, l_local)
        q += 0.5 * ρR * Rloc^2 + 0.5 * ρT * Tloc^2
    end
    sub["quad_rt_expr"] = q
    sub["ρR_cached"] = ρR
    sub["ρT_cached"] = ρT
    return q
end

function _build_quad_branch_pq_expr!(sub::Dict{String, Any}, ρP::Float64, ρQ::Float64,
    zone_edges::Vector{Int})
    m = sub["model"]
    q = zero(JuMP.QuadExpr)
    @inbounds for e in zone_edges
        l_local = get(sub["edge_to_local"], e, nothing)
        l_local === nothing && continue
        P_loc, Q_loc = _socp_flow_from_exprs_existing(m, l_local)
        q += 0.5 * ρP * P_loc^2 + 0.5 * ρQ * Q_loc^2
    end
    sub["quad_pq_expr"] = q
    sub["ρP_cached"] = ρP
    sub["ρQ_cached"] = ρQ
    return q
end

# ----- termos lineares (dependem de λ e consenso) -----
function _build_linear_bus_expr_socp(sub::Dict{String, Any},
    shared_buses_sorted::Vector{Int}, b2pos::Dict{Int, Int},
    v2bar_vec::Vector{Float64}, λ_vec_zone::Vector{Float64}, ρV2::Float64)
    m = sub["model"]
    V2 = m[:v_magnitude_sq]
    lin = zero(JuMP.AffExpr)
    for b in sub["shared_buses_local"]
        j = b2pos[b]
        lin += (λ_vec_zone[j] - ρV2 * v2bar_vec[j]) * V2[b, 1]
    end
    return lin
end

function _build_linear_branch_rt_expr(sub::Dict{String, Any},
    zone_edges::Vector{Int}, zR::Dict{Int, Float64}, zT::Dict{Int, Float64},
    λR_zone::Dict{Int, Float64}, λT_zone::Dict{Int, Float64}, ρR::Float64, ρT::Float64)
    m = sub["model"]
    lin = zero(JuMP.AffExpr)
    @inbounds for e in zone_edges
        l_local = get(sub["edge_to_local"], e, nothing)
        l_local === nothing && continue
        Rloc, Tloc = _rt_from_exprs_socp(m, l_local)
        lin += (get(λR_zone, e, 0.0) - ρR * get(zR, e, 0.0)) * Rloc
        lin += (get(λT_zone, e, 0.0) - ρT * get(zT, e, 0.0)) * Tloc
    end
    return lin
end

function _build_linear_branch_pq_expr(sub::Dict{String, Any},
    zone_edges::Vector{Int}, zP::Dict{Int, Float64}, zQ::Dict{Int, Float64},
    λP_zone::Dict{Int, Float64}, λQ_zone::Dict{Int, Float64}, ρP::Float64, ρQ::Float64)
    m = sub["model"]
    lin = zero(JuMP.AffExpr)
    @inbounds for e in zone_edges
        l_local = get(sub["edge_to_local"], e, nothing)
        l_local === nothing && continue
        P_loc, Q_loc = _socp_flow_from_exprs_existing(m, l_local)
        lin += (get(λP_zone, e, 0.0) - ρP * get(zP, e, 0.0)) * P_loc
        lin += (get(λQ_zone, e, 0.0) - ρQ * get(zQ, e, 0.0)) * Q_loc
    end
    return lin
end

# ----------------- Consensos (forma densa p/ BUS) -----------------
function prepare_dense_maps_socp(shared_buses::Vector{Int})
    shared_buses_sorted = sort(shared_buses)
    b2pos = Dict{Int, Int}(b => i for (i, b) in enumerate(shared_buses_sorted))
    return shared_buses_sorted, b2pos
end

function compute_v2_consensus_vec!(v2bar_vec::Vector{Float64},
    partitions::Dict{Int, Dict{String, Any}},
    shared_buses_sorted::Vector{Int}, b2pos::Dict{Int, Int})
    fill!(v2bar_vec, 0.0)
    counts = zeros(Int, length(shared_buses_sorted))
    for (_, sub) in partitions
        v2s = snapshot_v2_shared!(sub)
        lbs = sub["shared_buses_local"]
        @inbounds for i in eachindex(lbs)
            j = b2pos[lbs[i]]
            v2bar_vec[j] += v2s[i]
            counts[j] += 1
        end
    end
    @inbounds for j in eachindex(v2bar_vec)
        counts[j] != 0 && (v2bar_vec[j] /= counts[j])
    end
end

# versões baseadas em Dict (mantidas para edges)
function compute_rt_consensus!(zR::Dict{Int, Float64}, zT::Dict{Int, Float64},
    partitions, cut_edges::Vector{Int})
    sumsR = Dict{Int, Float64}()
    sumsT = Dict{Int, Float64}()
    cnts = Dict{Int, Int}()
    for (_, sub) in partitions
        @inbounds for e in cut_edges
            l = get(sub["edge_to_local"], e, nothing)
            l === nothing && continue
            R, T = _rt_from_exprs_socp(sub["model"], l)
            sumsR[e] = get(sumsR, e, 0.0) + value(R)
            sumsT[e] = get(sumsT, e, 0.0) + value(T)
            cnts[e] = get(cnts, e, 0) + 1
        end
    end
    @inbounds for e in cut_edges
        n = get(cnts, e, 0)
        n > 0 && (zR[e] = sumsR[e] / n; zT[e] = sumsT[e] / n)
    end
end

function compute_pq_consensus!(zP::Dict{Int, Float64}, zQ::Dict{Int, Float64},
    partitions, cut_edges::Vector{Int})
    sumsP = Dict{Int, Float64}()
    sumsQ = Dict{Int, Float64}()
    cnts = Dict{Int, Int}()
    for (_, sub) in partitions
        @inbounds for e in cut_edges
            l = get(sub["edge_to_local"], e, nothing)
            l === nothing && continue
            P, Q = _socp_flow_from_exprs_existing(sub["model"], l)
            sumsP[e] = get(sumsP, e, 0.0) + value(P)
            sumsQ[e] = get(sumsQ, e, 0.0) + value(Q)
            cnts[e] = get(cnts, e, 0) + 1
        end
    end
    @inbounds for e in cut_edges
        n = get(cnts, e, 0)
        n > 0 && (zP[e] = sumsP[e] / n; zQ[e] = sumsQ[e] / n)
    end
end

# ----------------- Duais (forma vetorial p/ BUS) -----------------
@inline function update_duals_v2_vec!(dual_vec_zone::Vector{Float64},
    sub::Dict{String, Any},
    v2bar_vec::Vector{Float64}, v2bar_prev_vec::Vector{Float64},
    ρ::Float64, α::Union{Nothing, Float64},
    shared_buses_sorted::Vector{Int}, b2pos::Dict{Int, Int})
    v2s = snapshot_v2_shared!(sub)
    lbs = sub["shared_buses_local"]
    if isnothing(α)
        @inbounds for i in eachindex(lbs)
            j = b2pos[lbs[i]]
            dual_vec_zone[j] += ρ * (v2s[i] - v2bar_vec[j])
        end
    else
        αv = something(α, 1.0)
        @inbounds for i in eachindex(lbs)
            j = b2pos[lbs[i]]
            v_rel = αv * v2s[i] + (1 - αv) * v2bar_prev_vec[j]
            dual_vec_zone[j] += ρ * (v_rel - v2bar_vec[j])
        end
    end
end

# ----------------- Resíduos / tolerâncias (bus em forma densa) -----------------
function primal_residual_bus_socp_vec(partitions, v2bar_vec, shared_buses_sorted, b2pos)
    acc = 0.0
    for (_, sub) in partitions
        v2buf = snapshot_v2!(sub)
        pos = sub["bus_local_index"]
        @inbounds for b in sub["shared_buses_local"]
            j = b2pos[b]
            acc += (v2buf[pos[b]] - v2bar_vec[j])^2
        end
    end
    return sqrt(acc)
end

dual_residual_bus_socp_vec(v2bar_vec, v2bar_prev_vec) =
    sqrt(sum((v2bar_vec[j] - v2bar_prev_vec[j])^2 for j in eachindex(v2bar_vec)))

function eps_primal_bus_socp_vec(
    partitions,
    v2bar_vec,
    shared_buses_sorted,
    b2pos,
    ϵ_abs,
    ϵ_rel,
)
    N = max(1, length(shared_buses_sorted))
    norm_x_sq = 0.0
    for (_, sub) in partitions
        v2buf = snapshot_v2!(sub)
        pos = sub["bus_local_index"]
        @inbounds for b in sub["shared_buses_local"]
            norm_x_sq += (v2buf[pos[b]])^2
        end
    end
    norm_z_sq = sum((v2bar_vec[j])^2 for j in eachindex(v2bar_vec))
    return sqrt(N) * ϵ_abs + ϵ_rel * max(sqrt(norm_x_sq), sqrt(norm_z_sq))
end

function eps_dual_bus_socp_vec(v2bar_vec, v2bar_prev_vec, ϵ_abs, ϵ_rel)
    N = max(1, length(v2bar_vec))
    return sqrt(N) * ϵ_abs +
           ϵ_rel * sqrt(sum((v2bar_vec[j] - v2bar_prev_vec[j])^2
                    for j in eachindex(v2bar_vec)))
end

# ----------------- Resíduos / tolerâncias (edges) -----------------
function primal_residual_branch_rt_socp(partitions, zR, zT, cut_edges)
    acc = 0.0
    for (_, sub) in partitions
        @inbounds for e in cut_edges
            l = get(sub["edge_to_local"], e, nothing)
            l === nothing && continue
            R, T = _rt_from_exprs_socp(sub["model"], l)
            acc += (value(R) - zR[e])^2 + (value(T) - zT[e])^2
        end
    end
    return sqrt(acc)
end

dual_residual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges) =
    sqrt(
        sum(
            (get(zR, e, 0.0) - get(zR_prev, e, 0.0))^2 +
            (get(zT, e, 0.0) - get(zT_prev, e, 0.0))^2
            for e in cut_edges
        ),
    )

function eps_primal_branch_rt_socp(partitions, zR, zT, cut_edges, ϵ_abs, ϵ_rel)
    dim = 0
    for (_, sub) in partitions
        dim += 2 * length(_zone_cut_edges(sub, cut_edges))
    end
    z_norm_sq = 0.0
    @inbounds for e in cut_edges
        z_norm_sq += get(zR, e, 0.0)^2 + get(zT, e, 0.0)^2
    end
    return sqrt(max(1, dim)) * ϵ_abs + ϵ_rel * max(sqrt(z_norm_sq), 1e-12)
end

function eps_dual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges, ϵ_abs, ϵ_rel)
    dim = 2 * max(1, length(cut_edges))
    s = dual_residual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges)
    return sqrt(dim) * ϵ_abs + ϵ_rel * s
end

function primal_residual_branch_pq_socp(partitions, zP, zQ, cut_edges)
    acc = 0.0
    for (_, sub) in partitions
        @inbounds for e in cut_edges
            l = get(sub["edge_to_local"], e, nothing)
            l === nothing && continue
            P, Q = _socp_flow_from_exprs_existing(sub["model"], l)
            acc += (value(P) - zP[e])^2 + (value(Q) - zQ[e])^2
        end
    end
    return sqrt(acc)
end

dual_residual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges) =
    sqrt(
        sum(
            (get(zP, e, 0.0) - get(zP_prev, e, 0.0))^2 +
            (get(zQ, e, 0.0) - get(zQ_prev, e, 0.0))^2
            for e in cut_edges
        ),
    )

function eps_primal_branch_pq_socp(partitions, zP, zQ, cut_edges, ϵ_abs, ϵ_rel)
    dim = 0
    for (_, sub) in partitions
        dim += 2 * length(_zone_cut_edges(sub, cut_edges))
    end
    z_norm_sq = 0.0
    @inbounds for e in cut_edges
        z_norm_sq += get(zP, e, 0.0)^2 + get(zQ, e, 0.0)^2
    end
    return sqrt(max(1, dim)) * ϵ_abs + ϵ_rel * max(sqrt(z_norm_sq), 1e-12)
end

function eps_dual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges, ϵ_abs, ϵ_rel)
    dim = 2 * max(1, length(cut_edges))
    s = dual_residual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges)
    return sqrt(dim) * ϵ_abs + ϵ_rel * s
end

# ---------- ρ adaptativo (paridade com DC) ----------
function adaptive_rho!(ρ::Base.RefValue{Float64}, r::Float64, s::Float64;
    μ_incr = 2.0, μ_decr = 2.0, τ_incr = 1.05, τ_decr = 0.98)
    if r > μ_incr * s
        ρ[] *= τ_incr
    elseif s > μ_decr * r
        ρ[] *= τ_decr
    end
end

# ---------- Plots (wrappers) ----------
function save_socp_v2rt_plots(r_hist, s_hist, eps_r_hist, eps_s_hist,
    rhoV2_hist, rhoR_hist, rhoT_hist; stem = "admm_socp_v2rt")
    plot_admm_convergence(r_hist, s_hist, eps_r_hist, eps_s_hist,
        "plot_convergencia_$(stem).html")
    PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(; x = 1:length(rhoV2_hist), y = rhoV2_hist, type = :scatter,
                    mode = "lines+markers", name = "ρ(V²)");
                p.layout.title = "Evolução de ρ(V²)";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "ρ(V²)";
                p.layout.yaxis.type = "log";
                p
            )
        )(), "plot_rho_v2_$(stem).html")
    PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(; x = 1:length(rhoR_hist), y = rhoR_hist, type = :scatter,
                    mode = "lines+markers", name = "ρR");
                p.layout.title = "Evolução de ρR";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "ρR";
                p.layout.yaxis.type = "log";
                p
            )
        )(), "plot_rhoR_$(stem).html")
    return PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(; x = 1:length(rhoT_hist), y = rhoT_hist, type = :scatter,
                    mode = "lines+markers", name = "ρT");
                p.layout.title = "Evolução de ρT";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "ρT";
                p.layout.yaxis.type = "log";
                p
            )
        )(), "plot_rhoT_$(stem).html")
end

function save_socp_branch_pq_plots(r_hist, s_hist, eps_r_hist, eps_s_hist,
    rhoP_hist, rhoQ_hist; stem = "admm_socp_branch")
    plot_admm_convergence(r_hist, s_hist, eps_r_hist, eps_s_hist,
        "plot_convergencia_$(stem).html")
    PlotlyLight.save(
        (
            () -> (
                p = Plot();
                p(; x = 1:length(rhoP_hist), y = rhoP_hist, type = :scatter,
                    mode = "lines+markers", name = "ρP");
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
                p(; x = 1:length(rhoQ_hist), y = rhoQ_hist, type = :scatter,
                    mode = "lines+markers", name = "ρQ");
                p.layout.title = "Evolução de ρQ";
                p.layout.xaxis.title = "Iterações";
                p.layout.yaxis.title = "ρQ";
                p.layout.yaxis.type = "log";
                p
            )
        )(), "plot_rhoQ_$(stem).html")
end

function save_socp_hybrid_plots(r_hist, s_hist, eps_r_hist, eps_s_hist,
    rhoV2_hist, rhoP_hist, rhoQ_hist, rhoR_hist, rhoT_hist;
    stem = "admm_socp_hybrid")
    plot_admm_convergence(r_hist, s_hist, eps_r_hist, eps_s_hist,
        "plot_convergencia_$(stem).html")
    for (lab, hist) in [("ρ(V²)", rhoV2_hist), ("ρP", rhoP_hist), ("ρQ", rhoQ_hist),
        ("ρR", rhoR_hist), ("ρT", rhoT_hist)]
        fname = replace(lab, "/" => "") |> x -> lowercase(replace(x, "ρ" => "rho"))
        PlotlyLight.save(
            (
                () -> (
                    p = Plot();
                    p(; x = 1:length(hist), y = hist, type = :scatter,
                        mode = "lines+markers", name = lab);
                    p.layout.title = "Evolução de $(lab)";
                    p.layout.xaxis.title = "Iterações";
                    p.layout.yaxis.title = lab;
                    p.layout.yaxis.type = "log";
                    p
                )
            )(), "plot_$(fname)_$(stem).html")
    end
end

# ================== RUNNERS (espelhando DC/AC) ==================

# ---- ADMM-SOCP: BUS (V²) + BRANCH (R,T) com estrutura DC-like ----
function run_admm_socp(config::DistributedOPFModel, data::Dict{String, Any})
    @assert config.base_model == SOCP "ADMM-SOCP: base_model deve ser SOCP."

    t0_ns = time_ns()
    iter_done = 0

    partitions, shared_buses = prepare_admm_partitions(data)
    cut_edges = _find_cut_edges(data)

    # estados de ρ
    ρV2 = Ref(config.rho_v)
    ρR = Ref(config.rho_flow_p)
    ρT = Ref(config.rho_flow_q)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual
    adaptive = _adaptive_on(config)
    accel = _accel_on(config)   # aceleração estilo DC: só no bloco BUS
    warm = _warm_on(config)

    # construtor de modelos
    optimizer_factory = resolve_socp_optimizer(config)

    # mapeamento denso (BUS)
    shared_buses_sorted, b2pos = prepare_dense_maps_socp(shared_buses)

    # build por zona
    for (z, sub) in partitions
        m = build_model(sub, SOCP; optimizer = optimizer_factory, build_only = true)
        sub["model"] = m
        sub["obj_base"] = m[:generation_cost]
        sub["v2_var"] = m[:v_magnitude_sq]
        _set_model_silent!(m)

        # filtra arestas da zona
        sub["zone_edges"] = _zone_cut_edges(sub, cut_edges)

        # shared locais + buffers
        busset = _bus_set(sub)
        sub["shared_buses_local"] = [b for b in shared_buses_sorted if b in busset]
        prepare_socp_hotpaths!(sub, shared_buses_sorted)

        # ---- warm-start de variáveis (apenas se warm==true) ----
        if warm
            for I in keys(sub["v2_var"])
                set_start_value(sub["v2_var"][I], 1.0)
            end
        end

        # cache quadráticos
        _build_quad_bus_expr_socp!(sub, ρV2[])
        _build_quad_branch_rt_expr!(sub, ρR[], ρT[], sub["zone_edges"])

        # objetivo "estático" = base + quads (lin será somado por iteração)
        sub["obj_static"] = sub["obj_base"] + sub["quad_bus_expr"] + sub["quad_rt_expr"]
        @objective(m, Min, sub["obj_static"])
        optimize!(m)
    end

    # vetores densos (BUS)
    v2bar_vec = zeros(Float64, length(shared_buses_sorted))
    v2bar_prev_vec = similar(v2bar_vec)
    v2bar_hat = similar(v2bar_vec)
    fill!(v2bar_hat, 0.0)
    fista = FISTAStateSOCP()

    # duais densos por zona (BUS)
    duals_v2_dense = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )
    λ_prev_real = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )
    λ_hat = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )
    λ_zone_buf = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )

    # consensos/duais (edges RT)
    zR = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zT = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zR_prev = Dict{Int, Float64}()
    zT_prev = Dict{Int, Float64}()
    dualsR =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))
    dualsT =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))

    # históricos
    r_hist = Float64[]
    s_hist = Float64[]
    eps_r_hist = Float64[]
    eps_s_hist = Float64[]
    rhoV2_hist = Float64[]
    rhoR_hist = Float64[]
    rhoT_hist = Float64[]

    # -------- warm-start ADMM (apenas se warm==true) --------
    if warm
        compute_v2_consensus_vec!(v2bar_vec, partitions, shared_buses_sorted, b2pos)
        v2bar_hat .= v2bar_vec
        for (z, sub) in partitions
            v2s = snapshot_v2_shared!(sub)
            lbs = sub["shared_buses_local"]
            @inbounds for i in eachindex(lbs)
                j = b2pos[lbs[i]]
                duals_v2_dense[z][j] = ρV2[] * (v2s[i] - v2bar_vec[j])
            end
            λ_prev_real[z] .= duals_v2_dense[z]
            λ_hat[z] .= duals_v2_dense[z]
        end
        compute_rt_consensus!(zR, zT, partitions, cut_edges)
        for (z, sub) in partitions
            m = sub["model"]
            for e in cut_edges
                l = get(sub["edge_to_local"], e, nothing)
                if l !== nothing
                    R, T = _rt_from_exprs_socp(m, l)
                    dualsR[z][e] = ρR[] * (value(R) - zR[e])
                    dualsT[z][e] = ρT[] * (value(T) - zT[e])
                end
            end
        end
    else
        # flat start seguro quando warm==false
        fill!(v2bar_vec, 1.0)
        v2bar_hat .= v2bar_vec
        # duais já estão zero; zR/zT já zero
        # λ_prev_real/λ_hat permanecem zero
    end

    zone_ids = collect(keys(partitions))
    prev_ρV2, prev_ρR, prev_ρT = ρV2[], ρR[], ρT[]

    for k in 1:config.max_iter
        _FAST_USE_LOGS && println("\n[ADMM-SOCP (V²+R,T)] Iteração $k")

        # se ρ mudou, só reconstroi os quadráticos (como no DC)
        if ρV2[] != prev_ρV2
            for (_, sub) in partitions
                _build_quad_bus_expr_socp!(sub, ρV2[])
                sub["obj_static"] =
                    sub["obj_base"] + sub["quad_bus_expr"] + sub["quad_rt_expr"]
            end
            prev_ρV2 = ρV2[]
        end
        if ρR[] != prev_ρR || ρT[] != prev_ρT
            for (_, sub) in partitions
                _build_quad_branch_rt_expr!(sub, ρR[], ρT[], sub["zone_edges"])
                sub["obj_static"] =
                    sub["obj_base"] + sub["quad_bus_expr"] + sub["quad_rt_expr"]
            end
            prev_ρR, prev_ρT = ρR[], ρT[]
        end

        α_used = accel ? nothing : _alpha_schedule(k, config.alpha_overrelax; k_on = 10)
        v2_target = accel ? v2bar_hat : v2bar_vec

        # resolve local (paralelizável)
        nz = length(zone_ids)
        results = Vector{Any}(undef, nz)
        @threads for i in 1:nz
            z = zone_ids[i]
            sub = partitions[z]
            m = sub["model"]

            # λ para V²: extrapolado quando aceleração ON
            λz = accel ? (λ_zone_buf[z] .= λ_hat[z]) : (λ_zone_buf[z] .= duals_v2_dense[z])

            lin_bus = _build_linear_bus_expr_socp(
                sub,
                shared_buses_sorted,
                b2pos,
                v2_target,
                λz,
                ρV2[],
            )
            lin_rt = _build_linear_branch_rt_expr(
                sub,
                sub["zone_edges"],
                zR,
                zT,
                dualsR[z],
                dualsT[z],
                ρR[],
                ρT[],
            )
            @objective(m, Min, sub["obj_static"] + (lin_bus + lin_rt))
            optimize!(m)

            gen_cost_expr = value(m[:generation_cost])
            real_cost = compute_real_generation_cost(m, sub["objective_params"])
            results[i] = (gen_cost_expr = gen_cost_expr, real_cost = real_cost)
        end

        for (i, z) in enumerate(zone_ids)
            sub = partitions[z]
            out = results[i]
            sub["generation_cost_model_expr"] = out.gen_cost_expr
            sub["real_generation_cost"] = out.real_cost
        end

        # ===== consenso + duais (reais) =====
        v2bar_prev_vec .= v2bar_vec
        zR_prev = copy(zR)
        zT_prev = copy(zT)

        compute_v2_consensus_vec!(v2bar_vec, partitions, shared_buses_sorted, b2pos)
        compute_rt_consensus!(zR, zT, partitions, cut_edges)

        # guardar λ reais antes do update (para FISTA)
        if accel
            for z in zone_ids
                λ_prev_real[z] .= duals_v2_dense[z]
            end
        end

        # duais reais
        for z in zone_ids
            update_duals_v2_vec!(duals_v2_dense[z], partitions[z], v2bar_vec,
                v2bar_prev_vec,
                ρV2[], α_used, shared_buses_sorted, b2pos)
        end
        for (z, sub) in partitions
            m = sub["model"]
            @inbounds for e in cut_edges
                l = get(sub["edge_to_local"], e, nothing)
                l === nothing && continue
                R, T = _rt_from_exprs_socp(m, l)
                dualsR[z][e] = get(dualsR[z], e, 0.0) + ρR[] * (value(R) - zR[e])
                dualsT[z][e] = get(dualsT[z], e, 0.0) + ρT[] * (value(T) - zT[e])
            end
        end

        # resíduos / tolerâncias
        r_bus =
            primal_residual_bus_socp_vec(partitions, v2bar_vec, shared_buses_sorted, b2pos)
        r_rt = primal_residual_branch_rt_socp(partitions, zR, zT, cut_edges)
        r = hypot(r_bus, r_rt)

        s_bus = dual_residual_bus_socp_vec(v2bar_vec, v2bar_prev_vec)
        s_rt = dual_residual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges)
        s = hypot(s_bus, s_rt)

        push!(r_hist, r)
        push!(s_hist, s)

        ε_pri = begin
            εb = eps_primal_bus_socp_vec(
                partitions,
                v2bar_vec,
                shared_buses_sorted,
                b2pos,
                ϵ_abs,
                ϵ_rel,
            )
            εrt = eps_primal_branch_rt_socp(partitions, zR, zT, cut_edges, ϵ_abs, ϵ_rel)
            hypot(εb, εrt)
        end
        ε_du = begin
            εb = eps_dual_bus_socp_vec(v2bar_vec, v2bar_prev_vec, ϵ_abs, ϵ_rel)
            εrt = eps_dual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges, ϵ_abs, ϵ_rel)
            hypot(εb, εrt)
        end
        push!(eps_r_hist, ε_pri)
        push!(eps_s_hist, ε_du)

        _FAST_USE_LOGS && println("‖r‖ = $r | ε_pri = $ε_pri")
        _FAST_USE_LOGS && println("‖s‖ = $s | ε_du  = $ε_du")

        push!(rhoV2_hist, ρV2[])
        push!(rhoR_hist, ρR[])
        push!(rhoT_hist, ρT[])

        iter_done = k
        if r < ε_pri && s < ε_du
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        # aceleração FISTA (só BUS), usando estados reais
        if accel
            metric_bus = _accel_metric(r_bus, s_bus)
            restarted = fista_update_socp!(fista; metric = metric_bus)
            if restarted
                v2bar_hat .= v2bar_vec
                for z in zone_ids
                    λ_hat[z] .= duals_v2_dense[z]
                end
            else
                _extrapolate_socp!(v2bar_hat, v2bar_vec, v2bar_prev_vec, fista.α)
                for z in zone_ids
                    _extrapolate_socp!(λ_hat[z], duals_v2_dense[z], λ_prev_real[z], fista.α)
                end
            end
        end

        # ρ adaptativo por bloco
        if adaptive
            adaptive_rho!(ρV2, r_bus, s_bus)
            adaptive_rho!(ρR, r_rt, s_rt)
            adaptive_rho!(ρT, r_rt, s_rt)
        end
    end

    save_socp_v2rt_plots(r_hist, s_hist, eps_r_hist, eps_s_hist,
        rhoV2_hist, rhoR_hist, rhoT_hist; stem = "admm_socp_no_theta")

    total_real_cost = sum(p["real_generation_cost"] for p in values(partitions))
    total_expr_cost = sum(p["generation_cost_model_expr"] for p in values(partitions))

    println("\n[ADMM-SOCP (V²+R,T)] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("ρ(V²)_final = $(ρV2[]) | ρR_final = $(ρR[]) | ρT_final = $(ρT[])")
    elapsed_s = (time_ns() - t0_ns) / 1e9
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (SOCP V²+RT): $(_pretty_time(elapsed_s))")

    return total_real_cost, partitions
end

# ---- ADMM-SOCP-branch (P,Q) ----
function run_admm_socp_branch(config::DistributedOPFModel, data::Dict{String, Any})
    @assert config.base_model == SOCP "ADMM-SOCP-branch: base_model deve ser SOCP."

    t0_ns = time_ns()
    iter_done = 0

    partitions, _ = prepare_admm_partitions(data)
    cut_edges = _find_cut_edges(data)

    ρP = Ref(config.rho_flow_p)
    ρQ = Ref(config.rho_flow_q)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual
    adaptive = _adaptive_on(config)
    warm = _warm_on(config)

    optimizer_factory = resolve_socp_optimizer(config)

    # build por zona
    for (z, sub) in partitions
        m = build_model(sub, SOCP; optimizer = optimizer_factory, build_only = true)
        sub["model"] = m
        sub["obj_base"] = m[:generation_cost]
        _set_model_silent!(m)
        sub["zone_edges"] = _zone_cut_edges(sub, cut_edges)

        # 🔧 ridge (apenas se Gurobi) — antes de definir obj_static
        if is_gurobi(m)
            _soft_anchor_branch_expr!(
                sub;
                κV2 = 1e-6,
                κPQ = 1e-8,
                use_bus = true,
                use_pq = true,
            )
        else
            sub["soft_anchor_expr"] = zero(JuMP.QuadExpr)
        end

        # quadrático do ADMM (P,Q)
        _build_quad_branch_pq_expr!(sub, ρP[], ρQ[], sub["zone_edges"])

        # 🔧 obj_static inclui SEMPRE o ridge (mesmo se zero)
        sub["obj_static"] = sub["obj_base"] + sub["quad_pq_expr"] + sub["soft_anchor_expr"]
        @objective(m, Min, sub["obj_static"])
        optimize!(m)
    end

    # consensos/duais
    zP = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zP_prev = Dict{Int, Float64}()
    zQ = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zQ_prev = Dict{Int, Float64}()
    dualsP =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))
    dualsQ =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))

    r_hist = Float64[]
    s_hist = Float64[]
    eps_r_hist = Float64[]
    eps_s_hist = Float64[]
    rhoP_hist = Float64[]
    rhoQ_hist = Float64[]

    prev_ρP, prev_ρQ = ρP[], ρQ[]

    # ---- warm-start ADMM (somente se warm==true) ----
    if warm
        compute_pq_consensus!(zP, zQ, partitions, cut_edges)
        for (z, sub) in partitions
            m = sub["model"]
            @inbounds for e in cut_edges
                l = get(sub["edge_to_local"], e, nothing)
                l === nothing && continue
                P, Q = _socp_flow_from_exprs_existing(m, l)
                dualsP[z][e] = ρP[] * (value(P) - zP[e])
                dualsQ[z][e] = ρQ[] * (value(Q) - zQ[e])
            end
        end
    else
        # zP, zQ já estão zero; duais permanecem zero (flat)
    end

    for k in 1:config.max_iter
        _FAST_USE_LOGS && println("\n[ADMM-SOCP-branch (P,Q)] Iteração $k")

        if ρP[] != prev_ρP || ρQ[] != prev_ρQ
            for (_, sub) in partitions
                _build_quad_branch_pq_expr!(sub, ρP[], ρQ[], sub["zone_edges"])
                sub["obj_static"] =
                    sub["obj_base"] + sub["quad_pq_expr"] + sub["soft_anchor_expr"]
            end
            prev_ρP, prev_ρQ = ρP[], ρQ[]
        end

        # locais (correto por zona)
        for (z, sub) in partitions
            m = sub["model"]
            lin = _build_linear_branch_pq_expr(sub, sub["zone_edges"], zP, zQ,
                dualsP[z], dualsQ[z], ρP[], ρQ[])
            @objective(m, Min, sub["obj_static"] + lin)
            optimize!(m)
            sub["generation_cost_model_expr"] = value(m[:generation_cost])
            sub["real_generation_cost"] =
                compute_real_generation_cost(m, sub["objective_params"])
        end

        zP_prev = copy(zP)
        zQ_prev = copy(zQ)
        compute_pq_consensus!(zP, zQ, partitions, cut_edges)

        for (z, sub) in partitions
            m = sub["model"]
            @inbounds for e in cut_edges
                l = get(sub["edge_to_local"], e, nothing)
                l === nothing && continue
                P, Q = _socp_flow_from_exprs_existing(m, l)
                dualsP[z][e] = get(dualsP[z], e, 0.0) + ρP[] * (value(P) - zP[e])
                dualsQ[z][e] = get(dualsQ[z], e, 0.0) + ρQ[] * (value(Q) - zQ[e])
            end
        end

        r = primal_residual_branch_pq_socp(partitions, zP, zQ, cut_edges)
        s = dual_residual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges)
        push!(r_hist, r)
        push!(s_hist, s)

        ε_pri = eps_primal_branch_pq_socp(partitions, zP, zQ, cut_edges, ϵ_abs, ϵ_rel)
        ε_du = eps_dual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges, ϵ_abs, ϵ_rel)
        push!(eps_r_hist, ε_pri)
        push!(eps_s_hist, ε_du)

        _FAST_USE_LOGS && println("‖r‖ = $r | ε_pri = $ε_pri")
        _FAST_USE_LOGS && println("‖s‖ = $s | ε_du  = $ε_du")

        push!(rhoP_hist, ρP[])
        push!(rhoQ_hist, ρQ[])

        iter_done = k
        if r < ε_pri && s < ε_du
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        if adaptive
            adaptive_rho!(ρP, r, s)
            adaptive_rho!(ρQ, r, s)
        end
    end

    save_socp_branch_pq_plots(r_hist, s_hist, eps_r_hist, eps_s_hist, rhoP_hist, rhoQ_hist;
        stem = "admm_socp_branch")

    total_real_cost = sum(p["real_generation_cost"] for p in values(partitions))
    total_expr_cost = sum(p["generation_cost_model_expr"] for p in values(partitions))
    println("\n[ADMM-SOCP-branch (P,Q)] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println("ρP_final = $(ρP[]) | ρQ_final = $(ρQ[])")
    elapsed_s = (time_ns() - t0_ns) / 1e9
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (SOCP P/Q): $(_pretty_time(elapsed_s))")

    return total_real_cost, partitions
end

# ---- ADMM-SOCP-hybrid (V² + P/Q + R/T) ----
function run_admm_socp_hybrid(config::DistributedOPFModel, data::Dict{String, Any})
    @assert config.base_model == SOCP "ADMM-SOCP-hybrid: base_model deve ser SOCP."

    t0_ns = time_ns()
    iter_done = 0

    partitions, shared_buses = prepare_admm_partitions(data)
    cut_edges = _find_cut_edges(data)

    ρV2 = Ref(config.rho_v)
    ρP = Ref(config.rho_flow_p)
    ρQ = Ref(config.rho_flow_q)
    ρR = Ref(config.rho_flow_p)
    ρT = Ref(config.rho_flow_q)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual
    adaptive = _adaptive_on(config)
    accel = _accel_on(config)    # aceleração só em BUS
    warm = _warm_on(config)

    optimizer_factory = resolve_socp_optimizer(config)

    shared_buses_sorted, b2pos = prepare_dense_maps_socp(shared_buses)

    for (z, sub) in partitions
        m = build_model(sub, SOCP; optimizer = optimizer_factory, build_only = true)
        sub["model"] = m
        sub["obj_base"] = m[:generation_cost]
        sub["v2_var"] = m[:v_magnitude_sq]
        _set_model_silent!(m)

        sub["zone_edges"] = _zone_cut_edges(sub, cut_edges)
        busset = _bus_set(sub)
        sub["shared_buses_local"] = [b for b in shared_buses_sorted if b in busset]
        prepare_socp_hotpaths!(sub, shared_buses_sorted)

        # warm-start de variáveis apenas se warm==true
        if warm
            for I in keys(sub["v2_var"])
                set_start_value(sub["v2_var"][I], 1.0)
            end
        end

        _build_quad_bus_expr_socp!(sub, ρV2[])
        _build_quad_branch_pq_expr!(sub, ρP[], ρQ[], sub["zone_edges"])
        _build_quad_branch_rt_expr!(sub, ρR[], ρT[], sub["zone_edges"])
        sub["obj_static"] =
            sub["obj_base"] + sub["quad_bus_expr"] + sub["quad_pq_expr"] +
            sub["quad_rt_expr"]
        @objective(m, Min, sub["obj_static"])
        optimize!(m)
    end

    # BUS (denso) + EDGES (dict)
    v2bar_vec = zeros(Float64, length(shared_buses_sorted))
    v2bar_prev_vec = similar(v2bar_vec)
    v2bar_hat = similar(v2bar_vec)
    fill!(v2bar_hat, 0.0)
    fista = FISTAStateSOCP()
    duals_v2_dense = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )
    λ_prev_real = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )
    λ_hat = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )
    λ_zone_buf = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )

    zP = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zP_prev = Dict{Int, Float64}()
    zQ = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zQ_prev = Dict{Int, Float64}()
    zR = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zR_prev = Dict{Int, Float64}()
    zT = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    zT_prev = Dict{Int, Float64}()

    dualsP =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))
    dualsQ =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))
    dualsR =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))
    dualsT =
        Dict{Int, Dict{Int, Float64}}(z => Dict{Int, Float64}() for z in keys(partitions))

    r_hist = Float64[]
    s_hist = Float64[]
    eps_r_hist = Float64[]
    eps_s_hist = Float64[]
    rhoV2_hist = Float64[]
    rhoP_hist = Float64[]
    rhoQ_hist = Float64[]
    rhoR_hist = Float64[]
    rhoT_hist = Float64[]

    # ---- warm-start ADMM (somente se warm==true) ----
    if warm
        compute_v2_consensus_vec!(v2bar_vec, partitions, shared_buses_sorted, b2pos)
        v2bar_hat .= v2bar_vec
        for (z, sub) in partitions
            v2s = snapshot_v2_shared!(sub)
            lbs = sub["shared_buses_local"]
            @inbounds for i in eachindex(lbs)
                j = b2pos[lbs[i]]
                duals_v2_dense[z][j] = ρV2[] * (v2s[i] - v2bar_vec[j])
            end
            λ_prev_real[z] .= duals_v2_dense[z]
            λ_hat[z] .= duals_v2_dense[z]
        end
        compute_pq_consensus!(zP, zQ, partitions, cut_edges)
        compute_rt_consensus!(zR, zT, partitions, cut_edges)
        for (z, sub) in partitions
            m = sub["model"]
            for e in cut_edges
                l = get(sub["edge_to_local"], e, nothing)
                if l !== nothing
                    P, Q = _socp_flow_from_exprs_existing(m, l)
                    R, T = _rt_from_exprs_socp(m, l)
                    dualsP[z][e] = ρP[] * (value(P) - zP[e])
                    dualsQ[z][e] = ρQ[] * (value(Q) - zQ[e])
                    dualsR[z][e] = ρR[] * (value(R) - zR[e])
                    dualsT[z][e] = ρT[] * (value(T) - zT[e])
                end
            end
        end
    else
        # flat start seguro quando warm==false
        fill!(v2bar_vec, 1.0)
        v2bar_hat .= v2bar_vec
        # zP/zQ/zR/zT e todos os duais ficam zerados
    end

    prev_ρV2, prev_ρP, prev_ρQ, prev_ρR, prev_ρT = ρV2[], ρP[], ρQ[], ρR[], ρT[]
    zone_ids = collect(keys(partitions))

    for k in 1:config.max_iter
        _FAST_USE_LOGS && println("\n[ADMM-SOCP-hybrid] Iteração $k")

        # rebuild de quads caso ρ mude
        if ρV2[] != prev_ρV2 || ρP[] != prev_ρP || ρQ[] != prev_ρQ || ρR[] != prev_ρR ||
           ρT[] != prev_ρT
            for (_, sub) in partitions
                (ρV2[] != prev_ρV2) && _build_quad_bus_expr_socp!(sub, ρV2[])
                (ρP[] != prev_ρP || ρQ[] != prev_ρQ) &&
                    _build_quad_branch_pq_expr!(sub, ρP[], ρQ[], sub["zone_edges"])
                (ρR[] != prev_ρR || ρT[] != prev_ρT) &&
                    _build_quad_branch_rt_expr!(sub, ρR[], ρT[], sub["zone_edges"])
                sub["obj_static"] =
                    sub["obj_base"] + sub["quad_bus_expr"] + sub["quad_pq_expr"] +
                    sub["quad_rt_expr"]
            end
            prev_ρV2, prev_ρP, prev_ρQ, prev_ρR, prev_ρT = ρV2[], ρP[], ρQ[], ρR[], ρT[]
        end

        α_used = accel ? nothing : _alpha_schedule(k, config.alpha_overrelax; k_on = 10)
        v2_target = accel ? v2bar_hat : v2bar_vec

        nz = length(zone_ids)
        @threads for i in 1:nz
            z = zone_ids[i]
            sub = partitions[z]
            m = sub["model"]
            λz = accel ? (λ_zone_buf[z] .= λ_hat[z]) : (λ_zone_buf[z] .= duals_v2_dense[z])

            lin_bus = _build_linear_bus_expr_socp(
                sub,
                shared_buses_sorted,
                b2pos,
                v2_target,
                λz,
                ρV2[],
            )
            lin_pq = _build_linear_branch_pq_expr(
                sub,
                sub["zone_edges"],
                zP,
                zQ,
                dualsP[z],
                dualsQ[z],
                ρP[],
                ρQ[],
            )
            lin_rt = _build_linear_branch_rt_expr(
                sub,
                sub["zone_edges"],
                zR,
                zT,
                dualsR[z],
                dualsT[z],
                ρR[],
                ρT[],
            )
            @objective(m, Min, sub["obj_static"] + (lin_bus + lin_pq + lin_rt))
            optimize!(m)

            sub["generation_cost_model_expr"] = value(m[:generation_cost])
            sub["real_generation_cost"] =
                compute_real_generation_cost(m, sub["objective_params"])
        end

        # consenso + duais reais
        v2bar_prev_vec .= v2bar_vec
        zP_prev = copy(zP)
        zQ_prev = copy(zQ)
        zR_prev = copy(zR)
        zT_prev = copy(zT)

        compute_v2_consensus_vec!(v2bar_vec, partitions, shared_buses_sorted, b2pos)
        compute_pq_consensus!(zP, zQ, partitions, cut_edges)
        compute_rt_consensus!(zR, zT, partitions, cut_edges)

        if accel
            for z in zone_ids
                λ_prev_real[z] .= duals_v2_dense[z]
            end
        end
        for z in zone_ids
            update_duals_v2_vec!(duals_v2_dense[z], partitions[z], v2bar_vec,
                v2bar_prev_vec,
                ρV2[], α_used, shared_buses_sorted, b2pos)
        end
        for (z, sub) in partitions
            m = sub["model"]
            @inbounds for e in cut_edges
                l = get(sub["edge_to_local"], e, nothing)
                l === nothing && continue
                P, Q = _socp_flow_from_exprs_existing(m, l)
                R, T = _rt_from_exprs_socp(m, l)
                dualsP[z][e] = get(dualsP[z], e, 0.0) + ρP[] * (value(P) - zP[e])
                dualsQ[z][e] = get(dualsQ[z], e, 0.0) + ρQ[] * (value(Q) - zQ[e])
                dualsR[z][e] = get(dualsR[z], e, 0.0) + ρR[] * (value(R) - zR[e])
                dualsT[z][e] = get(dualsT[z], e, 0.0) + ρT[] * (value(T) - zT[e])
            end
        end

        # resíduos / tolerâncias (agregados)
        r_bus =
            primal_residual_bus_socp_vec(partitions, v2bar_vec, shared_buses_sorted, b2pos)
        r_pq = primal_residual_branch_pq_socp(partitions, zP, zQ, cut_edges)
        r_rt = primal_residual_branch_rt_socp(partitions, zR, zT, cut_edges)
        r = sqrt(r_bus^2 + r_pq^2 + r_rt^2)

        s_bus = dual_residual_bus_socp_vec(v2bar_vec, v2bar_prev_vec)
        s_pq = dual_residual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges)
        s_rt = dual_residual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges)
        s = sqrt(s_bus^2 + s_pq^2 + s_rt^2)

        push!(r_hist, r)
        push!(s_hist, s)

        ε_pri = begin
            εb = eps_primal_bus_socp_vec(
                partitions,
                v2bar_vec,
                shared_buses_sorted,
                b2pos,
                ϵ_abs,
                ϵ_rel,
            )
            εpq = eps_primal_branch_pq_socp(partitions, zP, zQ, cut_edges, ϵ_abs, ϵ_rel)
            εrt = eps_primal_branch_rt_socp(partitions, zR, zT, cut_edges, ϵ_abs, ϵ_rel)
            sqrt(εb^2 + εpq^2 + εrt^2)
        end
        ε_du = begin
            εb = eps_dual_bus_socp_vec(v2bar_vec, v2bar_prev_vec, ϵ_abs, ϵ_rel)
            εpq = eps_dual_branch_pq_socp(zP, zQ, zP_prev, zQ_prev, cut_edges, ϵ_abs, ϵ_rel)
            εrt = eps_dual_branch_rt_socp(zR, zT, zR_prev, zT_prev, cut_edges, ϵ_abs, ϵ_rel)
            sqrt(εb^2 + εpq^2 + εrt^2)
        end
        push!(eps_r_hist, ε_pri)
        push!(eps_s_hist, ε_du)

        _FAST_USE_LOGS && println("‖r‖ = $r | ε_pri = $ε_pri")
        _FAST_USE_LOGS && println("‖s‖ = $s | ε_du  = $ε_du")

        push!(rhoV2_hist, ρV2[])
        push!(rhoP_hist, ρP[])
        push!(rhoQ_hist, ρQ[])
        push!(rhoR_hist, ρR[])
        push!(rhoT_hist, ρT[])

        iter_done = k
        if r < ε_pri && s < ε_du
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        if accel
            metric_bus = _accel_metric(r_bus, s_bus)
            restarted = fista_update_socp!(fista; metric = metric_bus)
            if restarted
                v2bar_hat .= v2bar_vec
                for z in zone_ids
                    λ_hat[z] .= duals_v2_dense[z]
                end
            else
                _extrapolate_socp!(v2bar_hat, v2bar_vec, v2bar_prev_vec, fista.α)
                for z in zone_ids
                    _extrapolate_socp!(λ_hat[z], duals_v2_dense[z], λ_prev_real[z], fista.α)
                end
            end
        end

        if adaptive
            adaptive_rho!(ρV2, r_bus, s_bus)
            adaptive_rho!(ρP, r_pq, s_pq)
            adaptive_rho!(ρQ, r_pq, s_pq)
            adaptive_rho!(ρR, r_rt, s_rt)
            adaptive_rho!(ρT, r_rt, s_rt)
        end
    end

    save_socp_hybrid_plots(r_hist, s_hist, eps_r_hist, eps_s_hist,
        rhoV2_hist, rhoP_hist, rhoQ_hist, rhoR_hist, rhoT_hist; stem = "admm_socp_hybrid")

    total_real_cost = sum(p["real_generation_cost"] for p in values(partitions))
    total_expr_cost = sum(p["generation_cost_model_expr"] for p in values(partitions))
    println("\n[ADMM-SOCP-hybrid] Resumo:")
    println("→ Total (value(:generation_cost))     = $total_expr_cost")
    println("→ Total (reavaliado via p_generation) = $total_real_cost")
    println(
        "ρ(V²)_final = $(ρV2[]) | ρP_final = $(ρP[]) | ρQ_final = $(ρQ[]) | ρR_final = $(ρR[]) | ρT_final = $(ρT[])",
    )
    elapsed_s = (time_ns() - t0_ns) / 1e9
    println("Iterações realizadas: $iter_done")
    println("Tempo ADMM (SOCP híbrido): $(_pretty_time(elapsed_s))")

    return total_real_cost, partitions
end
