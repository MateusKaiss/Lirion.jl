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
    model::JuMP.Model,
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

# ---------- Alpha schedule / aceleração ----------
@inline function _alpha_schedule_ac(k::Int, α_cfg::Union{Nothing, Float64}; k_on::Int = 10)
    if α_cfg === nothing || !(1.0 < α_cfg < 2.0)
        return nothing
    end
    return k >= k_on ? α_cfg : nothing
end

mutable struct FISTAStateAC
    β::Float64
    α::Float64
    prev_metric::Float64
end
FISTAStateAC() = FISTAStateAC(1.0, 0.0, Inf)

function fista_update!(st::FISTAStateAC; metric::Float64)
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

@inline function _extrapolateac!(
    y_hat::AbstractVector{T},
    y::AbstractVector{T},
    y_prev::AbstractVector{T},
    γ::T,
) where {T <: Real}
    @inbounds @simd for i in eachindex(y)
        y_hat[i] = y[i] + γ * (y[i] - y_prev[i])
    end
    return nothing
end

# ---------- Warm start (resolve base e puxa θ e |V|) ----------
function _warm_start_ac_base!(sub::Dict{String, Any})
    m = sub["model"]
    θ = m[:theta]
    V = m[:v_magnitude]

    @objective(m, Min, sub["obj_base"])
    optimize!(m)

    if !JuMP.has_values(m)
        @warn "Warm start (AC) ignorado: optimize! não retornou valores (status=$(JuMP.termination_status(m)))."
        return
    end
    θ_vals = Dict{typeof(first(keys(θ))), Float64}()
    for I in keys(θ)
        θ_vals[I] = JuMP.value(θ[I])
    end
    V_vals = Dict{typeof(first(keys(V))), Float64}()
    for I in keys(V)
        V_vals[I] = JuMP.value(V[I])
    end

    for I in keys(θ)
        set_start_value(θ[I], θ_vals[I])
    end
    for I in keys(V)
        set_start_value(V[I], V_vals[I])
    end
    return nothing
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
    sub["rho_bus_cached"] = (ρθ, ρV)
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
    sub["rho_branch_cached"] = (ρP, ρQ)
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

function compute_ac_branch_consensus!(
    zP::Dict{Int, Float64},
    zQ::Dict{Int, Float64},
    partitions::Dict{Int, Dict{String, Any}},
    cut_edges::Vector{Int},
)
    sumsP = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    sumsQ = Dict{Int, Float64}(e => 0.0 for e in cut_edges)
    cntP = Dict{Int, Int}(e => 0 for e in cut_edges)
    cntQ = Dict{Int, Int}(e => 0 for e in cut_edges)

    for (_, sub) in partitions
        m = sub["model"]
        for e in get(sub, "zone_edges", Int[])
            l_local = get(sub["edge_to_local"], e, nothing)
            l_local === nothing && continue
            P_loc, Q_loc = _ac_flow_from_exprs_existing(m, sub, l_local)
            pval = JuMP.value(P_loc)
            qval = JuMP.value(Q_loc)
            if isfinite(pval)
                sumsP[e] += pval
                cntP[e] += 1
            end
            if isfinite(qval)
                sumsQ[e] += qval
                cntQ[e] += 1
            end
        end
    end

    for e in cut_edges
        if cntP[e] > 0
            zP[e] = sumsP[e] / cntP[e]
        end
        if cntQ[e] > 0
            zQ[e] = sumsQ[e] / cntQ[e]
        end
    end
    return nothing
end

@inline function update_duals_bus_ac_vec!(
    dualsθ_vec_zone::Vector{Float64},
    dualsV_vec_zone::Vector{Float64},
    sub::Dict{String, Any},
    θ̄_vec::Vector{Float64},
    V̄_vec::Vector{Float64},
    ρθ::Float64,
    ρV::Float64,
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int};
    α::Union{Nothing, Float64} = nothing,
    θ̄_prev::Union{Nothing, Vector{Float64}} = nothing,
    V̄_prev::Union{Nothing, Vector{Float64}} = nothing,
)
    θs = snapshot_theta_ac_shared!(sub)
    Vs = snapshot_vmag_shared!(sub)
    lbs = sub["shared_buses_local"]

    if α === nothing
        @inbounds for i in eachindex(lbs)
            j = b2pos[lbs[i]]
            dθ = _wrap_angle_diff(θs[i] - θ̄_vec[j])
            dV = (Vs[i] - V̄_vec[j])
            dualsθ_vec_zone[j] += ρθ * dθ
            dualsV_vec_zone[j] += ρV * dV
        end
    else
        αv = something(α, 1.0)
        @inbounds for i in eachindex(lbs)
            j = b2pos[lbs[i]]
            θrel = _wrap_angle_diff(αv * θs[i] + (1 - αv) * θ̄_prev[j] - θ̄_vec[j])
            Vrel = αv * Vs[i] + (1 - αv) * V̄_prev[j] - V̄_vec[j]
            dualsθ_vec_zone[j] += ρθ * θrel
            dualsV_vec_zone[j] += ρV * Vrel
        end
    end
    return nothing
end

function update_duals_branch_ac!(
    dualsP::Dict{Int, Dict{Int, Float64}},
    dualsQ::Dict{Int, Dict{Int, Float64}},
    partitions::Dict{Int, Dict{String, Any}},
    zP::Dict{Int, Float64},
    zQ::Dict{Int, Float64},
    cut_edges::Vector{Int},
    ρP::Float64,
    ρQ::Float64;
    α::Union{Nothing, Float64} = nothing,
    zP_prev::Union{Nothing, Dict{Int, Float64}} = nothing,
    zQ_prev::Union{Nothing, Dict{Int, Float64}} = nothing,
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

            if α === nothing
                if isfinite(pval)
                    λPz[e] = get(λPz, e, 0.0) + ρP * (pval - get(zP, e, 0.0))
                end
                if isfinite(qval)
                    λQz[e] = get(λQz, e, 0.0) + ρQ * (qval - get(zQ, e, 0.0))
                end
            else
                αv = something(α, 1.0)
                if isfinite(pval)
                    p_rel = αv * pval + (1 - αv) * get(zP_prev, e, get(zP, e, 0.0))
                    λPz[e] = get(λPz, e, 0.0) + ρP * (p_rel - get(zP, e, 0.0))
                end
                if isfinite(qval)
                    q_rel = αv * qval + (1 - αv) * get(zQ_prev, e, get(zQ, e, 0.0))
                    λQz[e] = get(λQz, e, 0.0) + ρQ * (q_rel - get(zQ, e, 0.0))
                end
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

function primal_residual_bus_ac_components(
    partitions::Dict{Int, Dict{String, Any}},
    θ̄_vec::Vector{Float64},
    V̄_vec::Vector{Float64},
    shared_buses_sorted::Vector{Int},
    b2pos::Dict{Int, Int},
)
    rθ = 0.0
    rV = 0.0
    for (_, sub) in partitions
        θbuf = snapshot_theta_ac!(sub)
        Vbuf = snapshot_vmag!(sub)
        pos = sub["bus_local_index"]
        @inbounds for b in sub["shared_buses_local"]
            j = b2pos[b]
            rθ += _wrap_angle_diff(θbuf[pos[b]] - θ̄_vec[j])^2
            rV += (Vbuf[pos[b]] - V̄_vec[j])^2
        end
    end
    return sqrt(rθ), sqrt(rV)
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

function dual_residual_bus_ac_components(
    θ̄_vec::Vector{Float64}, V̄_vec::Vector{Float64},
    θ̄_prev::Vector{Float64}, V̄_prev::Vector{Float64},
)
    sθ = 0.0
    sV = 0.0
    @inbounds for j in eachindex(θ̄_vec)
        sθ += _wrap_angle_diff(θ̄_vec[j] - θ̄_prev[j])^2
        sV += (V̄_vec[j] - V̄_prev[j])^2
    end
    return sqrt(sθ), sqrt(sV)
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
    dim = 0
    for (_, sub) in partitions
        dim += 2 * length(get(sub, "zone_edges", Int[]))
    end
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
            ))(),
        "plot_rho_theta_$(stem).html")
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
            ))(),
        "plot_rho_vmag_$(stem).html")
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
            ))(),
        "plot_rhoP_$(stem).html")
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
            ))(),
        "plot_rhoQ_$(stem).html")
end

# ============================ ADMM: BUS (AC) ============================
function run_admm_ac(config::DistributedOPFModel, data::Dict{String, Any})
    t0 = time_ns()
    @assert config.base_model == AC "ADMM-AC: base_model deve ser AC."
    partitions, shared_buses = prepare_admm_partitions(data)

    ρθ = Ref(config.rho_theta)
    ρV = Ref(config.rho_v)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual
    adaptive = _adaptive_on(config)
    accel = _accel_on(config)

    if accel
        adaptive = false
    end

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

    θ̄_hat = similar(θ̄_vec)
    fill!(θ̄_hat, 0.0)
    V̄_hat = similar(V̄_vec)
    fill!(V̄_hat, 1.0)
    fista = FISTAStateAC()
    λθ_prev_real = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )
    λV_prev_real = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )
    λθ_hat = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )
    λV_hat = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )

    res_pri_hist = Float64[]
    res_du_hist = Float64[]
    eps_pri_hist = Float64[]
    eps_du_hist = Float64[]
    rhoθ_hist = Float64[]
    rhoV_hist = Float64[]
    θ_avg_hist = Dict{Int, Vector{Float64}}()

    optimizer_factory = resolve_optimizer_ac(config)

    warm = try
        getfield(config, :warm_start)
    catch
        false
    end

    # build inicial + hotpaths + pré-quad
    for (zone_id, sub) in partitions
        m = build_model(sub, AC; optimizer = optimizer_factory, build_only = true)
        sub["model"] = m
        sub["obj_base"] = m[:generation_cost]
        sub["theta_var"] = m[:theta]
        sub["vmag_var"] = m[:v_magnitude]
        sub["has_solution"] = false
        _set_model_silent!(m)

        for I in keys(sub["theta_var"])
            set_start_value(sub["theta_var"][I], 0.0)
        end
        for I in keys(sub["vmag_var"])
            set_start_value(sub["vmag_var"][I], 1.0)
        end

        prepare_hotpaths_ac!(sub, shared_buses)
        _build_quad_bus_expr_ac!(sub, ρθ[], ρV[])
        @objective(m, Min, sub["obj_base"] + sub["quad_bus_expr_ac"])
    end

    if warm
        for (_, sub) in partitions
            _warm_start_ac_base!(sub)
            @objective(sub["model"], Min, sub["obj_base"] + sub["quad_bus_expr_ac"])
        end
        compute_ac_bus_consensus_vec!(
            θ̄_vec,
            V̄_vec,
            partitions,
            shared_buses_sorted,
            b2pos,
        )
        θ̄_hat .= θ̄_vec
        V̄_hat .= V̄_vec
        for z in keys(partitions)
            θs = snapshot_theta_ac_shared!(partitions[z])
            Vs = snapshot_vmag_shared!(partitions[z])
            lbs = partitions[z]["shared_buses_local"]
            @inbounds for i in eachindex(lbs)
                j = b2pos[lbs[i]]
                dualsθ_dense[z][j] = ρθ[] * _wrap_angle_diff(θs[i] - θ̄_vec[j])
                dualsV_dense[z][j] = ρV[] * (Vs[i] - V̄_vec[j])
            end
            λθ_prev_real[z] .= dualsθ_dense[z]
            λθ_hat[z] .= dualsθ_dense[z]
            λV_prev_real[z] .= dualsV_dense[z]
            λV_hat[z] .= dualsV_dense[z]
        end
    end

    prev_ρθ, prev_ρV = ρθ[], ρV[]
    zone_ids = collect(keys(partitions))
    λθ_zone_buf = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in zone_ids
    )
    λV_zone_buf = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in zone_ids
    )

    for k in 1:config.max_iter
        _FAST_USE_LOGS && println("\n[ADMM-AC] Iteração $k")

        if ρθ[] != prev_ρθ || ρV[] != prev_ρV
            for (_, sub) in partitions
                _build_quad_bus_expr_ac!(sub, ρθ[], ρV[])
            end
            prev_ρθ, prev_ρV = ρθ[], ρV[]
        end

        α_used = accel ? nothing : _alpha_schedule_ac(k, config.alpha_overrelax; k_on = 10)

        θ̄_target = accel ? θ̄_hat : θ̄_vec
        V̄_target = accel ? V̄_hat : V̄_vec

        nz = length(zone_ids)
        results = Vector{Any}(undef, nz)
        Threads.@threads for i in 1:nz
            z = zone_ids[i]
            sub = partitions[z]
            λθ = accel ? (λθ_zone_buf[z] .= λθ_hat[z]) : (λθ_zone_buf[z] .= dualsθ_dense[z])
            λV = accel ? (λV_zone_buf[z] .= λV_hat[z]) : (λV_zone_buf[z] .= dualsV_dense[z])
            lin = _build_linear_bus_expr_ac!(sub, shared_buses_sorted, b2pos,
                θ̄_target, V̄_target, λθ, λV, ρθ[], ρV[])
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

        if accel
            for z in zone_ids
                λθ_prev_real[z] .= dualsθ_dense[z]
                λV_prev_real[z] .= dualsV_dense[z]
            end
        end

        for z in zone_ids
            update_duals_bus_ac_vec!(dualsθ_dense[z], dualsV_dense[z], partitions[z],
                θ̄_vec, V̄_vec, ρθ[], ρV[], shared_buses_sorted, b2pos;
                α = α_used, θ̄_prev = θ̄_prev, V̄_prev = V̄_prev)
        end

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

        rθ, rV = primal_residual_bus_ac_components(
            partitions,
            θ̄_vec,
            V̄_vec,
            shared_buses_sorted,
            b2pos,
        )
        sθ, sV = dual_residual_bus_ac_components(θ̄_vec, V̄_vec, θ̄_prev, V̄_prev)

        if r < εp && s < εd
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        if accel
            metric = _accel_metric(r, s)
            restarted = fista_update!(fista; metric = metric)
            if restarted
                θ̄_hat .= θ̄_vec
                V̄_hat .= V̄_vec
                for z in zone_ids
                    λθ_hat[z] .= dualsθ_dense[z]
                    λV_hat[z] .= dualsV_dense[z]
                end
            else
                _extrapolateac!(θ̄_hat, θ̄_vec, θ̄_prev, fista.α)
                _extrapolateac!(V̄_hat, V̄_vec, V̄_prev, fista.α)
                for z in zone_ids
                    _extrapolateac!(λθ_hat[z], dualsθ_dense[z], λθ_prev_real[z], fista.α)
                    _extrapolateac!(λV_hat[z], dualsV_dense[z], λV_prev_real[z], fista.α)
                end
            end
        end

        if adaptive
            adaptive_rho_ac!(ρθ, rθ, sθ; μ_incr = config.mu_incr, μ_decr = config.mu_decr,
                τ_incr = config.tau_incr, τ_decr = config.tau_decr)
            adaptive_rho_ac!(ρV, rV, sV; μ_incr = config.mu_incr, μ_decr = config.mu_decr,
                τ_incr = config.tau_incr, τ_decr = config.tau_decr)
        end
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
function run_admm_ac_branch(
    config::DistributedOPFModel,
    data::Dict{String, Any};
    warm_start::Bool = false,
)
    t0 = time_ns()
    @assert config.base_model == AC "ADMM-AC-branch: base_model deve ser AC."
    partitions, _ = prepare_admm_partitions(data)
    cut_edges = _find_cut_edges(data)

    ρP = Ref(config.rho_flow_p)
    ρQ = Ref(config.rho_flow_q)
    ϵ_abs, ϵ_rel = config.ε_primal, config.ε_dual
    adaptive = _adaptive_on(config)

    warm_cfg = try
        getfield(config, :warm_start)
    catch
        false
    end
    warm = warm_start || warm_cfg

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

        for I in keys(sub["theta_var"])
            set_start_value(sub["theta_var"][I], 0.0)
        end
        for I in keys(sub["vmag_var"])
            set_start_value(sub["vmag_var"][I], 1.0)
        end

        _build_quad_branch_expr_ac!(sub, ρP[], ρQ[], sub["zone_edges"])
        @objective(m, Min, sub["obj_base"] + sub["quad_branch_expr_ac"])
    end

    if warm
        for (_, sub) in partitions
            _warm_start_ac_base!(sub)
            @objective(sub["model"], Min, sub["obj_base"] + sub["quad_branch_expr_ac"])
        end
        compute_ac_branch_consensus!(zP, zQ, partitions, cut_edges)
        for (z, sub) in partitions
            m = sub["model"]
            for e in cut_edges
                l_local = get(sub["edge_to_local"], e, nothing)
                l_local === nothing && continue
                P_loc, Q_loc = _ac_flow_from_exprs_existing(m, sub, l_local)
                pval = JuMP.value(P_loc)
                qval = JuMP.value(Q_loc)
                if isfinite(pval)
                    dualsP[z][e] = ρP[] * (pval - get(zP, e, 0.0))
                end
                if isfinite(qval)
                    dualsQ[z][e] = ρQ[] * (qval - get(zQ, e, 0.0))
                end
            end
        end
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

        α_used = _alpha_schedule_ac(k, config.alpha_overrelax; k_on = 10)

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
        update_duals_branch_ac!(dualsP, dualsQ, partitions, zP, zQ, cut_edges, ρP[], ρQ[];
            α = α_used, zP_prev = zP_prev, zQ_prev = zQ_prev)

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

        if adaptive
            adaptive_rho_ac!(ρP, r, s; μ_incr = config.mu_incr, μ_decr = config.mu_decr,
                τ_incr = config.tau_incr, τ_decr = config.tau_decr)
            adaptive_rho_ac!(ρQ, r, s; μ_incr = config.mu_incr, μ_decr = config.mu_decr,
                τ_incr = config.tau_incr, τ_decr = config.tau_decr)
        end
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
    adaptive = _adaptive_on(config)
    accel = _accel_on(config)

    adaptive_bus = adaptive && !accel
    adaptive_br = adaptive

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

    θ̄_hat = similar(θ̄_vec)
    fill!(θ̄_hat, 0.0)
    V̄_hat = similar(V̄_vec)
    fill!(V̄_hat, 1.0)
    fista = FISTAStateAC()
    λθ_prev_real = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )
    λV_prev_real = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )
    λθ_hat = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )
    λV_hat = Dict{Int, Vector{Float64}}(
        z => zeros(Float64, length(shared_buses_sorted)) for z in keys(partitions)
    )

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

    warm = try
        getfield(config, :warm_start)
    catch
        false
    end

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

        for I in keys(sub["theta_var"])
            set_start_value(sub["theta_var"][I], 0.0)
        end
        for I in keys(sub["vmag_var"])
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

    if warm
        for (_, sub) in partitions
            _warm_start_ac_base!(sub)
            @objective(
                sub["model"],
                Min,
                sub["obj_base"] + sub["quad_bus_expr_ac"] + sub["quad_branch_expr_ac"]
            )
        end
        compute_ac_bus_consensus_vec!(
            θ̄_vec,
            V̄_vec,
            partitions,
            shared_buses_sorted,
            b2pos,
        )
        θ̄_hat .= θ̄_vec
        V̄_hat .= V̄_vec
        for z in keys(partitions)
            θs = snapshot_theta_ac_shared!(partitions[z])
            Vs = snapshot_vmag_shared!(partitions[z])
            lbs = partitions[z]["shared_buses_local"]
            @inbounds for i in eachindex(lbs)
                j = b2pos[lbs[i]]
                dualsθ_dense[z][j] = ρθ[] * _wrap_angle_diff(θs[i] - θ̄_vec[j])
                dualsV_dense[z][j] = ρV[] * (Vs[i] - V̄_vec[j])
            end
            λθ_prev_real[z] .= dualsθ_dense[z]
            λθ_hat[z] .= dualsθ_dense[z]
            λV_prev_real[z] .= dualsV_dense[z]
            λV_hat[z] .= dualsV_dense[z]
        end
        compute_ac_branch_consensus!(zP, zQ, partitions, cut_edges)
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

        α_used = accel ? nothing : _alpha_schedule_ac(k, config.alpha_overrelax; k_on = 10)
        θ̄_target = accel ? θ̄_hat : θ̄_vec
        V̄_target = accel ? V̄_hat : V̄_vec

        nz = length(zone_ids)
        results = Vector{Any}(undef, nz)
        Threads.@threads for i in 1:nz
            z = zone_ids[i]
            sub = partitions[z]
            λθ = accel ? (λθ_zone_buf[z] .= λθ_hat[z]) : (λθ_zone_buf[z] .= dualsθ_dense[z])
            λV = accel ? (λV_zone_buf[z] .= λV_hat[z]) : (λV_zone_buf[z] .= dualsV_dense[z])

            lin_bus = _build_linear_bus_expr_ac!(
                sub,
                shared_buses_sorted,
                b2pos,
                θ̄_target,
                V̄_target,
                λθ,
                λV,
                ρθ[],
                ρV[],
            )
            lin_br = _build_linear_branch_expr_ac!(
                sub,
                sub["zone_edges"],
                zP,
                zQ,
                dualsP[z],
                dualsQ[z],
                ρP[],
                ρQ[],
            )

            @objective(
                sub["model"],
                Min,
                sub["obj_base"] + sub["quad_bus_expr_ac"] + sub["quad_branch_expr_ac"] +
                (lin_bus + lin_br)
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
                θ̄_vec, V̄_vec, ρθ[], ρV[], shared_buses_sorted, b2pos;
                α = α_used, θ̄_prev = θ̄_prev, V̄_prev = V̄_prev)
        end
        update_duals_branch_ac!(dualsP, dualsQ, partitions, zP, zQ, cut_edges, ρP[], ρQ[];
            α = α_used, zP_prev = zP_prev, zQ_prev = zQ_prev)

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

        rθ, rV = primal_residual_bus_ac_components(
            partitions,
            θ̄_vec,
            V̄_vec,
            shared_buses_sorted,
            b2pos,
        )
        sθ, sV = dual_residual_bus_ac_components(θ̄_vec, V̄_vec, θ̄_prev, V̄_prev)

        if r_tot < ε_pri && s_tot < ε_du
            _FAST_USE_LOGS && println("Stopping criterion met.")
            break
        end

        if accel
            metric_bus = _accel_metric(r_bus, s_bus)
            restarted = fista_update!(fista; metric = metric_bus)
            if restarted
                θ̄_hat .= θ̄_vec
                V̄_hat .= V̄_vec
                for z in zone_ids
                    λθ_hat[z] .= dualsθ_dense[z]
                    λV_hat[z] .= dualsV_dense[z]
                end
            else
                _extrapolateac!(θ̄_hat, θ̄_vec, θ̄_prev, fista.α)
                _extrapolateac!(V̄_hat, V̄_vec, V̄_prev, fista.α)
                for z in zone_ids
                    _extrapolateac!(λθ_hat[z], dualsθ_dense[z], λθ_prev_real[z], fista.α)
                    _extrapolateac!(λV_hat[z], dualsV_dense[z], λV_prev_real[z], fista.α)
                end
            end
        end

        if adaptive_bus
            adaptive_rho_ac!(ρθ, rθ, sθ; μ_incr = config.mu_incr,
                μ_decr = config.mu_decr,
                τ_incr = config.tau_incr, τ_decr = config.tau_decr)
            adaptive_rho_ac!(ρV, rV, sV; μ_incr = config.mu_incr,
                μ_decr = config.mu_decr,
                τ_incr = config.tau_incr, τ_decr = config.tau_decr)
        end
        if adaptive_br
            adaptive_rho_ac!(ρP, r_br, s_br; μ_incr = config.mu_incr,
                μ_decr = config.mu_decr,
                τ_incr = config.tau_incr, τ_decr = config.tau_decr)
            adaptive_rho_ac!(ρQ, r_br, s_br; μ_incr = config.mu_incr,
                μ_decr = config.mu_decr,
                τ_incr = config.tau_incr, τ_decr = config.tau_decr)
        end
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
