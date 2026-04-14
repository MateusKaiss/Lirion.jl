# ============================================================
# ALADIN-SOCWR — Phase 2 kernels (dual update, residuals,
# tolerances, adaptive ρ/μ, logging, Hessian refresh hook)

function dual_update!(
    ws::SOCWRStateWorkspaceALADIN,
    nu::@NamedTuple{
        nu_v2::Vector{Float64},
        nu_r::Vector{Float64},
        nu_t::Vector{Float64},
    },
)
    copyto!(ws.lambda_v2_prev, ws.lambda_v2)
    copyto!(ws.lambda_r_prev, ws.lambda_r)
    copyto!(ws.lambda_t_prev, ws.lambda_t)
    a3 = ws.alpha3
    @inbounds for k in eachindex(ws.lambda_v2)
        ws.lambda_v2[k] += a3 * (nu.nu_v2[k] - ws.lambda_v2[k])
    end
    @inbounds for k in eachindex(ws.lambda_r)
        ws.lambda_r[k] += a3 * (nu.nu_r[k] - ws.lambda_r[k])
        ws.lambda_t[k] += a3 * (nu.nu_t[k] - ws.lambda_t[k])
    end
    return nothing
end

function dual_update!(
    ws::SOCWRPowerWorkspaceALADIN,
    nu::@NamedTuple{
        nu_p::Vector{Float64},
        nu_q::Vector{Float64},
    },
)
    copyto!(ws.lambda_p_prev, ws.lambda_p)
    copyto!(ws.lambda_q_prev, ws.lambda_q)
    a3 = ws.alpha3
    @inbounds for k in eachindex(ws.lambda_p)
        ws.lambda_p[k] += a3 * (nu.nu_p[k] - ws.lambda_p[k])
        ws.lambda_q[k] += a3 * (nu.nu_q[k] - ws.lambda_q[k])
    end
    return nothing
end

function dual_update!(
    ws::SOCWRHybridWorkspaceALADIN,
    nu::@NamedTuple{
        nu_v2::Vector{Float64},
        nu_r::Vector{Float64},
        nu_t::Vector{Float64},
        nu_p::Vector{Float64},
        nu_q::Vector{Float64},
    },
)
    copyto!(ws.lambda_v2_prev, ws.lambda_v2)
    copyto!(ws.lambda_r_prev, ws.lambda_r)
    copyto!(ws.lambda_t_prev, ws.lambda_t)
    copyto!(ws.lambda_p_prev, ws.lambda_p)
    copyto!(ws.lambda_q_prev, ws.lambda_q)
    a3 = ws.alpha3
    @inbounds for k in eachindex(ws.lambda_v2)
        ws.lambda_v2[k] += a3 * (nu.nu_v2[k] - ws.lambda_v2[k])
    end
    @inbounds for k in eachindex(ws.lambda_r)
        ws.lambda_r[k] += a3 * (nu.nu_r[k] - ws.lambda_r[k])
        ws.lambda_t[k] += a3 * (nu.nu_t[k] - ws.lambda_t[k])
        ws.lambda_p[k] += a3 * (nu.nu_p[k] - ws.lambda_p[k])
        ws.lambda_q[k] += a3 * (nu.nu_q[k] - ws.lambda_q[k])
    end
    return nothing
end

# ------------------------------------------------------------
# Residuals (HFD 2016 Theorem 4.1)
# ------------------------------------------------------------

function _progress_sq_socwr(zone::AladinZoneDataSOCWR)::Float64
    s = 0.0
    y = zone.y
    x = zone.x_prev
    @inbounds for k in 1:zone.n_vars
        d = y[k] - x[k]
        s += d * d
    end
    return s
end

@inline function _r_pri_bus_sq_socwr(ws)
    s = 0.0
    @inbounds for pair in ws.bus_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        d = zi.y[pair.v2_col_i] - zj.y[pair.v2_col_j]
        s += d * d
    end
    return s
end

@inline function _r_pri_rt_sq_socwr(ws)
    s = 0.0
    @inbounds for pair in ws.edge_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        ri = zi.y[zi.edge_r_cols[pair.edge_local_i]]
        rj = zj.y[zj.edge_r_cols[pair.edge_local_j]]
        ti = zi.y[zi.edge_t_cols[pair.edge_local_i]]
        tj = zj.y[zj.edge_t_cols[pair.edge_local_j]]
        s += (ri - rj)^2 + (ti - tj)^2
    end
    return s
end

@inline function _r_pri_pq_sq_socwr(ws)
    s = 0.0
    @inbounds for pair in ws.edge_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pfi = _socwr_flow_value_at_y(zi, pair.edge_local_i, :p)
        pfj = _socwr_flow_value_at_y(zj, pair.edge_local_j, :p)
        qfi = _socwr_flow_value_at_y(zi, pair.edge_local_i, :q)
        qfj = _socwr_flow_value_at_y(zj, pair.edge_local_j, :q)
        s += (pfi - pfj)^2 + (qfi - qfj)^2
    end
    return s
end

@inline function _sum_progress_sq_socwr(ws)
    s = 0.0
    @inbounds for zone in ws.zones
        s += _progress_sq_socwr(zone)
    end
    return s
end

function compute_residuals!(ws::SOCWRStateWorkspaceALADIN)
    ws.r_primal = sqrt(_r_pri_bus_sq_socwr(ws) + _r_pri_rt_sq_socwr(ws))
    ws.r_progress = sqrt(_sum_progress_sq_socwr(ws))

    lam_sq = 0.0
    @inbounds for k in eachindex(ws.lambda_v2)
        d = ws.lambda_v2[k] - ws.lambda_v2_prev[k]
        lam_sq += d * d
    end
    @inbounds for k in eachindex(ws.lambda_r)
        dr = ws.lambda_r[k] - ws.lambda_r_prev[k]
        dt = ws.lambda_t[k] - ws.lambda_t_prev[k]
        lam_sq += dr * dr + dt * dt
    end
    ws.lambda_step = sqrt(lam_sq)
    return nothing
end

function compute_residuals!(ws::SOCWRPowerWorkspaceALADIN)
    ws.r_primal = sqrt(_r_pri_pq_sq_socwr(ws))
    ws.r_progress = sqrt(_sum_progress_sq_socwr(ws))

    lam_sq = 0.0
    @inbounds for k in eachindex(ws.lambda_p)
        dp = ws.lambda_p[k] - ws.lambda_p_prev[k]
        dq = ws.lambda_q[k] - ws.lambda_q_prev[k]
        lam_sq += dp * dp + dq * dq
    end
    ws.lambda_step = sqrt(lam_sq)
    return nothing
end

function compute_residuals!(ws::SOCWRHybridWorkspaceALADIN)
    ws.r_primal = sqrt(
        _r_pri_bus_sq_socwr(ws) +
        _r_pri_rt_sq_socwr(ws) +
        _r_pri_pq_sq_socwr(ws),
    )
    ws.r_progress = sqrt(_sum_progress_sq_socwr(ws))

    lam_sq = 0.0
    @inbounds for k in eachindex(ws.lambda_v2)
        d = ws.lambda_v2[k] - ws.lambda_v2_prev[k]
        lam_sq += d * d
    end
    @inbounds for k in eachindex(ws.lambda_r)
        dr = ws.lambda_r[k] - ws.lambda_r_prev[k]
        dt = ws.lambda_t[k] - ws.lambda_t_prev[k]
        dp = ws.lambda_p[k] - ws.lambda_p_prev[k]
        dq = ws.lambda_q[k] - ws.lambda_q_prev[k]
        lam_sq += dr * dr + dt * dt + dp * dp + dq * dq
    end
    ws.lambda_step = sqrt(lam_sq)
    return nothing
end

# ------------------------------------------------------------
# Tolerances (Boyd-style)
# ------------------------------------------------------------

const _SOCWRPhase2Workspace = Union{
    SOCWRStateWorkspaceALADIN,
    SOCWRPowerWorkspaceALADIN,
    SOCWRHybridWorkspaceALADIN,
}

function _progress_tol_socwr(
    ws::_SOCWRPhase2Workspace,
    eps_abs::Float64,
    eps_rel::Float64,
)::Float64
    n = 0
    ny_sq = 0.0
    nx_sq = 0.0
    for zone in ws.zones
        n += zone.n_vars
        @inbounds for k in 1:zone.n_vars
            ny_sq += zone.y[k]^2
            nx_sq += zone.x_prev[k]^2
        end
    end
    return sqrt(max(n, 1)) * eps_abs + eps_rel * max(sqrt(ny_sq), sqrt(nx_sq))
end

@inline function _norms_bus_v2_socwr(ws)
    ni = 0.0
    nj = 0.0
    @inbounds for pair in ws.bus_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        vi = zi.y[pair.v2_col_i]
        vj = zj.y[pair.v2_col_j]
        ni += vi * vi
        nj += vj * vj
    end
    return ni, nj
end

@inline function _norms_edge_rt_socwr(ws)
    ni = 0.0
    nj = 0.0
    @inbounds for pair in ws.edge_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        ri = zi.y[zi.edge_r_cols[pair.edge_local_i]]
        rj = zj.y[zj.edge_r_cols[pair.edge_local_j]]
        ti = zi.y[zi.edge_t_cols[pair.edge_local_i]]
        tj = zj.y[zj.edge_t_cols[pair.edge_local_j]]
        ni += ri * ri + ti * ti
        nj += rj * rj + tj * tj
    end
    return ni, nj
end

@inline function _norms_edge_pq_socwr(ws)
    ni = 0.0
    nj = 0.0
    @inbounds for pair in ws.edge_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = _socwr_flow_value_at_y(zi, pair.edge_local_i, :p)
        pj_val = _socwr_flow_value_at_y(zj, pair.edge_local_j, :p)
        qi_val = _socwr_flow_value_at_y(zi, pair.edge_local_i, :q)
        qj_val = _socwr_flow_value_at_y(zj, pair.edge_local_j, :q)
        ni += pi_val * pi_val + qi_val * qi_val
        nj += pj_val * pj_val + qj_val * qj_val
    end
    return ni, nj
end

function compute_tolerances!(
    ws::SOCWRStateWorkspaceALADIN,
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = length(ws.bus_pairs) + 2 * length(ws.edge_pairs)
    ni_b, nj_b = _norms_bus_v2_socwr(ws)
    ni_rt, nj_rt = _norms_edge_rt_socwr(ws)
    ni = ni_b + ni_rt
    nj = nj_b + nj_rt
    ws.eps_primal =
        sqrt(max(n, 1)) * eps_abs +
        eps_rel * max(sqrt(ni), sqrt(nj))
    ws.eps_progress = _progress_tol_socwr(ws, eps_abs, eps_rel)
    return nothing
end

function compute_tolerances!(
    ws::SOCWRPowerWorkspaceALADIN,
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = 2 * length(ws.edge_pairs)
    ni, nj = _norms_edge_pq_socwr(ws)
    ws.eps_primal =
        sqrt(max(n, 1)) * eps_abs +
        eps_rel * max(sqrt(ni), sqrt(nj))
    ws.eps_progress = _progress_tol_socwr(ws, eps_abs, eps_rel)
    return nothing
end

function compute_tolerances!(
    ws::SOCWRHybridWorkspaceALADIN,
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = length(ws.bus_pairs) + 4 * length(ws.edge_pairs)
    ni_b, nj_b = _norms_bus_v2_socwr(ws)
    ni_rt, nj_rt = _norms_edge_rt_socwr(ws)
    ni_pq, nj_pq = _norms_edge_pq_socwr(ws)
    ni = ni_b + ni_rt + ni_pq
    nj = nj_b + nj_rt + nj_pq
    ws.eps_primal =
        sqrt(max(n, 1)) * eps_abs +
        eps_rel * max(sqrt(ni), sqrt(nj))
    ws.eps_progress = _progress_tol_socwr(ws, eps_abs, eps_rel)
    return nothing
end

# ------------------------------------------------------------
# Adaptive ρ/μ
# ------------------------------------------------------------

function adapt_rho_mu!(
    ws::_SOCWRPhase2Workspace,
    params::ALADINParams,
)
    params.adaptive || return nothing
    _grow_ref!(ws.rho, params.rho_growth, params.rho_upper)
    _grow_ref!(ws.mu, params.mu_growth, params.mu_upper)
    return nothing
end

# ------------------------------------------------------------
# Sensitivity refresh hook — called from aladin/core/local_solve.jl
# ------------------------------------------------------------

@inline function refresh_lagrangian_hessian!(
    zone::AladinZoneDataSOCWR,
    ws::_SOCWRPhase2Workspace,
)
    _refresh_hess_diag!(
        zone,
        ws.rho[];
        sigma_theta = ws.sigma_theta,
        sigma_v2 = ws.sigma_v2,
        sigma_r = ws.sigma_r,
        sigma_t = ws.sigma_t,
        sigma_pg = ws.sigma_pg,
        sigma_qg = ws.sigma_qg,
        sigma_other = ws.sigma_other,
    )
    return nothing
end

# ------------------------------------------------------------
# Logging tags
# ------------------------------------------------------------

_aladin_print_tag(::SOCWRStateWorkspaceALADIN) = "[ALADIN-state/SOCWR]"
_aladin_print_tag(::SOCWRPowerWorkspaceALADIN) = "[ALADIN-power/SOCWR]"
_aladin_print_tag(::SOCWRHybridWorkspaceALADIN) = "[ALADIN-hybrid/SOCWR]"

# ------------------------------------------------------------
# Armijo line search on the HFD 2016 merit function
# ------------------------------------------------------------
# HFD 2016 Section 5 globalizes ALADIN by line-searching the primal
# step on the augmented Lagrangian merit function
#
#     Φ(x, λ) = Σ f_i(x_i) + λᵀ r(x) + (μ/2) ‖r(x)‖²
#
# where r(x) = Σ A_i x_i is the consensus violation (the same quantity
# that gives ws.r_primal). We search along the Newton direction
# d_i = Δy_i obtained from the coordinator QP, holding α₁ = α₃ = 1
# and backtracking α₂ ∈ (0, 1] until the Armijo sufficient-decrease
# condition
#
#     Φ(y + α·d) ≤ Φ(y) + c₁·α·Φ'(y;d)
#
# is satisfied. The directional derivative Φ'(y;d) is estimated by a
# forward difference at ε = 1e-6, which is exact for the quadratic
# pieces (cost + μ/2 ‖r‖²) and accurate-to-machine for the affine
# piece (λᵀr).
#
# With the cone active-face rows now present in the active Jacobian,
# the canonical Newton direction is descent and α₂ = 1 is typically
# accepted on the first try. The line search thus acts as a safety
# net for tight cases rather than as a primary globalization, and
# removes the need for a user-specified damping hyperparameter.

@inline function _socwr_cost_at_trial(zone::AladinZoneDataSOCWR, α::Float64)
    acc = 0.0
    y = zone.y
    delta = zone.delta
    @inbounds for k in eachindex(zone.pg_cols)
        col = zone.pg_cols[k]
        xk = y[col] + α * delta[col]
        acc += 0.5 * zone.pg_c2_sbase2[k] * xk * xk + zone.pg_c1_sbase[k] * xk
    end
    return acc
end

@inline function _socwr_flow_at_trial(
    zone::AladinZoneDataSOCWR,
    le::Int,
    α::Float64,
    which::Symbol,
)::Float64
    cols = zone.edge_flow_cols[le]
    coefs = which === :p ? zone.edge_pflow_coefs[le] : zone.edge_qflow_coefs[le]
    const_val = which === :p ? zone.edge_pflow_const[le] : zone.edge_qflow_const[le]
    acc = const_val
    y = zone.y
    delta = zone.delta
    @inbounds for i in eachindex(cols)
        c = cols[i]
        acc += coefs[i] * (y[c] + α * delta[c])
    end
    return acc
end

function _socwr_merit(ws::SOCWRStateWorkspaceALADIN, α::Float64)
    Φ = 0.0
    for zone in ws.zones
        Φ += _socwr_cost_at_trial(zone, α)
    end

    r_sq = 0.0
    λr = 0.0
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        xi = zi.y[pair.v2_col_i] + α * zi.delta[pair.v2_col_i]
        xj = zj.y[pair.v2_col_j] + α * zj.delta[pair.v2_col_j]
        d = xi - xj
        r_sq += d * d
        λr += ws.lambda_v2[k] * d
    end
    @inbounds for (k, pair) in enumerate(ws.edge_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        ri_col = zi.edge_r_cols[pair.edge_local_i]
        rj_col = zj.edge_r_cols[pair.edge_local_j]
        ti_col = zi.edge_t_cols[pair.edge_local_i]
        tj_col = zj.edge_t_cols[pair.edge_local_j]
        dr = (zi.y[ri_col] + α * zi.delta[ri_col]) -
             (zj.y[rj_col] + α * zj.delta[rj_col])
        dt = (zi.y[ti_col] + α * zi.delta[ti_col]) -
             (zj.y[tj_col] + α * zj.delta[tj_col])
        r_sq += dr * dr + dt * dt
        λr += ws.lambda_r[k] * dr + ws.lambda_t[k] * dt
    end

    mu = ws.mu[]
    return Φ + λr + 0.5 * mu * r_sq
end

function _socwr_merit(ws::SOCWRPowerWorkspaceALADIN, α::Float64)
    Φ = 0.0
    for zone in ws.zones
        Φ += _socwr_cost_at_trial(zone, α)
    end

    r_sq = 0.0
    λr = 0.0
    @inbounds for (k, pair) in enumerate(ws.edge_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = _socwr_flow_at_trial(zi, pair.edge_local_i, α, :p)
        pj_val = _socwr_flow_at_trial(zj, pair.edge_local_j, α, :p)
        qi_val = _socwr_flow_at_trial(zi, pair.edge_local_i, α, :q)
        qj_val = _socwr_flow_at_trial(zj, pair.edge_local_j, α, :q)
        dp = pi_val - pj_val
        dq = qi_val - qj_val
        r_sq += dp * dp + dq * dq
        λr += ws.lambda_p[k] * dp + ws.lambda_q[k] * dq
    end

    mu = ws.mu[]
    return Φ + λr + 0.5 * mu * r_sq
end

function _socwr_merit(ws::SOCWRHybridWorkspaceALADIN, α::Float64)
    Φ = 0.0
    for zone in ws.zones
        Φ += _socwr_cost_at_trial(zone, α)
    end

    r_sq = 0.0
    λr = 0.0
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        xi = zi.y[pair.v2_col_i] + α * zi.delta[pair.v2_col_i]
        xj = zj.y[pair.v2_col_j] + α * zj.delta[pair.v2_col_j]
        d = xi - xj
        r_sq += d * d
        λr += ws.lambda_v2[k] * d
    end
    @inbounds for (k, pair) in enumerate(ws.edge_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        ri_col = zi.edge_r_cols[pair.edge_local_i]
        rj_col = zj.edge_r_cols[pair.edge_local_j]
        ti_col = zi.edge_t_cols[pair.edge_local_i]
        tj_col = zj.edge_t_cols[pair.edge_local_j]
        dr = (zi.y[ri_col] + α * zi.delta[ri_col]) -
             (zj.y[rj_col] + α * zj.delta[rj_col])
        dt = (zi.y[ti_col] + α * zi.delta[ti_col]) -
             (zj.y[tj_col] + α * zj.delta[tj_col])

        pi_val = _socwr_flow_at_trial(zi, pair.edge_local_i, α, :p)
        pj_val = _socwr_flow_at_trial(zj, pair.edge_local_j, α, :p)
        qi_val = _socwr_flow_at_trial(zi, pair.edge_local_i, α, :q)
        qj_val = _socwr_flow_at_trial(zj, pair.edge_local_j, α, :q)
        dp = pi_val - pj_val
        dq = qi_val - qj_val

        r_sq += dr * dr + dt * dt + dp * dp + dq * dq
        λr +=
            ws.lambda_r[k] * dr + ws.lambda_t[k] * dt +
            ws.lambda_p[k] * dp + ws.lambda_q[k] * dq
    end

    mu = ws.mu[]
    return Φ + λr + 0.5 * mu * r_sq
end

# ------------------------------------------------------------
# Line search hook — called between QP solve and primal_update!
# ------------------------------------------------------------

function _line_search_step!(
    ws::_SOCWRPhase2Workspace;
    c1::Float64 = 1e-4,
    ε_fd::Float64 = 1e-6,
    max_backtrack::Int = 25,
    α_min::Float64 = 1e-3,
)
    ws.line_search || return nothing

    Φ_0 = _socwr_merit(ws, 0.0)
    Φ_ε = _socwr_merit(ws, ε_fd)
    dΦ_0 = (Φ_ε - Φ_0) / ε_fd

    # If the QP direction is not descent, fall back to α_min to keep
    # the iteration making progress via the strong prox pull.
    if !isfinite(dΦ_0) || dΦ_0 >= 0.0
        ws.alpha2 = α_min
        return nothing
    end

    α = ws.alpha2_init
    @inbounds for _ in 1:max_backtrack
        Φ_α = _socwr_merit(ws, α)
        if isfinite(Φ_α) && Φ_α <= Φ_0 + c1 * α * dΦ_0
            ws.alpha2 = α
            return nothing
        end
        α *= 0.5
        if α < α_min
            break
        end
    end
    ws.alpha2 = α_min
    return nothing
end
