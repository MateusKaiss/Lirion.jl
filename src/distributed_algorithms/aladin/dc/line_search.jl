# ============================================================
# ALADIN-DC — Armijo line search on the HFD 2016 merit function
# ============================================================
#
# Implements Algorithm 3 of Houska-Frasch-Diehl 2016 (p. 1112-1113).
# DC OPF is convex, so the local Lemma 3 / Theorem 4 hypotheses are
# satisfied without globalisation; line search is included for
# robustness against non-quadratic descent and to bring DC parity with
# SOCWR. Set `line_search = false` in the workspace to disable.
#
# Φ(α) = Σ cost_i(y + αδ) + λᵀ·r(α) + (μ/2)·‖r(α)‖²
#
# r(α) is the consensus residual (affine in α since DC fluxos = b·Δθ),
# so dΦ/dα|_{α=0} = Σ ∇cost_i(y)·δ + (λ + μ·r₀)ᵀ·Δr is exact.

const _DCPhase2Workspace = Union{
    DCStateWorkspaceALADIN,
    DCPowerWorkspaceALADIN,
    DCHybridWorkspaceALADIN,
}

# ------------------------------------------------------------
# Cost evaluation at trial point
# ------------------------------------------------------------

@inline function _dc_cost_at_trial(zone::AladinZoneData, α::Float64)
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

@inline function _dc_cost_dir_deriv_0(zone::AladinZoneData)
    acc = 0.0
    y = zone.y
    delta = zone.delta
    @inbounds for k in eachindex(zone.pg_cols)
        col = zone.pg_cols[k]
        acc += (zone.pg_c2_sbase2[k] * y[col] + zone.pg_c1_sbase[k]) * delta[col]
    end
    return acc
end

# ------------------------------------------------------------
# Helpers for affine flow consensus
# ------------------------------------------------------------

@inline function _dc_pflow_at(zone::AladinZoneData, le::Int, α::Float64)
    fcol = zone.edge_from_cols[le]
    tcol = zone.edge_to_cols[le]
    bval = zone.edge_b[le]
    yf = zone.y[fcol] + α * zone.delta[fcol]
    yt = zone.y[tcol] + α * zone.delta[tcol]
    return bval * (yf - yt)
end

@inline function _dc_pflow_dir_deriv(zone::AladinZoneData, le::Int)
    fcol = zone.edge_from_cols[le]
    tcol = zone.edge_to_cols[le]
    bval = zone.edge_b[le]
    return bval * (zone.delta[fcol] - zone.delta[tcol])
end

@inline function _dc_dummy_pflow_at(
    zone::AladinZoneData,
    fcol::Int, tcol::Int, bval::Float64,
    α::Float64,
)
    yf = zone.y[fcol] + α * zone.delta[fcol]
    yt = zone.y[tcol] + α * zone.delta[tcol]
    return bval * (yf - yt)
end

@inline function _dc_dummy_pflow_dir_deriv(
    zone::AladinZoneData,
    fcol::Int, tcol::Int, bval::Float64,
)
    return bval * (zone.delta[fcol] - zone.delta[tcol])
end

# ------------------------------------------------------------
# Merit function Φ(α)
# ------------------------------------------------------------

function _dc_merit(ws::DCStateWorkspaceALADIN, α::Float64)
    Φ = 0.0
    for zone in ws.zones
        Φ += _dc_cost_at_trial(zone, α)
    end
    r_sq = 0.0
    λr = 0.0
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        xi = zi.y[pair.col_i] + α * zi.delta[pair.col_i]
        xj = zj.y[pair.col_j] + α * zj.delta[pair.col_j]
        d = xi - xj
        r_sq += d * d
        λr += ws.vars.theta.lambda[k] * d
    end
    return Φ + λr + 0.5 * ws.mu[] * r_sq
end

function _dc_merit(ws::DCPowerWorkspaceALADIN, α::Float64)
    Φ = 0.0
    for zone in ws.zones
        Φ += _dc_cost_at_trial(zone, α)
    end
    r_sq = 0.0
    λr = 0.0
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = _dc_pflow_at(zi, pair.line_local_i, α)
        pj_val = _dc_pflow_at(zj, pair.line_local_j, α)
        d = pi_val - pj_val
        r_sq += d * d
        λr += ws.vars.p_flow.lambda[k] * d
    end
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = _dc_dummy_pflow_at(
            zi, pair.from_col_i, pair.to_col_i, pair.b_i, α,
        )
        pj_val = _dc_dummy_pflow_at(
            zj, pair.from_col_j, pair.to_col_j, pair.b_j, α,
        )
        d = pi_val - pj_val
        r_sq += d * d
        λr += ws.vars.p_dummy.lambda[k] * d
    end
    return Φ + λr + 0.5 * ws.mu[] * r_sq
end

function _dc_merit(ws::DCHybridWorkspaceALADIN, α::Float64)
    Φ = 0.0
    for zone in ws.zones
        Φ += _dc_cost_at_trial(zone, α)
    end
    r_sq = 0.0
    λr = 0.0
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        xi = zi.y[pair.col_i] + α * zi.delta[pair.col_i]
        xj = zj.y[pair.col_j] + α * zj.delta[pair.col_j]
        d = xi - xj
        r_sq += d * d
        λr += ws.vars.theta.lambda[k] * d
    end
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = _dc_pflow_at(zi, pair.line_local_i, α)
        pj_val = _dc_pflow_at(zj, pair.line_local_j, α)
        d = pi_val - pj_val
        r_sq += d * d
        λr += ws.vars.p_flow.lambda[k] * d
    end
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = _dc_dummy_pflow_at(
            zi, pair.from_col_i, pair.to_col_i, pair.b_i, α,
        )
        pj_val = _dc_dummy_pflow_at(
            zj, pair.from_col_j, pair.to_col_j, pair.b_j, α,
        )
        d = pi_val - pj_val
        r_sq += d * d
        λr += ws.vars.p_dummy.lambda[k] * d
    end
    @inbounds for (k, pair) in enumerate(ws.split_bus_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        # split_bus consensus is x_i + x_j = 0 (sign-flipped sum convention).
        xi = zi.y[pair.pinj_col_i] + α * zi.delta[pair.pinj_col_i]
        xj = zj.y[pair.pinj_col_j] + α * zj.delta[pair.pinj_col_j]
        d = xi + xj
        r_sq += d * d
        λr += ws.vars.p_inj.lambda[k] * d
    end
    return Φ + λr + 0.5 * ws.mu[] * r_sq
end

# ------------------------------------------------------------
# Analytical directional derivative dΦ/dα at α = 0
# ------------------------------------------------------------

function _dc_dΦ0(ws::DCStateWorkspaceALADIN)
    dΦ = 0.0
    for zone in ws.zones
        dΦ += _dc_cost_dir_deriv_0(zone)
    end
    mu = ws.mu[]
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        r0 = zi.y[pair.col_i] - zj.y[pair.col_j]
        Δr = zi.delta[pair.col_i] - zj.delta[pair.col_j]
        dΦ += (ws.vars.theta.lambda[k] + mu * r0) * Δr
    end
    return dΦ
end

function _dc_dΦ0(ws::DCPowerWorkspaceALADIN)
    dΦ = 0.0
    for zone in ws.zones
        dΦ += _dc_cost_dir_deriv_0(zone)
    end
    mu = ws.mu[]
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        r0 =
            _dc_pflow_at(zi, pair.line_local_i, 0.0) -
            _dc_pflow_at(zj, pair.line_local_j, 0.0)
        Δr =
            _dc_pflow_dir_deriv(zi, pair.line_local_i) -
            _dc_pflow_dir_deriv(zj, pair.line_local_j)
        dΦ += (ws.vars.p_flow.lambda[k] + mu * r0) * Δr
    end
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        r0 =
            _dc_dummy_pflow_at(zi, pair.from_col_i, pair.to_col_i, pair.b_i, 0.0) -
            _dc_dummy_pflow_at(zj, pair.from_col_j, pair.to_col_j, pair.b_j, 0.0)
        Δr =
            _dc_dummy_pflow_dir_deriv(zi, pair.from_col_i, pair.to_col_i, pair.b_i) -
            _dc_dummy_pflow_dir_deriv(zj, pair.from_col_j, pair.to_col_j, pair.b_j)
        dΦ += (ws.vars.p_dummy.lambda[k] + mu * r0) * Δr
    end
    return dΦ
end

function _dc_dΦ0(ws::DCHybridWorkspaceALADIN)
    dΦ = 0.0
    for zone in ws.zones
        dΦ += _dc_cost_dir_deriv_0(zone)
    end
    mu = ws.mu[]
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        r0 = zi.y[pair.col_i] - zj.y[pair.col_j]
        Δr = zi.delta[pair.col_i] - zj.delta[pair.col_j]
        dΦ += (ws.vars.theta.lambda[k] + mu * r0) * Δr
    end
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        r0 =
            _dc_pflow_at(zi, pair.line_local_i, 0.0) -
            _dc_pflow_at(zj, pair.line_local_j, 0.0)
        Δr =
            _dc_pflow_dir_deriv(zi, pair.line_local_i) -
            _dc_pflow_dir_deriv(zj, pair.line_local_j)
        dΦ += (ws.vars.p_flow.lambda[k] + mu * r0) * Δr
    end
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        r0 =
            _dc_dummy_pflow_at(zi, pair.from_col_i, pair.to_col_i, pair.b_i, 0.0) -
            _dc_dummy_pflow_at(zj, pair.from_col_j, pair.to_col_j, pair.b_j, 0.0)
        Δr =
            _dc_dummy_pflow_dir_deriv(zi, pair.from_col_i, pair.to_col_i, pair.b_i) -
            _dc_dummy_pflow_dir_deriv(zj, pair.from_col_j, pair.to_col_j, pair.b_j)
        dΦ += (ws.vars.p_dummy.lambda[k] + mu * r0) * Δr
    end
    @inbounds for (k, pair) in enumerate(ws.split_bus_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        r0 = zi.y[pair.pinj_col_i] + zj.y[pair.pinj_col_j]
        Δr = zi.delta[pair.pinj_col_i] + zj.delta[pair.pinj_col_j]
        dΦ += (ws.vars.p_inj.lambda[k] + mu * r0) * Δr
    end
    return dΦ
end

# ------------------------------------------------------------
# Armijo backtracking
# ------------------------------------------------------------

function _line_search_step!(
    ws::_DCPhase2Workspace;
    c1::Float64 = 1e-4,
    max_backtrack::Int = 25,
    α_min::Float64 = 1e-3,
)
    ws.line_search || return nothing

    Φ_0 = _dc_merit(ws, 0.0)
    dΦ_0 = _dc_dΦ0(ws)

    if !isfinite(dΦ_0) || dΦ_0 >= 0.0
        ws.alpha2 = α_min
        return nothing
    end

    α = ws.alpha2_init
    @inbounds for _ in 1:max_backtrack
        Φ_α = _dc_merit(ws, α)
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
