# ============================================================
# ALADIN-DC — Phase 2 line search (Armijo backtracking)
# ============================================================

const _DCPhase2Workspace = Union{
    DCStateWorkspaceALADIN,
    DCPowerWorkspaceALADIN,
    DCHybridWorkspaceALADIN,
}

# ============================================================
# Per-zone cost trial evaluations
# ============================================================

@inline function _dc_cost_at_trial(zone::AladinZoneData, alpha::Float64)
    acc = 0.0
    y = zone.y
    delta = zone.delta
    @inbounds for k in eachindex(zone.pg_cols)
        col = zone.pg_cols[k]
        xk = y[col] + alpha * delta[col]
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

# ============================================================
# Interzonal line flow trial evaluations
# ============================================================

@inline function _dc_pflow_at(zone::AladinZoneData, le::Int, alpha::Float64)
    fcol = zone.edge_from_cols[le]
    tcol = zone.edge_to_cols[le]
    bval = zone.edge_b[le]
    yf = zone.y[fcol] + alpha * zone.delta[fcol]
    yt = zone.y[tcol] + alpha * zone.delta[tcol]
    return bval * (yf - yt)
end

@inline function _dc_pflow_dir_deriv(zone::AladinZoneData, le::Int)
    fcol = zone.edge_from_cols[le]
    tcol = zone.edge_to_cols[le]
    bval = zone.edge_b[le]
    return bval * (zone.delta[fcol] - zone.delta[tcol])
end

# ============================================================
# Dummy line flow trial evaluations
# ============================================================

@inline function _dc_dummy_pflow_at(
    zone::AladinZoneData,
    fcol::Int, tcol::Int, bval::Float64,
    alpha::Float64,
)
    yf = zone.y[fcol] + alpha * zone.delta[fcol]
    yt = zone.y[tcol] + alpha * zone.delta[tcol]
    return bval * (yf - yt)
end

@inline function _dc_dummy_pflow_dir_deriv(
    zone::AladinZoneData,
    fcol::Int, tcol::Int, bval::Float64,
)
    return bval * (zone.delta[fcol] - zone.delta[tcol])
end

# ============================================================
# Merit function phi(alpha)
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function _dc_merit(ws::DCStateWorkspaceALADIN, alpha::Float64)
    phi = 0.0
    for zone in ws.zones
        phi += _dc_cost_at_trial(zone, alpha)
    end
    r_sq = 0.0
    lambda_r = 0.0
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        xi = zi.y[pair.col_i] + alpha * zi.delta[pair.col_i]
        xj = zj.y[pair.col_j] + alpha * zj.delta[pair.col_j]
        d = xi - xj
        r_sq += d * d
        lambda_r += ws.vars.theta.lambda[k] * d
    end
    return phi + lambda_r + 0.5 * ws.mu[] * r_sq
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function _dc_merit(ws::DCPowerWorkspaceALADIN, alpha::Float64)
    phi = 0.0
    for zone in ws.zones
        phi += _dc_cost_at_trial(zone, alpha)
    end
    r_sq = 0.0
    lambda_r = 0.0
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = _dc_pflow_at(zi, pair.line_local_i, alpha)
        pj_val = _dc_pflow_at(zj, pair.line_local_j, alpha)
        d = pi_val - pj_val
        r_sq += d * d
        lambda_r += ws.vars.p_flow.lambda[k] * d
    end
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = _dc_dummy_pflow_at(
            zi, pair.from_col_i, pair.to_col_i, pair.b_i, alpha,
        )
        pj_val = _dc_dummy_pflow_at(
            zj, pair.from_col_j, pair.to_col_j, pair.b_j, alpha,
        )
        d = pi_val - pj_val
        r_sq += d * d
        lambda_r += ws.vars.p_dummy.lambda[k] * d
    end
    return phi + lambda_r + 0.5 * ws.mu[] * r_sq
end

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

function _dc_merit(ws::DCHybridWorkspaceALADIN, alpha::Float64)
    phi = 0.0
    for zone in ws.zones
        phi += _dc_cost_at_trial(zone, alpha)
    end
    r_sq = 0.0
    lambda_r = 0.0
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        xi = zi.y[pair.col_i] + alpha * zi.delta[pair.col_i]
        xj = zj.y[pair.col_j] + alpha * zj.delta[pair.col_j]
        d = xi - xj
        r_sq += d * d
        lambda_r += ws.vars.theta.lambda[k] * d
    end
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = _dc_pflow_at(zi, pair.line_local_i, alpha)
        pj_val = _dc_pflow_at(zj, pair.line_local_j, alpha)
        d = pi_val - pj_val
        r_sq += d * d
        lambda_r += ws.vars.p_flow.lambda[k] * d
    end
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = _dc_dummy_pflow_at(
            zi, pair.from_col_i, pair.to_col_i, pair.b_i, alpha,
        )
        pj_val = _dc_dummy_pflow_at(
            zj, pair.from_col_j, pair.to_col_j, pair.b_j, alpha,
        )
        d = pi_val - pj_val
        r_sq += d * d
        lambda_r += ws.vars.p_dummy.lambda[k] * d
    end
    @inbounds for (k, pair) in enumerate(ws.split_bus_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        xi = zi.y[pair.p_inj_col_i] + alpha * zi.delta[pair.p_inj_col_i]
        xj = zj.y[pair.p_inj_col_j] + alpha * zj.delta[pair.p_inj_col_j]
        d = xi + xj
        r_sq += d * d
        lambda_r += ws.vars.p_inj.lambda[k] * d
    end
    return phi + lambda_r + 0.5 * ws.mu[] * r_sq
end

# ============================================================
# Directional derivative dphi/dalpha at alpha = 0
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function _dc_dphi_0(ws::DCStateWorkspaceALADIN)
    dphi = 0.0
    for zone in ws.zones
        dphi += _dc_cost_dir_deriv_0(zone)
    end
    mu = ws.mu[]
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        r0 = zi.y[pair.col_i] - zj.y[pair.col_j]
        dr = zi.delta[pair.col_i] - zj.delta[pair.col_j]
        dphi += (ws.vars.theta.lambda[k] + mu * r0) * dr
    end
    return dphi
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function _dc_dphi_0(ws::DCPowerWorkspaceALADIN)
    dphi = 0.0
    for zone in ws.zones
        dphi += _dc_cost_dir_deriv_0(zone)
    end
    mu = ws.mu[]
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        r0 =
            _dc_pflow_at(zi, pair.line_local_i, 0.0) -
            _dc_pflow_at(zj, pair.line_local_j, 0.0)
        dr =
            _dc_pflow_dir_deriv(zi, pair.line_local_i) -
            _dc_pflow_dir_deriv(zj, pair.line_local_j)
        dphi += (ws.vars.p_flow.lambda[k] + mu * r0) * dr
    end
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        r0 =
            _dc_dummy_pflow_at(zi, pair.from_col_i, pair.to_col_i, pair.b_i, 0.0) -
            _dc_dummy_pflow_at(zj, pair.from_col_j, pair.to_col_j, pair.b_j, 0.0)
        dr =
            _dc_dummy_pflow_dir_deriv(zi, pair.from_col_i, pair.to_col_i, pair.b_i) -
            _dc_dummy_pflow_dir_deriv(zj, pair.from_col_j, pair.to_col_j, pair.b_j)
        dphi += (ws.vars.p_dummy.lambda[k] + mu * r0) * dr
    end
    return dphi
end

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

function _dc_dphi_0(ws::DCHybridWorkspaceALADIN)
    dphi = 0.0
    for zone in ws.zones
        dphi += _dc_cost_dir_deriv_0(zone)
    end
    mu = ws.mu[]
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        r0 = zi.y[pair.col_i] - zj.y[pair.col_j]
        dr = zi.delta[pair.col_i] - zj.delta[pair.col_j]
        dphi += (ws.vars.theta.lambda[k] + mu * r0) * dr
    end
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        r0 =
            _dc_pflow_at(zi, pair.line_local_i, 0.0) -
            _dc_pflow_at(zj, pair.line_local_j, 0.0)
        dr =
            _dc_pflow_dir_deriv(zi, pair.line_local_i) -
            _dc_pflow_dir_deriv(zj, pair.line_local_j)
        dphi += (ws.vars.p_flow.lambda[k] + mu * r0) * dr
    end
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        r0 =
            _dc_dummy_pflow_at(zi, pair.from_col_i, pair.to_col_i, pair.b_i, 0.0) -
            _dc_dummy_pflow_at(zj, pair.from_col_j, pair.to_col_j, pair.b_j, 0.0)
        dr =
            _dc_dummy_pflow_dir_deriv(zi, pair.from_col_i, pair.to_col_i, pair.b_i) -
            _dc_dummy_pflow_dir_deriv(zj, pair.from_col_j, pair.to_col_j, pair.b_j)
        dphi += (ws.vars.p_dummy.lambda[k] + mu * r0) * dr
    end
    @inbounds for (k, pair) in enumerate(ws.split_bus_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        r0 = zi.y[pair.p_inj_col_i] + zj.y[pair.p_inj_col_j]
        dr = zi.delta[pair.p_inj_col_i] + zj.delta[pair.p_inj_col_j]
        dphi += (ws.vars.p_inj.lambda[k] + mu * r0) * dr
    end
    return dphi
end

# ============================================================
# Armijo backtracking step
# ============================================================

function _line_search_step!(
    ws::_DCPhase2Workspace;
    c1::Float64 = 1e-4,
    max_backtrack::Int = 25,
    alpha_min::Float64 = 1e-3,
)
    ws.line_search || return nothing

    phi_0 = _dc_merit(ws, 0.0)
    dphi_0 = _dc_dphi_0(ws)

    if !isfinite(dphi_0) || dphi_0 >= 0.0
        ws.alpha2 = alpha_min
        return nothing
    end

    alpha = ws.alpha2_init
    @inbounds for _ in 1:max_backtrack
        phi_alpha = _dc_merit(ws, alpha)
        if isfinite(phi_alpha) && phi_alpha <= phi_0 + c1 * alpha * dphi_0
            ws.alpha2 = alpha
            return nothing
        end
        alpha *= 0.5
        if alpha < alpha_min
            break
        end
    end
    ws.alpha2 = alpha_min
    return nothing
end
