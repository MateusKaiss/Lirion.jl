# ============================================================
# ALADIN-SOCWR — Phase 2 line search (Armijo backtracking)
# ============================================================

const _SOCWRPhase2Workspace = Union{
    SOCWRStateWorkspaceALADIN,
    SOCWRPowerWorkspaceALADIN,
    SOCWRHybridWorkspaceALADIN,
}

# ============================================================
# Per-zone cost trial evaluations
# ============================================================

@inline function _socwr_cost_at_trial(zone::AladinZoneDataSOCWR, alpha::Float64)
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

@inline function _socwr_cost_dir_deriv_0(zone::AladinZoneDataSOCWR)
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
# SOCWR flow trial evaluations
# ============================================================

@inline function _socwr_flow_at_trial(
    zone::AladinZoneDataSOCWR,
    le::Int,
    alpha::Float64,
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
        acc += coefs[i] * (y[c] + alpha * delta[c])
    end
    return acc
end

@inline function _socwr_flow_dir_deriv(
    zone::AladinZoneDataSOCWR,
    le::Int,
    which::Symbol,
)::Float64
    cols = zone.edge_flow_cols[le]
    coefs = which === :p ? zone.edge_pflow_coefs[le] : zone.edge_qflow_coefs[le]
    acc = 0.0
    delta = zone.delta
    @inbounds for i in eachindex(cols)
        acc += coefs[i] * delta[cols[i]]
    end
    return acc
end

# ============================================================
# Merit function phi(alpha)
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function _socwr_merit(ws::SOCWRStateWorkspaceALADIN, alpha::Float64)
    phi = 0.0
    for zone in ws.zones
        phi += _socwr_cost_at_trial(zone, alpha)
    end

    r_sq = 0.0
    lambda_r = 0.0
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        xi = zi.y[pair.v2_col_i] + alpha * zi.delta[pair.v2_col_i]
        xj = zj.y[pair.v2_col_j] + alpha * zj.delta[pair.v2_col_j]
        d = xi - xj
        r_sq += d * d
        lambda_r += ws.vars.v2.lambda[k] * d
    end
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        ri_col = zi.edge_r_cols[pair.line_local_i]
        rj_col = zj.edge_r_cols[pair.line_local_j]
        ti_col = zi.edge_t_cols[pair.line_local_i]
        tj_col = zj.edge_t_cols[pair.line_local_j]
        dr = (zi.y[ri_col] + alpha * zi.delta[ri_col]) -
             (zj.y[rj_col] + alpha * zj.delta[rj_col])
        dt = (zi.y[ti_col] + alpha * zi.delta[ti_col]) -
             (zj.y[tj_col] + alpha * zj.delta[tj_col])
        r_sq += dr * dr + dt * dt
        lambda_r += ws.vars.r.lambda[k] * dr + ws.vars.t.lambda[k] * dt
    end

    mu = ws.mu[]
    return phi + lambda_r + 0.5 * mu * r_sq
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function _socwr_merit(ws::SOCWRPowerWorkspaceALADIN, alpha::Float64)
    phi = 0.0
    for zone in ws.zones
        phi += _socwr_cost_at_trial(zone, alpha)
    end

    r_sq = 0.0
    lambda_r = 0.0
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = _socwr_flow_at_trial(zi, pair.line_local_i, alpha, :p)
        pj_val = _socwr_flow_at_trial(zj, pair.line_local_j, alpha, :p)
        qi_val = _socwr_flow_at_trial(zi, pair.line_local_i, alpha, :q)
        qj_val = _socwr_flow_at_trial(zj, pair.line_local_j, alpha, :q)
        dp = pi_val - pj_val
        dq = qi_val - qj_val
        r_sq += dp * dp + dq * dq
        lambda_r += ws.vars.p_flow.lambda[k] * dp + ws.vars.q_flow.lambda[k] * dq
    end

    mu = ws.mu[]
    return phi + lambda_r + 0.5 * mu * r_sq
end

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

function _socwr_merit(ws::SOCWRHybridWorkspaceALADIN, alpha::Float64)
    phi = 0.0
    for zone in ws.zones
        phi += _socwr_cost_at_trial(zone, alpha)
    end

    r_sq = 0.0
    lambda_r = 0.0
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        xi = zi.y[pair.v2_col_i] + alpha * zi.delta[pair.v2_col_i]
        xj = zj.y[pair.v2_col_j] + alpha * zj.delta[pair.v2_col_j]
        d = xi - xj
        r_sq += d * d
        lambda_r += ws.vars.v2.lambda[k] * d
    end
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        ri_col = zi.edge_r_cols[pair.line_local_i]
        rj_col = zj.edge_r_cols[pair.line_local_j]
        ti_col = zi.edge_t_cols[pair.line_local_i]
        tj_col = zj.edge_t_cols[pair.line_local_j]
        dr = (zi.y[ri_col] + alpha * zi.delta[ri_col]) -
             (zj.y[rj_col] + alpha * zj.delta[rj_col])
        dt = (zi.y[ti_col] + alpha * zi.delta[ti_col]) -
             (zj.y[tj_col] + alpha * zj.delta[tj_col])

        pi_val = _socwr_flow_at_trial(zi, pair.line_local_i, alpha, :p)
        pj_val = _socwr_flow_at_trial(zj, pair.line_local_j, alpha, :p)
        qi_val = _socwr_flow_at_trial(zi, pair.line_local_i, alpha, :q)
        qj_val = _socwr_flow_at_trial(zj, pair.line_local_j, alpha, :q)
        dp = pi_val - pj_val
        dq = qi_val - qj_val

        r_sq += dr * dr + dt * dt + dp * dp + dq * dq
        lambda_r +=
            ws.vars.r.lambda[k] * dr + ws.vars.t.lambda[k] * dt +
            ws.vars.p_flow.lambda[k] * dp + ws.vars.q_flow.lambda[k] * dq
    end

    mu = ws.mu[]
    return phi + lambda_r + 0.5 * mu * r_sq
end

# ============================================================
# Directional derivative dphi/dalpha at alpha = 0
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function _socwr_dphi_0(ws::SOCWRStateWorkspaceALADIN)
    dphi = 0.0
    for zone in ws.zones
        dphi += _socwr_cost_dir_deriv_0(zone)
    end
    mu = ws.mu[]

    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        r0 = zi.y[pair.v2_col_i] - zj.y[pair.v2_col_j]
        dr = zi.delta[pair.v2_col_i] - zj.delta[pair.v2_col_j]
        dphi += (ws.vars.v2.lambda[k] + mu * r0) * dr
    end
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        ri_col = zi.edge_r_cols[pair.line_local_i]
        rj_col = zj.edge_r_cols[pair.line_local_j]
        ti_col = zi.edge_t_cols[pair.line_local_i]
        tj_col = zj.edge_t_cols[pair.line_local_j]

        r0_r = zi.y[ri_col] - zj.y[rj_col]
        dr_r = zi.delta[ri_col] - zj.delta[rj_col]
        r0_t = zi.y[ti_col] - zj.y[tj_col]
        dr_t = zi.delta[ti_col] - zj.delta[tj_col]

        dphi += (ws.vars.r.lambda[k] + mu * r0_r) * dr_r
        dphi += (ws.vars.t.lambda[k] + mu * r0_t) * dr_t
    end
    return dphi
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function _socwr_dphi_0(ws::SOCWRPowerWorkspaceALADIN)
    dphi = 0.0
    for zone in ws.zones
        dphi += _socwr_cost_dir_deriv_0(zone)
    end
    mu = ws.mu[]

    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]

        pi_val = _socwr_flow_at_trial(zi, pair.line_local_i, 0.0, :p)
        pj_val = _socwr_flow_at_trial(zj, pair.line_local_j, 0.0, :p)
        qi_val = _socwr_flow_at_trial(zi, pair.line_local_i, 0.0, :q)
        qj_val = _socwr_flow_at_trial(zj, pair.line_local_j, 0.0, :q)
        dpi = _socwr_flow_dir_deriv(zi, pair.line_local_i, :p)
        dpj = _socwr_flow_dir_deriv(zj, pair.line_local_j, :p)
        dqi = _socwr_flow_dir_deriv(zi, pair.line_local_i, :q)
        dqj = _socwr_flow_dir_deriv(zj, pair.line_local_j, :q)

        r0_p = pi_val - pj_val
        dr_p = dpi - dpj
        r0_q = qi_val - qj_val
        dr_q = dqi - dqj

        dphi += (ws.vars.p_flow.lambda[k] + mu * r0_p) * dr_p
        dphi += (ws.vars.q_flow.lambda[k] + mu * r0_q) * dr_q
    end
    return dphi
end

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

function _socwr_dphi_0(ws::SOCWRHybridWorkspaceALADIN)
    dphi = 0.0
    for zone in ws.zones
        dphi += _socwr_cost_dir_deriv_0(zone)
    end
    mu = ws.mu[]

    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        r0 = zi.y[pair.v2_col_i] - zj.y[pair.v2_col_j]
        dr = zi.delta[pair.v2_col_i] - zj.delta[pair.v2_col_j]
        dphi += (ws.vars.v2.lambda[k] + mu * r0) * dr
    end

    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        ri_col = zi.edge_r_cols[pair.line_local_i]
        rj_col = zj.edge_r_cols[pair.line_local_j]
        ti_col = zi.edge_t_cols[pair.line_local_i]
        tj_col = zj.edge_t_cols[pair.line_local_j]

        r0_r = zi.y[ri_col] - zj.y[rj_col]
        dr_r = zi.delta[ri_col] - zj.delta[rj_col]
        r0_t = zi.y[ti_col] - zj.y[tj_col]
        dr_t = zi.delta[ti_col] - zj.delta[tj_col]

        pi_val = _socwr_flow_at_trial(zi, pair.line_local_i, 0.0, :p)
        pj_val = _socwr_flow_at_trial(zj, pair.line_local_j, 0.0, :p)
        qi_val = _socwr_flow_at_trial(zi, pair.line_local_i, 0.0, :q)
        qj_val = _socwr_flow_at_trial(zj, pair.line_local_j, 0.0, :q)
        dpi = _socwr_flow_dir_deriv(zi, pair.line_local_i, :p)
        dpj = _socwr_flow_dir_deriv(zj, pair.line_local_j, :p)
        dqi = _socwr_flow_dir_deriv(zi, pair.line_local_i, :q)
        dqj = _socwr_flow_dir_deriv(zj, pair.line_local_j, :q)

        r0_p = pi_val - pj_val
        dr_p = dpi - dpj
        r0_q = qi_val - qj_val
        dr_q = dqi - dqj

        dphi += (ws.vars.r.lambda[k] + mu * r0_r) * dr_r
        dphi += (ws.vars.t.lambda[k] + mu * r0_t) * dr_t
        dphi += (ws.vars.p_flow.lambda[k] + mu * r0_p) * dr_p
        dphi += (ws.vars.q_flow.lambda[k] + mu * r0_q) * dr_q
    end
    return dphi
end

# ============================================================
# Armijo backtracking step
# ============================================================

function _line_search_step!(
    ws::_SOCWRPhase2Workspace;
    c1::Float64 = 1e-4,
    max_backtrack::Int = 25,
    alpha_min::Float64 = 1e-3,
)
    ws.line_search || return nothing

    phi_0 = _socwr_merit(ws, 0.0)
    dphi_0 = _socwr_dphi_0(ws)

    if !isfinite(dphi_0) || dphi_0 >= 0.0
        ws.alpha2 = alpha_min
        return nothing
    end

    alpha = ws.alpha2_init
    @inbounds for _ in 1:max_backtrack
        phi_alpha = _socwr_merit(ws, alpha)
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
