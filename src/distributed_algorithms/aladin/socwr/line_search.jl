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
# Dummy-half flow trial evaluations (DummyBus sectioning)
# ============================================================

@inline function _socwr_dummy_side(pair::DummyLinePairSOCWR, side::Symbol, which::Symbol)
    if side === :i
        cols = pair.cols_i
        coefs = which === :p ? pair.pcoef_i : pair.qcoef_i
        cst = which === :p ? pair.pconst_i : pair.qconst_i
    else
        cols = pair.cols_j
        coefs = which === :p ? pair.pcoef_j : pair.qcoef_j
        cst = which === :p ? pair.pconst_j : pair.qconst_j
    end
    return cols, coefs, cst
end

@inline function _socwr_dummy_flow_at_trial(
    pair::DummyLinePairSOCWR,
    zone::AladinZoneDataSOCWR,
    side::Symbol,
    alpha::Float64,
    which::Symbol,
)::Float64
    cols, coefs, cst = _socwr_dummy_side(pair, side, which)
    acc = cst
    y = zone.y
    delta = zone.delta
    @inbounds for i in eachindex(cols)
        c = cols[i]
        acc += coefs[i] * (y[c] + alpha * delta[c])
    end
    return acc
end

@inline function _socwr_dummy_flow_dir_deriv(
    pair::DummyLinePairSOCWR,
    zone::AladinZoneDataSOCWR,
    side::Symbol,
    which::Symbol,
)::Float64
    cols, coefs, _ = _socwr_dummy_side(pair, side, which)
    acc = 0.0
    delta = zone.delta
    @inbounds for i in eachindex(cols)
        acc += coefs[i] * delta[cols[i]]
    end
    return acc
end

@inline function _socwr_dummy_flow_pt(
    pair::DummyLinePairSOCWR,
    zone::AladinZoneDataSOCWR,
    side::Symbol,
    which::Symbol,
    acc::F,
)::Float64 where {F}
    cols, coefs, cst = _socwr_dummy_side(pair, side, which)
    val = cst
    @inbounds for i in eachindex(cols)
        val += coefs[i] * acc(zone, cols[i])
    end
    return val
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
        dr =
            (zi.y[ri_col] + alpha * zi.delta[ri_col]) -
            (zj.y[rj_col] + alpha * zj.delta[rj_col])
        dt =
            (zi.y[ti_col] + alpha * zi.delta[ti_col]) -
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
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = _socwr_dummy_flow_at_trial(pair, zi, :i, alpha, :p)
        pj_val = _socwr_dummy_flow_at_trial(pair, zj, :j, alpha, :p)
        qi_val = _socwr_dummy_flow_at_trial(pair, zi, :i, alpha, :q)
        qj_val = _socwr_dummy_flow_at_trial(pair, zj, :j, alpha, :q)
        dp = pi_val - pj_val
        dq = qi_val - qj_val
        r_sq += dp * dp + dq * dq
        lambda_r += ws.vars.p_dummy.lambda[k] * dp + ws.vars.q_dummy.lambda[k] * dq
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
        dr =
            (zi.y[ri_col] + alpha * zi.delta[ri_col]) -
            (zj.y[rj_col] + alpha * zj.delta[rj_col])
        dt =
            (zi.y[ti_col] + alpha * zi.delta[ti_col]) -
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
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = _socwr_dummy_flow_at_trial(pair, zi, :i, alpha, :p)
        pj_val = _socwr_dummy_flow_at_trial(pair, zj, :j, alpha, :p)
        qi_val = _socwr_dummy_flow_at_trial(pair, zi, :i, alpha, :q)
        qj_val = _socwr_dummy_flow_at_trial(pair, zj, :j, alpha, :q)
        dp = pi_val - pj_val
        dq = qi_val - qj_val
        r_sq += dp * dp + dq * dq
        lambda_r += ws.vars.p_dummy.lambda[k] * dp + ws.vars.q_dummy.lambda[k] * dq
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
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]

        pi_val = _socwr_dummy_flow_at_trial(pair, zi, :i, 0.0, :p)
        pj_val = _socwr_dummy_flow_at_trial(pair, zj, :j, 0.0, :p)
        qi_val = _socwr_dummy_flow_at_trial(pair, zi, :i, 0.0, :q)
        qj_val = _socwr_dummy_flow_at_trial(pair, zj, :j, 0.0, :q)
        dpi = _socwr_dummy_flow_dir_deriv(pair, zi, :i, :p)
        dpj = _socwr_dummy_flow_dir_deriv(pair, zj, :j, :p)
        dqi = _socwr_dummy_flow_dir_deriv(pair, zi, :i, :q)
        dqj = _socwr_dummy_flow_dir_deriv(pair, zj, :j, :q)

        r0_p = pi_val - pj_val
        dr_p = dpi - dpj
        r0_q = qi_val - qj_val
        dr_q = dqi - dqj

        dphi += (ws.vars.p_dummy.lambda[k] + mu * r0_p) * dr_p
        dphi += (ws.vars.q_dummy.lambda[k] + mu * r0_q) * dr_q
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
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]

        pi_val = _socwr_dummy_flow_at_trial(pair, zi, :i, 0.0, :p)
        pj_val = _socwr_dummy_flow_at_trial(pair, zj, :j, 0.0, :p)
        qi_val = _socwr_dummy_flow_at_trial(pair, zi, :i, 0.0, :q)
        qj_val = _socwr_dummy_flow_at_trial(pair, zj, :j, 0.0, :q)
        dpi = _socwr_dummy_flow_dir_deriv(pair, zi, :i, :p)
        dpj = _socwr_dummy_flow_dir_deriv(pair, zj, :j, :p)
        dqi = _socwr_dummy_flow_dir_deriv(pair, zi, :i, :q)
        dqj = _socwr_dummy_flow_dir_deriv(pair, zj, :j, :q)

        r0_p = pi_val - pj_val
        dr_p = dpi - dpj
        r0_q = qi_val - qj_val
        dr_q = dqi - dqj

        dphi += (ws.vars.p_dummy.lambda[k] + mu * r0_p) * dr_p
        dphi += (ws.vars.q_dummy.lambda[k] + mu * r0_q) * dr_q
    end
    return dphi
end

# ============================================================
# HFD Algorithm 3 globalization — L1 exact penalty Φ, steps (a)/(b)/(c)
# ============================================================

@inline function _socwr_cost_pt(zone::AladinZoneDataSOCWR, acc::F) where {F}
    total = 0.0
    @inbounds for k in eachindex(zone.pg_cols)
        col = zone.pg_cols[k]
        xk = acc(zone, col)
        total += 0.5 * zone.pg_c2_sbase2[k] * xk * xk + zone.pg_c1_sbase[k] * xk
    end
    return total
end

_socwr_cost_sum(ws::_SOCWRPhase2Workspace, acc::F) where {F} =
    sum(_socwr_cost_pt(zone, acc) for zone in ws.zones)

@inline function _socwr_flow_pt(
    zone::AladinZoneDataSOCWR, le::Int, which::Symbol, acc::F,
)::Float64 where {F}
    cols = zone.edge_flow_cols[le]
    coefs = which === :p ? zone.edge_pflow_coefs[le] : zone.edge_qflow_coefs[le]
    acc_val = which === :p ? zone.edge_pflow_const[le] : zone.edge_qflow_const[le]
    @inbounds for i in eachindex(cols)
        acc_val += coefs[i] * acc(zone, cols[i])
    end
    return acc_val
end

# ℓ1 consensus residual ‖Ax‖₁ at the point defined by `acc`.
function _socwr_resid_l1(ws::SOCWRStateWorkspaceALADIN, acc::F) where {F}
    r = 0.0
    @inbounds for pair in ws.bus_pairs
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r += abs(acc(zi, pair.v2_col_i) - acc(zj, pair.v2_col_j))
    end
    @inbounds for pair in ws.interzonal_line_pairs
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r += abs(acc(zi, zi.edge_r_cols[pair.line_local_i]) -
                 acc(zj, zj.edge_r_cols[pair.line_local_j]))
        r += abs(acc(zi, zi.edge_t_cols[pair.line_local_i]) -
                 acc(zj, zj.edge_t_cols[pair.line_local_j]))
    end
    return r
end

function _socwr_resid_l1(ws::SOCWRPowerWorkspaceALADIN, acc::F) where {F}
    r = 0.0
    @inbounds for pair in ws.interzonal_line_pairs
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r += abs(_socwr_flow_pt(zi, pair.line_local_i, :p, acc) -
                 _socwr_flow_pt(zj, pair.line_local_j, :p, acc))
        r += abs(_socwr_flow_pt(zi, pair.line_local_i, :q, acc) -
                 _socwr_flow_pt(zj, pair.line_local_j, :q, acc))
    end
    @inbounds for pair in ws.dummy_line_pairs
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r += abs(_socwr_dummy_flow_pt(pair, zi, :i, :p, acc) -
                 _socwr_dummy_flow_pt(pair, zj, :j, :p, acc))
        r += abs(_socwr_dummy_flow_pt(pair, zi, :i, :q, acc) -
                 _socwr_dummy_flow_pt(pair, zj, :j, :q, acc))
    end
    return r
end

function _socwr_resid_l1(ws::SOCWRHybridWorkspaceALADIN, acc::F) where {F}
    r = 0.0
    @inbounds for pair in ws.bus_pairs
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r += abs(acc(zi, pair.v2_col_i) - acc(zj, pair.v2_col_j))
    end
    @inbounds for pair in ws.interzonal_line_pairs
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r += abs(acc(zi, zi.edge_r_cols[pair.line_local_i]) -
                 acc(zj, zj.edge_r_cols[pair.line_local_j]))
        r += abs(acc(zi, zi.edge_t_cols[pair.line_local_i]) -
                 acc(zj, zj.edge_t_cols[pair.line_local_j]))
        r += abs(_socwr_flow_pt(zi, pair.line_local_i, :p, acc) -
                 _socwr_flow_pt(zj, pair.line_local_j, :p, acc))
        r += abs(_socwr_flow_pt(zi, pair.line_local_i, :q, acc) -
                 _socwr_flow_pt(zj, pair.line_local_j, :q, acc))
    end
    @inbounds for pair in ws.split_bus_pairs
        dp = 0.0; dq = 0.0
        for c in eachindex(pair.zones)
            zone = ws.zones[pair.zones[c]]
            dp += acc(zone, pair.p_inj_cols[c])
            dq += acc(zone, pair.q_inj_cols[c])
        end
        r += abs(dp) + abs(dq)
    end
    @inbounds for pair in ws.dummy_line_pairs
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r += abs(_socwr_dummy_flow_pt(pair, zi, :i, :p, acc) -
                 _socwr_dummy_flow_pt(pair, zj, :j, :p, acc))
        r += abs(_socwr_dummy_flow_pt(pair, zi, :i, :q, acc) -
                 _socwr_dummy_flow_pt(pair, zj, :j, :q, acc))
    end
    return r
end

# Σ (ρ/2)‖yᵢ − xᵢ‖²_Σᵢ — value of the proximal quadratic at d = y − x.
function _socwr_prox(ws::_SOCWRPhase2Workspace)
    s = 0.0
    @inbounds for zone in ws.zones
        v2c = zone.var_to_col
        y = zone.y; x = zone.x_prev
        for (pair, a) in zone.quad_static.terms
            i = v2c[pair.a]; j = v2c[pair.b]
            s += a * (y[i] - x[i]) * (y[j] - x[j])
        end
    end
    return s
end

function _socwr_local_infeas_l1(ws::_SOCWRPhase2Workspace, acc::F) where {F}
    total = 0.0
    for zone in ws.zones
        model = zone.local_model
        v2c = zone.var_to_col
        valfn = v -> acc(zone, v2c[v])
        for (Fc, S) in JuMP.list_of_constraint_types(model)
            Fc <: JuMP.VariableRef && continue
            if S <: MOI.RotatedSecondOrderCone || S <: MOI.SecondOrderCone
                is_rsoc = S <: MOI.RotatedSecondOrderCone
                for ci in JuMP.all_constraints(model, Fc, S)
                    cobj = JuMP.constraint_object(ci)
                    fvals = try
                        Float64[JuMP.value(valfn, f) for f in cobj.func]
                    catch
                        return Inf
                    end
                    m = length(fvals)
                    m >= 2 || continue
                    if is_rsoc                      # 2 t u ≥ ‖x₃:‖²
                        q = -2.0 * fvals[1] * fvals[2]
                        @inbounds for k in 3:m
                            q += fvals[k] * fvals[k]
                        end
                    else                            # t ≥ ‖x₂:‖  (squared)
                        q = -fvals[1] * fvals[1]
                        @inbounds for k in 2:m
                            q += fvals[k] * fvals[k]
                        end
                    end
                    q > 0.0 && (total += q)
                end
                continue
            end
            for ci in JuMP.all_constraints(model, Fc, S)
                cobj = JuMP.constraint_object(ci)
                g = try
                    JuMP.value(valfn, cobj.func)
                catch
                    return Inf
                end
                s = cobj.set
                if s isa MOI.EqualTo{Float64}
                    total += abs(g - s.value)
                elseif s isa MOI.LessThan{Float64}
                    total += max(0.0, g - s.upper)
                elseif s isa MOI.GreaterThan{Float64}
                    total += max(0.0, s.lower - g)
                elseif s isa MOI.Interval{Float64}
                    total += max(0.0, g - s.upper) + max(0.0, s.lower - g)
                end
            end
        end
        lb = zone.lb; ub = zone.ub
        @inbounds for c in 1:zone.n_vars
            v = acc(zone, c)
            isnan(lb[c]) || (total += max(0.0, lb[c] - v))
            isnan(ub[c]) || (total += max(0.0, v - ub[c]))
        end
    end
    return total
end

function _socwr_max_local_dual(ws::_SOCWRPhase2Workspace)
    maxkap = 0.0
    for zone in ws.zones
        model = zone.local_model
        (try; JuMP.has_duals(model); catch; false; end) || continue
        for (Fc, S) in JuMP.list_of_constraint_types(model)
            Fc <: JuMP.VariableRef && continue
            for ci in JuMP.all_constraints(model, Fc, S)
                d = try
                    JuMP.dual(ci)
                catch
                    continue
                end
                if d isa Real
                    ad = abs(d); ad > maxkap && (maxkap = ad)
                elseif d isa AbstractVector
                    for e in d
                        ae = abs(e); ae > maxkap && (maxkap = ae)
                    end
                end
            end
        end
    end
    return maxkap
end

function _alg3_step!(ws::_SOCWRPhase2Workspace; gamma::Float64 = 1e-4)
    accx = (z, c) -> z.x_prev[c]
    accy = (z, c) -> z.y[c]
    accyd = (z, c) -> z.y[c] + z.delta[c]
    maxlam = 0.0
    for state in values(ws.vars)
        @inbounds for v in state.lambda
            av = abs(v); av > maxlam && (maxlam = av)
        end
    end
    use_viol = _ALG3_PHI_LOCAL_VIOL[]
    use_viol && (maxlam = max(maxlam, _socwr_max_local_dual(ws)))
    if ws.lam_bar <= 0.0 || ws.lam_bar < 1.5 * maxlam
        ws.lam_bar = max(10.0 * maxlam, 1.0)
    end
    lam_bar = ws.lam_bar

    infeas = use_viol ?
        (acc -> lam_bar * _socwr_local_infeas_l1(ws, acc)) : (acc -> 0.0)

    rhs = gamma * (_socwr_prox(ws) + lam_bar * _socwr_resid_l1(ws, accy))
    phi_x = _socwr_cost_sum(ws, accx) + lam_bar * _socwr_resid_l1(ws, accx) +
        infeas(accx)

    phi_a = _socwr_cost_sum(ws, accyd) + lam_bar * _socwr_resid_l1(ws, accyd) +
        infeas(accyd)
    if isfinite(phi_a) && phi_x - phi_a >= rhs
        ws.alpha1 = 1.0; ws.alpha2 = 1.0; ws.alpha3 = 1.0
        return nothing
    end

    phi_b = _socwr_cost_sum(ws, accy) + lam_bar * _socwr_resid_l1(ws, accy) +
        infeas(accy)
    if isfinite(phi_b) && phi_x - phi_b >= rhs
        ws.alpha1 = 1.0; ws.alpha2 = 0.0; ws.alpha3 = 0.0
        return nothing
    end

    return _branch_c_dual_linesearch!(ws)
end

# ============================================================
# HFD Algorithm 3 — :alg3_original (verbatim port, SOCWR)
# ============================================================

const _ALG3_ORIG_SAFETY = Ref(2.0)      # majorisation margin
const _ALG3_ORIG_FLOOR = Ref(1.0)       # λ̄, κ̄ ≥ floor (avoid 0 at λ = 0)

# Euclidean distance from v = (t, x₁..xₙ) to the SOC {t ≥ ‖x‖}.
@inline function _soc_dist(v::AbstractVector{Float64})::Float64
    t = v[1]
    n2 = 0.0
    @inbounds for k in 2:length(v)
        n2 += v[k] * v[k]
    end
    n = sqrt(n2)
    n <= t && return 0.0
    n <= -t && return sqrt(t * t + n2)
    return (n - t) / sqrt(2.0)
end

@inline function _rsoc_dist(v::AbstractVector{Float64})::Float64
    s2 = sqrt(2.0)
    n = length(v)
    w = Vector{Float64}(undef, n)          # (t', u', x₁..xₙ₋₂): same length as v
    w[1] = (v[1] + v[2]) / s2              # t'
    w[2] = (v[1] - v[2]) / s2              # u'  (first entry of the SOC "x" block)
    for k in 3:n
        w[k] = v[k]
    end
    return _soc_dist(w)
end

# Σ over cone constraints of dist(v(acc), K) — the exact conic violation.
function _socwr_cone_dist_l1(ws::_SOCWRPhase2Workspace, acc::F) where {F}
    total = 0.0
    for zone in ws.zones
        model = zone.local_model
        v2c = zone.var_to_col
        valfn = v -> acc(zone, v2c[v])
        for (Fc, S) in JuMP.list_of_constraint_types(model)
            Fc <: JuMP.VariableRef && continue
            is_rsoc = S <: MOI.RotatedSecondOrderCone
            is_soc = S <: MOI.SecondOrderCone
            (is_rsoc || is_soc) || continue
            for ci in JuMP.all_constraints(model, Fc, S)
                cobj = JuMP.constraint_object(ci)
                fvals = try
                    Float64[JuMP.value(valfn, f) for f in cobj.func]
                catch
                    return Inf
                end
                length(fvals) >= 2 || continue
                total += is_rsoc ? _rsoc_dist(fvals) : _soc_dist(fvals)
            end
        end
    end
    return total
end

function _socwr_scalar_infeas_l1(ws::_SOCWRPhase2Workspace, acc::F) where {F}
    total = 0.0
    for zone in ws.zones
        model = zone.local_model
        v2c = zone.var_to_col
        valfn = v -> acc(zone, v2c[v])
        for (Fc, S) in JuMP.list_of_constraint_types(model)
            Fc <: JuMP.VariableRef && continue
            (S <: MOI.RotatedSecondOrderCone || S <: MOI.SecondOrderCone) && continue
            for ci in JuMP.all_constraints(model, Fc, S)
                cobj = JuMP.constraint_object(ci)
                if cobj.func isa AbstractVector
                    vals = try
                        Float64[JuMP.value(valfn, f) for f in cobj.func]
                    catch
                        return Inf
                    end
                    if S <: MOI.Zeros
                        total += sum(abs, vals; init = 0.0)
                    elseif S <: MOI.Nonnegatives
                        total += sum(v -> max(0.0, -v), vals; init = 0.0)
                    elseif S <: MOI.Nonpositives
                        total += sum(v -> max(0.0, v), vals; init = 0.0)
                    else
                        error(":alg3_original Φ (SOCWR): unsupported vector " *
                              "constraint set $(S) in the local model.")
                    end
                    continue
                end
                g = try
                    JuMP.value(valfn, cobj.func)
                catch
                    return Inf
                end
                s = cobj.set
                if s isa MOI.EqualTo{Float64}
                    total += abs(g - s.value)
                elseif s isa MOI.LessThan{Float64}
                    total += max(0.0, g - s.upper)
                elseif s isa MOI.GreaterThan{Float64}
                    total += max(0.0, s.lower - g)
                elseif s isa MOI.Interval{Float64}
                    total += max(0.0, g - s.upper) + max(0.0, s.lower - g)
                end
            end
        end
        lb = zone.lb; ub = zone.ub
        @inbounds for c in 1:zone.n_vars
            v = acc(zone, c)
            isnan(lb[c]) || (total += max(0.0, lb[c] - v))
            isnan(ub[c]) || (total += max(0.0, v - ub[c]))
        end
    end
    return total
end

# Exact local penalty of Φ: scalar ℓ1 part + conic distance part.
_socwr_local_infeas_exact(ws::_SOCWRPhase2Workspace, acc::F) where {F} =
    _socwr_scalar_infeas_l1(ws, acc) + _socwr_cone_dist_l1(ws, acc)

function _socwr_max_local_dual_exact(ws::_SOCWRPhase2Workspace)
    maxkap = 0.0
    for zone in ws.zones
        model = zone.local_model
        (try; JuMP.has_duals(model); catch; false; end) || continue
        for (Fc, S) in JuMP.list_of_constraint_types(model)
            Fc <: JuMP.VariableRef && continue
            for ci in JuMP.all_constraints(model, Fc, S)
                d = try
                    JuMP.dual(ci)
                catch
                    continue
                end
                if d isa Real
                    ad = abs(d); ad > maxkap && (maxkap = ad)
                elseif d isa AbstractVector
                    nd = LinearAlgebra.norm(d); nd > maxkap && (maxkap = nd)
                end
            end
        end
        # Box-bound multipliers (JuMP exposes them as variable bound duals).
        @inbounds for v in zone.all_vars
            if JuMP.has_lower_bound(v)
                d = try; abs(JuMP.dual(JuMP.LowerBoundRef(v))); catch; 0.0; end
                d > maxkap && (maxkap = d)
            end
            if JuMP.has_upper_bound(v)
                d = try; abs(JuMP.dual(JuMP.UpperBoundRef(v))); catch; 0.0; end
                d > maxkap && (maxkap = d)
            end
        end
    end
    return maxkap
end

# ------------------------------------------------------------
# Danskin slope hook for branch (c): dᵀ (A y − b)
# ------------------------------------------------------------

@inline _socwr_live(z::AladinZoneDataSOCWR, c::Int) = JuMP.value(z.all_vars[c])

function _socwr_consensus_dot(ws::SOCWRStateWorkspaceALADIN, d, acc::F) where {F}
    s = 0.0
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        s += d.v2[k] * (acc(zi, pair.v2_col_i) - acc(zj, pair.v2_col_j))
    end
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        s += d.r[k] * (acc(zi, zi.edge_r_cols[pair.line_local_i]) -
                       acc(zj, zj.edge_r_cols[pair.line_local_j]))
        s += d.t[k] * (acc(zi, zi.edge_t_cols[pair.line_local_i]) -
                       acc(zj, zj.edge_t_cols[pair.line_local_j]))
    end
    return s
end

function _socwr_consensus_dot(ws::SOCWRPowerWorkspaceALADIN, d, acc::F) where {F}
    s = 0.0
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        s += d.p_flow[k] * (_socwr_flow_pt(zi, pair.line_local_i, :p, acc) -
                            _socwr_flow_pt(zj, pair.line_local_j, :p, acc))
        s += d.q_flow[k] * (_socwr_flow_pt(zi, pair.line_local_i, :q, acc) -
                            _socwr_flow_pt(zj, pair.line_local_j, :q, acc))
    end
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        s += d.p_dummy[k] * (_socwr_dummy_flow_pt(pair, zi, :i, :p, acc) -
                             _socwr_dummy_flow_pt(pair, zj, :j, :p, acc))
        s += d.q_dummy[k] * (_socwr_dummy_flow_pt(pair, zi, :i, :q, acc) -
                             _socwr_dummy_flow_pt(pair, zj, :j, :q, acc))
    end
    return s
end

function _socwr_consensus_dot(ws::SOCWRHybridWorkspaceALADIN, d, acc::F) where {F}
    s = 0.0
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        s += d.v2[k] * (acc(zi, pair.v2_col_i) - acc(zj, pair.v2_col_j))
    end
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        s += d.r[k] * (acc(zi, zi.edge_r_cols[pair.line_local_i]) -
                       acc(zj, zj.edge_r_cols[pair.line_local_j]))
        s += d.t[k] * (acc(zi, zi.edge_t_cols[pair.line_local_i]) -
                       acc(zj, zj.edge_t_cols[pair.line_local_j]))
        s += d.p_flow[k] * (_socwr_flow_pt(zi, pair.line_local_i, :p, acc) -
                            _socwr_flow_pt(zj, pair.line_local_j, :p, acc))
        s += d.q_flow[k] * (_socwr_flow_pt(zi, pair.line_local_i, :q, acc) -
                            _socwr_flow_pt(zj, pair.line_local_j, :q, acc))
    end
    @inbounds for (k, pair) in enumerate(ws.split_bus_pairs)
        dp = 0.0; dq = 0.0
        for c in eachindex(pair.zones)
            zone = ws.zones[pair.zones[c]]
            dp += acc(zone, pair.p_inj_cols[c])
            dq += acc(zone, pair.q_inj_cols[c])
        end
        s += d.p_inj[k] * dp + d.q_inj[k] * dq
    end
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        s += d.p_dummy[k] * (_socwr_dummy_flow_pt(pair, zi, :i, :p, acc) -
                             _socwr_dummy_flow_pt(pair, zj, :j, :p, acc))
        s += d.q_dummy[k] * (_socwr_dummy_flow_pt(pair, zi, :i, :q, acc) -
                             _socwr_dummy_flow_pt(pair, zj, :j, :q, acc))
    end
    return s
end

_dual_slope(ws::_SOCWRPhase2Workspace, dirs) =
    _socwr_consensus_dot(ws, dirs, _socwr_live)

function _socwr_consensus_rows(ws::SOCWRStateWorkspaceALADIN, acc::F) where {F}
    r_v2 = Vector{Float64}(undef, length(ws.bus_pairs))
    r_r = Vector{Float64}(undef, length(ws.interzonal_line_pairs))
    r_t = Vector{Float64}(undef, length(ws.interzonal_line_pairs))
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r_v2[k] = acc(zi, pair.v2_col_i) - acc(zj, pair.v2_col_j)
    end
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r_r[k] = acc(zi, zi.edge_r_cols[pair.line_local_i]) -
                 acc(zj, zj.edge_r_cols[pair.line_local_j])
        r_t[k] = acc(zi, zi.edge_t_cols[pair.line_local_i]) -
                 acc(zj, zj.edge_t_cols[pair.line_local_j])
    end
    return NamedTuple{keys(ws.vars)}((v2 = r_v2, r = r_r, t = r_t))
end

function _socwr_consensus_rows(ws::SOCWRPowerWorkspaceALADIN, acc::F) where {F}
    n_e = length(ws.interzonal_line_pairs)
    n_d = length(ws.dummy_line_pairs)
    r_p = Vector{Float64}(undef, n_e)
    r_q = Vector{Float64}(undef, n_e)
    r_pd = Vector{Float64}(undef, n_d)
    r_qd = Vector{Float64}(undef, n_d)
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r_p[k] = _socwr_flow_pt(zi, pair.line_local_i, :p, acc) -
                 _socwr_flow_pt(zj, pair.line_local_j, :p, acc)
        r_q[k] = _socwr_flow_pt(zi, pair.line_local_i, :q, acc) -
                 _socwr_flow_pt(zj, pair.line_local_j, :q, acc)
    end
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r_pd[k] = _socwr_dummy_flow_pt(pair, zi, :i, :p, acc) -
                  _socwr_dummy_flow_pt(pair, zj, :j, :p, acc)
        r_qd[k] = _socwr_dummy_flow_pt(pair, zi, :i, :q, acc) -
                  _socwr_dummy_flow_pt(pair, zj, :j, :q, acc)
    end
    return NamedTuple{keys(ws.vars)}((
        p_flow = r_p, q_flow = r_q, p_dummy = r_pd, q_dummy = r_qd,
    ))
end

function _socwr_consensus_rows(ws::SOCWRHybridWorkspaceALADIN, acc::F) where {F}
    n_e = length(ws.interzonal_line_pairs)
    n_s = length(ws.split_bus_pairs)
    n_d = length(ws.dummy_line_pairs)
    r_v2 = Vector{Float64}(undef, length(ws.bus_pairs))
    r_r = Vector{Float64}(undef, n_e)
    r_t = Vector{Float64}(undef, n_e)
    r_p = Vector{Float64}(undef, n_e)
    r_q = Vector{Float64}(undef, n_e)
    r_pi = Vector{Float64}(undef, n_s)
    r_qi = Vector{Float64}(undef, n_s)
    r_pd = Vector{Float64}(undef, n_d)
    r_qd = Vector{Float64}(undef, n_d)
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r_v2[k] = acc(zi, pair.v2_col_i) - acc(zj, pair.v2_col_j)
    end
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r_r[k] = acc(zi, zi.edge_r_cols[pair.line_local_i]) -
                 acc(zj, zj.edge_r_cols[pair.line_local_j])
        r_t[k] = acc(zi, zi.edge_t_cols[pair.line_local_i]) -
                 acc(zj, zj.edge_t_cols[pair.line_local_j])
        r_p[k] = _socwr_flow_pt(zi, pair.line_local_i, :p, acc) -
                 _socwr_flow_pt(zj, pair.line_local_j, :p, acc)
        r_q[k] = _socwr_flow_pt(zi, pair.line_local_i, :q, acc) -
                 _socwr_flow_pt(zj, pair.line_local_j, :q, acc)
    end
    @inbounds for (k, pair) in enumerate(ws.split_bus_pairs)
        dp = 0.0; dq = 0.0
        for c in eachindex(pair.zones)
            zone = ws.zones[pair.zones[c]]
            dp += acc(zone, pair.p_inj_cols[c])
            dq += acc(zone, pair.q_inj_cols[c])
        end
        r_pi[k] = dp
        r_qi[k] = dq
    end
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r_pd[k] = _socwr_dummy_flow_pt(pair, zi, :i, :p, acc) -
                  _socwr_dummy_flow_pt(pair, zj, :j, :p, acc)
        r_qd[k] = _socwr_dummy_flow_pt(pair, zi, :i, :q, acc) -
                  _socwr_dummy_flow_pt(pair, zj, :j, :q, acc)
    end
    return NamedTuple{keys(ws.vars)}((
        v2 = r_v2, p_flow = r_p, q_flow = r_q, r = r_r, t = r_t,
        p_inj = r_pi, q_inj = r_qi, p_dummy = r_pd, q_dummy = r_qd,
    ))
end

_dual_residual(ws::_SOCWRPhase2Workspace) =
    _socwr_consensus_rows(ws, _socwr_live)

# ------------------------------------------------------------
# λ̄ / κ̄ (paper: "sufficiently large" fixed constants)
# ------------------------------------------------------------

function _alg3_original_update_weights!(ws::_SOCWRPhase2Workspace)
    safety = _ALG3_ORIG_SAFETY[]
    floor_ = _ALG3_ORIG_FLOOR[]
    if ws.alg3_lam_bar > 0.0
        # Paper literal: fixed, user-chosen "sufficiently large" λ̄ = κ̄.
        ws.lam_bar = ws.alg3_lam_bar
        ws.kap_bar = ws.alg3_lam_bar
        maxlam = 0.0
        for state in values(ws.vars)
            @inbounds for k in eachindex(state.lambda)
                maxlam = max(maxlam, abs(state.lambda[k]))
            end
        end
        if maxlam > ws.lam_bar
            ws.hist.meta[:lam_bar_violations] =
                Int(get(ws.hist.meta, :lam_bar_violations, 0)) + 1
            @warn "ALADIN :alg3_original (SOCWR) — |λ|∞ = $(maxlam) exceeds the fixed " *
                  "L1 penalty weight λ̄ = $(ws.lam_bar) (iter $(ws.iter))." maxlog = 1
        end
        maxkap = _socwr_max_local_dual_exact(ws)
        if maxkap > ws.kap_bar
            ws.hist.meta[:kap_bar_violations] =
                Int(get(ws.hist.meta, :kap_bar_violations, 0)) + 1
            @warn "ALADIN :alg3_original (SOCWR) — max‖κ_local‖ = $(maxkap) exceeds the " *
                  "fixed local-penalty weight κ̄ = $(ws.kap_bar) (iter $(ws.iter))." maxlog = 1
        end
        return nothing
    end
    maxlam = 0.0
    for state in values(ws.vars)
        @inbounds for k in eachindex(state.lambda)
            maxlam = max(maxlam, abs(state.lambda[k]), abs(state.nu[k]))
        end
    end
    if ws.lam_bar < max(maxlam, floor_)
        ws.lam_bar = max(safety * maxlam, floor_)
    end
    maxkap = _socwr_max_local_dual_exact(ws)
    if ws.kap_bar < max(maxkap, floor_)
        ws.kap_bar = max(safety * maxkap, floor_)
    end
    return nothing
end

# Φ(x) = Σ fᵢ(xᵢ) + λ̄ ‖Ax − b‖₁ + κ̄ · [scalar ℓ1 viol + Σ dist(·, K)]
_socwr_phi_exact(ws::_SOCWRPhase2Workspace, acc::F) where {F} =
    _socwr_cost_sum(ws, acc) + ws.lam_bar * _socwr_resid_l1(ws, acc) +
    ws.kap_bar * _socwr_local_infeas_exact(ws, acc)

function _alg3_original_step!(ws::_SOCWRPhase2Workspace; gamma::Float64 = 1e-4)
    accx = (z, c) -> z.x_prev[c]
    accy = (z, c) -> z.y[c]
    accyd = (z, c) -> z.y[c] + z.delta[c]

    _alg3_original_update_weights!(ws)

    # (6.2) right-hand side: γ ( Σ (ρ/2)‖y − x‖²_Σ + λ̄ ‖Ay − b‖₁ )
    rhs = gamma * (_socwr_prox(ws) + ws.lam_bar * _socwr_resid_l1(ws, accy))
    phi_x = _socwr_phi_exact(ws, accx)

    # (a) α₁ = α₂ = α₃ = 1, x⁺ = y + Δy
    phi_a = _socwr_phi_exact(ws, accyd)
    if isfinite(phi_a) && phi_x - phi_a >= rhs
        ws.alpha1 = 1.0; ws.alpha2 = 1.0; ws.alpha3 = 1.0
        ws.alg3_branch = :a
        _socwr_alg3o_record!(ws, :a, :none)
        return nothing
    end

    # (b) α₁ = 1, α₂ = α₃ = 0, x⁺ = y
    phi_b = _socwr_phi_exact(ws, accy)
    if isfinite(phi_b) && phi_x - phi_b >= rhs
        ws.alpha1 = 1.0; ws.alpha2 = 0.0; ws.alpha3 = 0.0
        ws.alg3_branch = :b
        _socwr_alg3o_record!(ws, :b, :none)
        return nothing
    end

    # (c) α₁ = α₂ = 0, α₃ from the Wolfe dual line search on V_ρ
    ws.alg3_branch = :c
    ws.alg3_c = _branch_c_dual_linesearch_wolfe!(ws)
    _socwr_alg3o_record!(ws, :c, ws.alg3_c)
    return nothing
end

@inline function _socwr_alg3o_record!(ws::_SOCWRPhase2Workspace, branch::Symbol, status_c::Symbol)
    log = get!(ws.hist.meta, :alg3_original,
        Vector{Tuple{Symbol, Float64, Float64, Symbol}}())
    push!(log::Vector{Tuple{Symbol, Float64, Float64, Symbol}},
        (branch, ws.alpha3, ws.lam_bar, status_c))
    return nothing
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
    ws.globalization === :alg3 && return _alg3_step!(ws)
    ws.globalization === :alg3_original && return _alg3_original_step!(ws)
    # :full_step — undamped α₁=α₂=1, α₃=alpha3_cfg (per-iter reset in run.jl).
    ws.globalization === :full_step && return nothing

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
