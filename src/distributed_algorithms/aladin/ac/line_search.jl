# ============================================================
# ALADIN-AC — Phase 2 line search (Armijo backtracking)
# ============================================================

const _ACPhase2Workspace = Union{
    ACStateWorkspaceALADIN,
    ACPowerWorkspaceALADIN,
    ACHybridWorkspaceALADIN,
}

# ============================================================
# Per-zone cost trial evaluations
# ============================================================

@inline function _ac_cost_at_trial(zone::AladinZoneDataAC, alpha::Float64)
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

@inline function _ac_cost_dir_deriv_0(zone::AladinZoneDataAC)
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
# Interzonal AC flow trial evaluations
# ============================================================

@inline function _ac_pf_at_trial(zone::AladinZoneDataAC, le::Int, alpha::Float64)
    fcol_theta = zone.edge_from_theta_col[le]
    tcol_theta = zone.edge_to_theta_col[le]
    fcol_v = zone.edge_from_vm_col[le]
    tcol_v = zone.edge_to_vm_col[le]
    theta_f = zone.y[fcol_theta] + alpha * zone.delta[fcol_theta]
    theta_t = zone.y[tcol_theta] + alpha * zone.delta[tcol_theta]
    vf = zone.y[fcol_v] + alpha * zone.delta[fcol_v]
    vt = zone.y[tcol_v] + alpha * zone.delta[tcol_v]
    g = zone.edge_g[le]
    b = zone.edge_b[le]
    tau = zone.edge_inv_tap[le]
    delta_theta = theta_f - theta_t - zone.edge_tap_angle[le]
    return g * (tau * vf)^2 -
           tau * vf * vt * (g * cos(delta_theta) + b * sin(delta_theta))
end

@inline function _ac_qf_at_trial(zone::AladinZoneDataAC, le::Int, alpha::Float64)
    fcol_theta = zone.edge_from_theta_col[le]
    tcol_theta = zone.edge_to_theta_col[le]
    fcol_v = zone.edge_from_vm_col[le]
    tcol_v = zone.edge_to_vm_col[le]
    theta_f = zone.y[fcol_theta] + alpha * zone.delta[fcol_theta]
    theta_t = zone.y[tcol_theta] + alpha * zone.delta[tcol_theta]
    vf = zone.y[fcol_v] + alpha * zone.delta[fcol_v]
    vt = zone.y[tcol_v] + alpha * zone.delta[tcol_v]
    g = zone.edge_g[le]
    b = zone.edge_b[le]
    b_sh = zone.edge_b_sh[le]
    tau = zone.edge_inv_tap[le]
    delta_theta = theta_f - theta_t - zone.edge_tap_angle[le]
    return -(b + b_sh / 2) * (tau * vf)^2 -
           tau * vf * vt * (-b * cos(delta_theta) + g * sin(delta_theta))
end

# ============================================================
# Interzonal AC directional derivatives at alpha = 0
# ============================================================

@inline function _ac_pf_dir_deriv(zone::AladinZoneDataAC, le::Int)
    fcol_theta = zone.edge_from_theta_col[le]
    tcol_theta = zone.edge_to_theta_col[le]
    fcol_v = zone.edge_from_vm_col[le]
    tcol_v = zone.edge_to_vm_col[le]
    theta_f = zone.y[fcol_theta]
    theta_t = zone.y[tcol_theta]
    vf = zone.y[fcol_v]
    vt = zone.y[tcol_v]
    dtheta_f = zone.delta[fcol_theta]
    dtheta_t = zone.delta[tcol_theta]
    dvf = zone.delta[fcol_v]
    dvt = zone.delta[tcol_v]
    g = zone.edge_g[le]
    b = zone.edge_b[le]
    tau = zone.edge_inv_tap[le]
    delta_theta = theta_f - theta_t - zone.edge_tap_angle[le]
    c_delta = cos(delta_theta)
    s_delta = sin(delta_theta)
    dpf_dvf = 2 * g * tau * tau * vf - tau * vt * (g * c_delta + b * s_delta)
    dpf_dvt = -tau * vf * (g * c_delta + b * s_delta)
    dpf_dtheta_f = tau * vf * vt * (g * s_delta - b * c_delta)
    return dpf_dvf * dvf + dpf_dvt * dvt + dpf_dtheta_f * (dtheta_f - dtheta_t)
end

@inline function _ac_qf_dir_deriv(zone::AladinZoneDataAC, le::Int)
    fcol_theta = zone.edge_from_theta_col[le]
    tcol_theta = zone.edge_to_theta_col[le]
    fcol_v = zone.edge_from_vm_col[le]
    tcol_v = zone.edge_to_vm_col[le]
    theta_f = zone.y[fcol_theta]
    theta_t = zone.y[tcol_theta]
    vf = zone.y[fcol_v]
    vt = zone.y[tcol_v]
    dtheta_f = zone.delta[fcol_theta]
    dtheta_t = zone.delta[tcol_theta]
    dvf = zone.delta[fcol_v]
    dvt = zone.delta[tcol_v]
    g = zone.edge_g[le]
    b = zone.edge_b[le]
    b_sh = zone.edge_b_sh[le]
    tau = zone.edge_inv_tap[le]
    delta_theta = theta_f - theta_t - zone.edge_tap_angle[le]
    c_delta = cos(delta_theta)
    s_delta = sin(delta_theta)
    dqf_dvf = -2 * (b + b_sh / 2) * tau * tau * vf -
              tau * vt * (-b * c_delta + g * s_delta)
    dqf_dvt = -tau * vf * (-b * c_delta + g * s_delta)
    dqf_dtheta_f = -tau * vf * vt * (b * s_delta + g * c_delta)
    return dqf_dvf * dvf + dqf_dvt * dvt + dqf_dtheta_f * (dtheta_f - dtheta_t)
end

# ============================================================
# Dummy line AC flow trial evaluations
# ============================================================

@inline function _ac_dummy_pf_at_trial(
    pair::DummyLinePairAC, zone::AladinZoneDataAC, side::Symbol, alpha::Float64,
)
    if side === :i
        fc_theta, tc_theta = pair.from_theta_col_i, pair.to_theta_col_i
        fc_v, tc_v = pair.from_vm_col_i, pair.to_vm_col_i
        g, b, tau, ta = pair.g_i, pair.b_i, pair.inv_tap_i, pair.tap_angle_i
    else
        fc_theta, tc_theta = pair.from_theta_col_j, pair.to_theta_col_j
        fc_v, tc_v = pair.from_vm_col_j, pair.to_vm_col_j
        g, b, tau, ta = pair.g_j, pair.b_j, pair.inv_tap_j, pair.tap_angle_j
    end
    theta_f = zone.y[fc_theta] + alpha * zone.delta[fc_theta]
    theta_t = zone.y[tc_theta] + alpha * zone.delta[tc_theta]
    vf = zone.y[fc_v] + alpha * zone.delta[fc_v]
    vt = zone.y[tc_v] + alpha * zone.delta[tc_v]
    delta_theta = theta_f - theta_t - ta
    c = cos(delta_theta)
    s = sin(delta_theta)
    if side === :i
        return -g * vt^2 + tau * vf * vt * (g * c - b * s)
    end
    return g * (tau * vf)^2 - tau * vf * vt * (g * c + b * s)
end

@inline function _ac_dummy_qf_at_trial(
    pair::DummyLinePairAC, zone::AladinZoneDataAC, side::Symbol, alpha::Float64,
)
    if side === :i
        fc_theta, tc_theta = pair.from_theta_col_i, pair.to_theta_col_i
        fc_v, tc_v = pair.from_vm_col_i, pair.to_vm_col_i
        g, b, b_sh = pair.g_i, pair.b_i, pair.b_sh_i
        tau, ta = pair.inv_tap_i, pair.tap_angle_i
    else
        fc_theta, tc_theta = pair.from_theta_col_j, pair.to_theta_col_j
        fc_v, tc_v = pair.from_vm_col_j, pair.to_vm_col_j
        g, b, b_sh = pair.g_j, pair.b_j, pair.b_sh_j
        tau, ta = pair.inv_tap_j, pair.tap_angle_j
    end
    theta_f = zone.y[fc_theta] + alpha * zone.delta[fc_theta]
    theta_t = zone.y[tc_theta] + alpha * zone.delta[tc_theta]
    vf = zone.y[fc_v] + alpha * zone.delta[fc_v]
    vt = zone.y[tc_v] + alpha * zone.delta[tc_v]
    delta_theta = theta_f - theta_t - ta
    c = cos(delta_theta)
    s = sin(delta_theta)
    if side === :i
        # half-A: TO-side reactive across D (see _ac_dummy_pf_at_trial).
        return (b + b_sh / 2) * vt^2 + tau * vf * vt * (-b * c - g * s)
    end
    return -(b + b_sh / 2) * (tau * vf)^2 - tau * vf * vt * (-b * c + g * s)
end

# ============================================================
# Dummy line AC directional derivatives at alpha = 0
# ============================================================

@inline function _ac_dummy_pf_dir_deriv(
    pair::DummyLinePairAC, zone::AladinZoneDataAC, side::Symbol,
)
    if side === :i
        fc_theta, tc_theta = pair.from_theta_col_i, pair.to_theta_col_i
        fc_v, tc_v = pair.from_vm_col_i, pair.to_vm_col_i
        g, b, tau, ta = pair.g_i, pair.b_i, pair.inv_tap_i, pair.tap_angle_i
    else
        fc_theta, tc_theta = pair.from_theta_col_j, pair.to_theta_col_j
        fc_v, tc_v = pair.from_vm_col_j, pair.to_vm_col_j
        g, b, tau, ta = pair.g_j, pair.b_j, pair.inv_tap_j, pair.tap_angle_j
    end
    theta_f = zone.y[fc_theta]
    theta_t = zone.y[tc_theta]
    vf = zone.y[fc_v]
    vt = zone.y[tc_v]
    dtheta_f = zone.delta[fc_theta]
    dtheta_t = zone.delta[tc_theta]
    dvf = zone.delta[fc_v]
    dvt = zone.delta[tc_v]
    delta_theta = theta_f - theta_t - ta
    c_delta = cos(delta_theta)
    s_delta = sin(delta_theta)
    if side === :i
        dpt_dvf = tau * vt * (g * c_delta - b * s_delta)
        dpt_dvt = -2 * g * vt + tau * vf * (g * c_delta - b * s_delta)
        dpt_dtheta_f = tau * vf * vt * (-g * s_delta - b * c_delta)
        return dpt_dvf * dvf + dpt_dvt * dvt + dpt_dtheta_f * (dtheta_f - dtheta_t)
    end
    dpf_dvf = 2 * g * tau * tau * vf - tau * vt * (g * c_delta + b * s_delta)
    dpf_dvt = -tau * vf * (g * c_delta + b * s_delta)
    dpf_dtheta_f = tau * vf * vt * (g * s_delta - b * c_delta)
    return dpf_dvf * dvf + dpf_dvt * dvt + dpf_dtheta_f * (dtheta_f - dtheta_t)
end

@inline function _ac_dummy_qf_dir_deriv(
    pair::DummyLinePairAC, zone::AladinZoneDataAC, side::Symbol,
)
    if side === :i
        fc_theta, tc_theta = pair.from_theta_col_i, pair.to_theta_col_i
        fc_v, tc_v = pair.from_vm_col_i, pair.to_vm_col_i
        g, b, b_sh = pair.g_i, pair.b_i, pair.b_sh_i
        tau, ta = pair.inv_tap_i, pair.tap_angle_i
    else
        fc_theta, tc_theta = pair.from_theta_col_j, pair.to_theta_col_j
        fc_v, tc_v = pair.from_vm_col_j, pair.to_vm_col_j
        g, b, b_sh = pair.g_j, pair.b_j, pair.b_sh_j
        tau, ta = pair.inv_tap_j, pair.tap_angle_j
    end
    theta_f = zone.y[fc_theta]
    theta_t = zone.y[tc_theta]
    vf = zone.y[fc_v]
    vt = zone.y[tc_v]
    dtheta_f = zone.delta[fc_theta]
    dtheta_t = zone.delta[tc_theta]
    dvf = zone.delta[fc_v]
    dvt = zone.delta[tc_v]
    delta_theta = theta_f - theta_t - ta
    c_delta = cos(delta_theta)
    s_delta = sin(delta_theta)
    if side === :i
        # half-A: derivatives of the TO-side reactive (see _ac_dummy_pf_dir_deriv).
        b_to = b + b_sh / 2
        dqt_dvf = tau * vt * (-b * c_delta - g * s_delta)
        dqt_dvt = 2 * b_to * vt + tau * vf * (-b * c_delta - g * s_delta)
        dqt_dtheta_f = tau * vf * vt * (b * s_delta - g * c_delta)
        return dqt_dvf * dvf + dqt_dvt * dvt + dqt_dtheta_f * (dtheta_f - dtheta_t)
    end
    dqf_dvf = -2 * (b + b_sh / 2) * tau * tau * vf -
              tau * vt * (-b * c_delta + g * s_delta)
    dqf_dvt = -tau * vf * (-b * c_delta + g * s_delta)
    dqf_dtheta_f = -tau * vf * vt * (b * s_delta + g * c_delta)
    return dqf_dvf * dvf + dqf_dvt * dvt + dqf_dtheta_f * (dtheta_f - dtheta_t)
end

# ============================================================
# Merit function phi(alpha)
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function _ac_merit(ws::ACStateWorkspaceALADIN, alpha::Float64)
    phi = 0.0
    for zone in ws.zones
        phi += _ac_cost_at_trial(zone, alpha)
    end
    r_sq = 0.0
    lambda_r = 0.0
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        theta_i = zi.y[pair.theta_col_i] + alpha * zi.delta[pair.theta_col_i]
        theta_j = zj.y[pair.theta_col_j] + alpha * zj.delta[pair.theta_col_j]
        v_i = zi.y[pair.vm_col_i] + alpha * zi.delta[pair.vm_col_i]
        v_j = zj.y[pair.vm_col_j] + alpha * zj.delta[pair.vm_col_j]
        dtheta = theta_i - theta_j
        dv = v_i - v_j
        r_sq += dtheta * dtheta + dv * dv
        lambda_r += ws.vars.theta.lambda[k] * dtheta + ws.vars.v_mag.lambda[k] * dv
    end
    return phi + lambda_r + 0.5 * ws.mu[] * r_sq
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function _ac_merit(ws::ACPowerWorkspaceALADIN, alpha::Float64)
    phi = 0.0
    for zone in ws.zones
        phi += _ac_cost_at_trial(zone, alpha)
    end
    r_sq = 0.0
    lambda_r = 0.0
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = _ac_pf_at_trial(zi, pair.line_local_i, alpha)
        pj_val = _ac_pf_at_trial(zj, pair.line_local_j, alpha)
        qi_val = _ac_qf_at_trial(zi, pair.line_local_i, alpha)
        qj_val = _ac_qf_at_trial(zj, pair.line_local_j, alpha)
        dp = pi_val - pj_val
        dq = qi_val - qj_val
        r_sq += dp * dp + dq * dq
        lambda_r += ws.vars.p_flow.lambda[k] * dp + ws.vars.q_flow.lambda[k] * dq
    end
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = _ac_dummy_pf_at_trial(pair, zi, :i, alpha)
        pj_val = _ac_dummy_pf_at_trial(pair, zj, :j, alpha)
        qi_val = _ac_dummy_qf_at_trial(pair, zi, :i, alpha)
        qj_val = _ac_dummy_qf_at_trial(pair, zj, :j, alpha)
        dp = pi_val - pj_val
        dq = qi_val - qj_val
        r_sq += dp * dp + dq * dq
        lambda_r += ws.vars.p_dummy.lambda[k] * dp + ws.vars.q_dummy.lambda[k] * dq
    end
    return phi + lambda_r + 0.5 * ws.mu[] * r_sq
end

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

function _ac_merit(ws::ACHybridWorkspaceALADIN, alpha::Float64)
    phi = 0.0
    for zone in ws.zones
        phi += _ac_cost_at_trial(zone, alpha)
    end
    r_sq = 0.0
    lambda_r = 0.0
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        theta_i = zi.y[pair.theta_col_i] + alpha * zi.delta[pair.theta_col_i]
        theta_j = zj.y[pair.theta_col_j] + alpha * zj.delta[pair.theta_col_j]
        v_i = zi.y[pair.vm_col_i] + alpha * zi.delta[pair.vm_col_i]
        v_j = zj.y[pair.vm_col_j] + alpha * zj.delta[pair.vm_col_j]
        dtheta = theta_i - theta_j
        dv = v_i - v_j
        r_sq += dtheta * dtheta + dv * dv
        lambda_r += ws.vars.theta.lambda[k] * dtheta + ws.vars.v_mag.lambda[k] * dv
    end
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = _ac_pf_at_trial(zi, pair.line_local_i, alpha)
        pj_val = _ac_pf_at_trial(zj, pair.line_local_j, alpha)
        qi_val = _ac_qf_at_trial(zi, pair.line_local_i, alpha)
        qj_val = _ac_qf_at_trial(zj, pair.line_local_j, alpha)
        dp = pi_val - pj_val
        dq = qi_val - qj_val
        r_sq += dp * dp + dq * dq
        lambda_r += ws.vars.p_flow.lambda[k] * dp + ws.vars.q_flow.lambda[k] * dq
    end
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = _ac_dummy_pf_at_trial(pair, zi, :i, alpha)
        pj_val = _ac_dummy_pf_at_trial(pair, zj, :j, alpha)
        qi_val = _ac_dummy_qf_at_trial(pair, zi, :i, alpha)
        qj_val = _ac_dummy_qf_at_trial(pair, zj, :j, alpha)
        dp = pi_val - pj_val
        dq = qi_val - qj_val
        r_sq += dp * dp + dq * dq
        lambda_r += ws.vars.p_dummy.lambda[k] * dp + ws.vars.q_dummy.lambda[k] * dq
    end
    @inbounds for (k, pair) in enumerate(ws.split_bus_pairs)
        dp = 0.0
        dq = 0.0
        for c in eachindex(pair.zones)
            zone = ws.zones[pair.zones[c]]
            pcol = pair.p_inj_cols[c]
            qcol = pair.q_inj_cols[c]
            dp += zone.y[pcol] + alpha * zone.delta[pcol]
            dq += zone.y[qcol] + alpha * zone.delta[qcol]
        end
        r_sq += dp * dp + dq * dq
        lambda_r += ws.vars.p_inj.lambda[k] * dp + ws.vars.q_inj.lambda[k] * dq
    end
    return phi + lambda_r + 0.5 * ws.mu[] * r_sq
end

# ============================================================
# Directional derivative dphi/dalpha at alpha = 0
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function _ac_dphi_0(ws::ACStateWorkspaceALADIN)
    dphi = 0.0
    for zone in ws.zones
        dphi += _ac_cost_dir_deriv_0(zone)
    end
    mu = ws.mu[]
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        r0_theta = zi.y[pair.theta_col_i] - zj.y[pair.theta_col_j]
        dr_theta = zi.delta[pair.theta_col_i] - zj.delta[pair.theta_col_j]
        r0_v = zi.y[pair.vm_col_i] - zj.y[pair.vm_col_j]
        dr_v = zi.delta[pair.vm_col_i] - zj.delta[pair.vm_col_j]
        dphi += (ws.vars.theta.lambda[k] + mu * r0_theta) * dr_theta
        dphi += (ws.vars.v_mag.lambda[k] + mu * r0_v) * dr_v
    end
    return dphi
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function _ac_dphi_0(ws::ACPowerWorkspaceALADIN)
    dphi = 0.0
    for zone in ws.zones
        dphi += _ac_cost_dir_deriv_0(zone)
    end
    mu = ws.mu[]
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        r0_p =
            _ac_pf_at_trial(zi, pair.line_local_i, 0.0) -
            _ac_pf_at_trial(zj, pair.line_local_j, 0.0)
        dr_p =
            _ac_pf_dir_deriv(zi, pair.line_local_i) -
            _ac_pf_dir_deriv(zj, pair.line_local_j)
        r0_q =
            _ac_qf_at_trial(zi, pair.line_local_i, 0.0) -
            _ac_qf_at_trial(zj, pair.line_local_j, 0.0)
        dr_q =
            _ac_qf_dir_deriv(zi, pair.line_local_i) -
            _ac_qf_dir_deriv(zj, pair.line_local_j)
        dphi += (ws.vars.p_flow.lambda[k] + mu * r0_p) * dr_p
        dphi += (ws.vars.q_flow.lambda[k] + mu * r0_q) * dr_q
    end
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        r0_p =
            _ac_dummy_pf_at_trial(pair, zi, :i, 0.0) -
            _ac_dummy_pf_at_trial(pair, zj, :j, 0.0)
        dr_p = _ac_dummy_pf_dir_deriv(pair, zi, :i) -
               _ac_dummy_pf_dir_deriv(pair, zj, :j)
        r0_q =
            _ac_dummy_qf_at_trial(pair, zi, :i, 0.0) -
            _ac_dummy_qf_at_trial(pair, zj, :j, 0.0)
        dr_q = _ac_dummy_qf_dir_deriv(pair, zi, :i) -
               _ac_dummy_qf_dir_deriv(pair, zj, :j)
        dphi += (ws.vars.p_dummy.lambda[k] + mu * r0_p) * dr_p
        dphi += (ws.vars.q_dummy.lambda[k] + mu * r0_q) * dr_q
    end
    return dphi
end

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

function _ac_dphi_0(ws::ACHybridWorkspaceALADIN)
    dphi = 0.0
    for zone in ws.zones
        dphi += _ac_cost_dir_deriv_0(zone)
    end
    mu = ws.mu[]
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        r0_theta = zi.y[pair.theta_col_i] - zj.y[pair.theta_col_j]
        dr_theta = zi.delta[pair.theta_col_i] - zj.delta[pair.theta_col_j]
        r0_v = zi.y[pair.vm_col_i] - zj.y[pair.vm_col_j]
        dr_v = zi.delta[pair.vm_col_i] - zj.delta[pair.vm_col_j]
        dphi += (ws.vars.theta.lambda[k] + mu * r0_theta) * dr_theta
        dphi += (ws.vars.v_mag.lambda[k] + mu * r0_v) * dr_v
    end
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        r0_p =
            _ac_pf_at_trial(zi, pair.line_local_i, 0.0) -
            _ac_pf_at_trial(zj, pair.line_local_j, 0.0)
        dr_p =
            _ac_pf_dir_deriv(zi, pair.line_local_i) -
            _ac_pf_dir_deriv(zj, pair.line_local_j)
        r0_q =
            _ac_qf_at_trial(zi, pair.line_local_i, 0.0) -
            _ac_qf_at_trial(zj, pair.line_local_j, 0.0)
        dr_q =
            _ac_qf_dir_deriv(zi, pair.line_local_i) -
            _ac_qf_dir_deriv(zj, pair.line_local_j)
        dphi += (ws.vars.p_flow.lambda[k] + mu * r0_p) * dr_p
        dphi += (ws.vars.q_flow.lambda[k] + mu * r0_q) * dr_q
    end
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        r0_p =
            _ac_dummy_pf_at_trial(pair, zi, :i, 0.0) -
            _ac_dummy_pf_at_trial(pair, zj, :j, 0.0)
        dr_p = _ac_dummy_pf_dir_deriv(pair, zi, :i) -
               _ac_dummy_pf_dir_deriv(pair, zj, :j)
        r0_q =
            _ac_dummy_qf_at_trial(pair, zi, :i, 0.0) -
            _ac_dummy_qf_at_trial(pair, zj, :j, 0.0)
        dr_q = _ac_dummy_qf_dir_deriv(pair, zi, :i) -
               _ac_dummy_qf_dir_deriv(pair, zj, :j)
        dphi += (ws.vars.p_dummy.lambda[k] + mu * r0_p) * dr_p
        dphi += (ws.vars.q_dummy.lambda[k] + mu * r0_q) * dr_q
    end
    @inbounds for (k, pair) in enumerate(ws.split_bus_pairs)
        r0_p = 0.0
        dr_p = 0.0
        r0_q = 0.0
        dr_q = 0.0
        for c in eachindex(pair.zones)
            zone = ws.zones[pair.zones[c]]
            pcol = pair.p_inj_cols[c]
            qcol = pair.q_inj_cols[c]
            r0_p += zone.y[pcol]
            dr_p += zone.delta[pcol]
            r0_q += zone.y[qcol]
            dr_q += zone.delta[qcol]
        end
        dphi += (ws.vars.p_inj.lambda[k] + mu * r0_p) * dr_p
        dphi += (ws.vars.q_inj.lambda[k] + mu * r0_q) * dr_q
    end
    return dphi
end

# ============================================================
# HFD Algorithm 3 globalization — L1 exact penalty Φ, steps (a)/(b)/(c)
# ============================================================

@inline function _ac_cost(zone::AladinZoneDataAC, acc::F) where {F}
    total = 0.0
    @inbounds for k in eachindex(zone.pg_cols)
        col = zone.pg_cols[k]
        xk = acc(zone, col)
        total += 0.5 * zone.pg_c2_sbase2[k] * xk * xk + zone.pg_c1_sbase[k] * xk
    end
    return total
end

@inline function _ac_pf(zone::AladinZoneDataAC, le::Int, acc::F) where {F}
    theta_f = acc(zone, zone.edge_from_theta_col[le])
    theta_t = acc(zone, zone.edge_to_theta_col[le])
    vf = acc(zone, zone.edge_from_vm_col[le])
    vt = acc(zone, zone.edge_to_vm_col[le])
    g = zone.edge_g[le]
    b = zone.edge_b[le]
    tau = zone.edge_inv_tap[le]
    dth = theta_f - theta_t - zone.edge_tap_angle[le]
    return g * (tau * vf)^2 - tau * vf * vt * (g * cos(dth) + b * sin(dth))
end

@inline function _ac_qf(zone::AladinZoneDataAC, le::Int, acc::F) where {F}
    theta_f = acc(zone, zone.edge_from_theta_col[le])
    theta_t = acc(zone, zone.edge_to_theta_col[le])
    vf = acc(zone, zone.edge_from_vm_col[le])
    vt = acc(zone, zone.edge_to_vm_col[le])
    g = zone.edge_g[le]
    b = zone.edge_b[le]
    b_sh = zone.edge_b_sh[le]
    tau = zone.edge_inv_tap[le]
    dth = theta_f - theta_t - zone.edge_tap_angle[le]
    return -(b + b_sh / 2) * (tau * vf)^2 - tau * vf * vt * (-b * cos(dth) + g * sin(dth))
end

@inline function _ac_dpf(pair::DummyLinePairAC, zone::AladinZoneDataAC, side::Symbol, acc::F) where {F}
    if side === :i
        fc, tc = pair.from_theta_col_i, pair.to_theta_col_i
        fv, tv = pair.from_vm_col_i, pair.to_vm_col_i
        g, b, tau, ta = pair.g_i, pair.b_i, pair.inv_tap_i, pair.tap_angle_i
    else
        fc, tc = pair.from_theta_col_j, pair.to_theta_col_j
        fv, tv = pair.from_vm_col_j, pair.to_vm_col_j
        g, b, tau, ta = pair.g_j, pair.b_j, pair.inv_tap_j, pair.tap_angle_j
    end
    theta_f = acc(zone, fc); theta_t = acc(zone, tc)
    vf = acc(zone, fv); vt = acc(zone, tv)
    dth = theta_f - theta_t - ta
    c = cos(dth); s = sin(dth)
    side === :i && return -g * vt^2 + tau * vf * vt * (g * c - b * s)
    return g * (tau * vf)^2 - tau * vf * vt * (g * c + b * s)
end

@inline function _ac_dqf(pair::DummyLinePairAC, zone::AladinZoneDataAC, side::Symbol, acc::F) where {F}
    if side === :i
        fc, tc = pair.from_theta_col_i, pair.to_theta_col_i
        fv, tv = pair.from_vm_col_i, pair.to_vm_col_i
        g, b, b_sh = pair.g_i, pair.b_i, pair.b_sh_i
        tau, ta = pair.inv_tap_i, pair.tap_angle_i
    else
        fc, tc = pair.from_theta_col_j, pair.to_theta_col_j
        fv, tv = pair.from_vm_col_j, pair.to_vm_col_j
        g, b, b_sh = pair.g_j, pair.b_j, pair.b_sh_j
        tau, ta = pair.inv_tap_j, pair.tap_angle_j
    end
    theta_f = acc(zone, fc); theta_t = acc(zone, tc)
    vf = acc(zone, fv); vt = acc(zone, tv)
    dth = theta_f - theta_t - ta
    c = cos(dth); s = sin(dth)
    side === :i && return (b + b_sh / 2) * vt^2 + tau * vf * vt * (-b * c - g * s)
    return -(b + b_sh / 2) * (tau * vf)^2 - tau * vf * vt * (-b * c + g * s)
end

_ac_cost_sum(ws::_ACPhase2Workspace, acc::F) where {F} =
    sum(_ac_cost(zone, acc) for zone in ws.zones)

# ℓ1 consensus residual ‖Ax − b‖₁ at the point defined by `acc` (b = 0).
function _ac_resid_l1(ws::ACStateWorkspaceALADIN, acc::F) where {F}
    r = 0.0
    @inbounds for pair in ws.bus_pairs
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r += abs(acc(zi, pair.theta_col_i) - acc(zj, pair.theta_col_j))
        r += abs(acc(zi, pair.vm_col_i) - acc(zj, pair.vm_col_j))
    end
    return r
end

function _ac_resid_l1(ws::ACPowerWorkspaceALADIN, acc::F) where {F}
    r = 0.0
    @inbounds for pair in ws.interzonal_line_pairs
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r += abs(_ac_pf(zi, pair.line_local_i, acc) - _ac_pf(zj, pair.line_local_j, acc))
        r += abs(_ac_qf(zi, pair.line_local_i, acc) - _ac_qf(zj, pair.line_local_j, acc))
    end
    @inbounds for pair in ws.dummy_line_pairs
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r += abs(_ac_dpf(pair, zi, :i, acc) - _ac_dpf(pair, zj, :j, acc))
        r += abs(_ac_dqf(pair, zi, :i, acc) - _ac_dqf(pair, zj, :j, acc))
    end
    return r
end

function _ac_resid_l1(ws::ACHybridWorkspaceALADIN, acc::F) where {F}
    r = 0.0
    @inbounds for pair in ws.bus_pairs
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r += abs(acc(zi, pair.theta_col_i) - acc(zj, pair.theta_col_j))
        r += abs(acc(zi, pair.vm_col_i) - acc(zj, pair.vm_col_j))
    end
    @inbounds for pair in ws.interzonal_line_pairs
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r += abs(_ac_pf(zi, pair.line_local_i, acc) - _ac_pf(zj, pair.line_local_j, acc))
        r += abs(_ac_qf(zi, pair.line_local_i, acc) - _ac_qf(zj, pair.line_local_j, acc))
    end
    @inbounds for pair in ws.dummy_line_pairs
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r += abs(_ac_dpf(pair, zi, :i, acc) - _ac_dpf(pair, zj, :j, acc))
        r += abs(_ac_dqf(pair, zi, :i, acc) - _ac_dqf(pair, zj, :j, acc))
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
    return r
end

_ac_flow_prox_value(::ACStateWorkspaceALADIN) = 0.0
function _ac_flow_prox_value(ws::Union{ACPowerWorkspaceALADIN, ACHybridWorkspaceALADIN})
    ws.proximal_mode === :pflow || return 0.0
    hp = 0.5 * ws.rho[] * ws.sigma_pflow
    hq = 0.5 * ws.rho[] * ws.sigma_qflow
    s = 0.0
    @inbounds for pair in ws.interzonal_line_pairs
        for (zp, le) in ((pair.zone_i, pair.line_local_i), (pair.zone_j, pair.line_local_j))
            zone = ws.zones[zp]
            p_prev, q_prev = _ac_flow_prev_interzonal(zone, le)
            s += hp * (_ac_pf_at(zone, le) - p_prev)^2 + hq * (_ac_qf_at(zone, le) - q_prev)^2
        end
    end
    @inbounds for pair in ws.dummy_line_pairs
        for (zp, side) in ((pair.zone_i, :i), (pair.zone_j, :j))
            zone = ws.zones[zp]
            p_prev, q_prev = _ac_flow_prev_dummy(zone, pair, side)
            s += hp * (_dummy_ac_pf_at(pair, zone.y, side) - p_prev)^2 +
                 hq * (_dummy_ac_qf_at(pair, zone.y, side) - q_prev)^2
        end
    end
    return s
end

# Σ (ρ/2)‖yᵢ − xᵢ‖²_Σᵢ = value of the proximal quadratic (½ dᵀ(ρΣ)d) at d = y − x
function _ac_prox(ws::_ACPhase2Workspace)
    s = 0.0
    @inbounds for zone in ws.zones
        v2c = zone.var_to_col
        y = zone.y; x = zone.x_prev
        for (pair, a) in zone.quad_static.terms
            i = v2c[pair.a]; j = v2c[pair.b]
            s += a * (y[i] - x[i]) * (y[j] - x[j])
        end
    end
    return s + _ac_flow_prox_value(ws)
end

const _ALG3_PHI_LOCAL_VIOL = Ref(true)

function _ac_local_infeas_l1(ws::_ACPhase2Workspace, acc::F) where {F}
    total = 0.0
    for zone in ws.zones
        model = zone.local_model
        v2c = zone.var_to_col
        valfn = v -> acc(zone, v2c[v])
        for (Fc, S) in JuMP.list_of_constraint_types(model)
            Fc <: JuMP.VariableRef && continue
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

function _ac_max_local_dual(ws::_ACPhase2Workspace; include_bounds::Bool = false)
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
                d isa Real || continue
                ad = abs(d)
                ad > maxkap && (maxkap = ad)
            end
        end
        include_bounds || continue
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

function _alg3_step!(ws::_ACPhase2Workspace; gamma::Float64 = 1e-4)
    accx = (z, c) -> z.x_prev[c]
    accy = (z, c) -> z.y[c]
    accyd = (z, c) -> z.y[c] + z.delta[c]
    maxlam = 0.0
    for state in values(ws.vars)
        @inbounds for v in state.lambda
            av = abs(v)
            av > maxlam && (maxlam = av)
        end
    end
    use_viol = _ALG3_PHI_LOCAL_VIOL[]
    use_viol && (maxlam = max(maxlam, _ac_max_local_dual(ws)))
    if ws.lam_bar <= 0.0 || ws.lam_bar < 1.5 * maxlam
        ws.lam_bar = max(10.0 * maxlam, 1.0)
    end
    lam_bar = ws.lam_bar

    infeas = use_viol ?
        (acc -> lam_bar * _ac_local_infeas_l1(ws, acc)) : (acc -> 0.0)

    rhs = gamma * (_ac_prox(ws) + lam_bar * _ac_resid_l1(ws, accy))
    phi_x = _ac_cost_sum(ws, accx) + lam_bar * _ac_resid_l1(ws, accx) +
        infeas(accx)

    phi_a = _ac_cost_sum(ws, accyd) + lam_bar * _ac_resid_l1(ws, accyd) +
        infeas(accyd)
    if isfinite(phi_a) && phi_x - phi_a >= rhs
        ws.alpha1 = 1.0; ws.alpha2 = 1.0; ws.alpha3 = 1.0
        return nothing
    end

    phi_b = _ac_cost_sum(ws, accy) + lam_bar * _ac_resid_l1(ws, accy) +
        infeas(accy)
    if isfinite(phi_b) && phi_x - phi_b >= rhs
        ws.alpha1 = 1.0; ws.alpha2 = 0.0; ws.alpha3 = 0.0
        return nothing
    end

    return _branch_c_dual_linesearch!(ws)
end

# ============================================================
# HFD Algorithm 3 as printed — :alg3_original
# ============================================================

const _ALG3O_LAMBAR_HYST = 1.5      # raise λ̄ when it no longer majorises 1.5·m
const _ALG3O_LAMBAR_FACTOR = 10.0   # new λ̄ = max(10·m, _ALG3O_LAMBAR_MIN)
const _ALG3O_LAMBAR_MIN = 1.0

function _alg3_original_lam_bar!(ws::_ACPhase2Workspace)
    if ws.lam_bar_fixed
        maxlam = 0.0
        for state in values(ws.vars)
            @inbounds for v in state.lambda
                av = abs(v)
                av > maxlam && (maxlam = av)
            end
        end
        if maxlam > ws.lam_bar
            ws.hist.meta[:lam_bar_violations] =
                Int(get(ws.hist.meta, :lam_bar_violations, 0)) + 1
            @warn "ALADIN :alg3_original (AC) — |λ|∞ = $(maxlam) exceeds the " *
                  "fixed L1 penalty weight λ̄ = $(ws.lam_bar) (iter $(ws.iter)); " *
                  "Φ is no longer guaranteed exact." maxlog = 1
        end
        maxkap = _ac_max_local_dual(ws; include_bounds = true)
        if maxkap > ws.lam_bar
            ws.hist.meta[:kap_bar_violations] =
                Int(get(ws.hist.meta, :kap_bar_violations, 0)) + 1
            @warn "ALADIN :alg3_original (AC) — max|κ_local| = $(maxkap) exceeds the " *
                  "fixed local-penalty weight κ̄ = λ̄ = $(ws.lam_bar) (iter $(ws.iter)); " *
                  "the local term of Φ is no longer guaranteed exact." maxlog = 1
        end
        return ws.lam_bar
    end
    m = 0.0
    for state in values(ws.vars)
        @inbounds for k in eachindex(state.lambda)
            m = max(m, abs(state.lambda[k]), abs(state.nu[k]))
        end
    end
    m = max(m, _ac_max_local_dual(ws; include_bounds = true))
    if ws.lam_bar < _ALG3O_LAMBAR_HYST * m || ws.lam_bar <= 0.0
        ws.lam_bar = max(_ALG3O_LAMBAR_FACTOR * m, _ALG3O_LAMBAR_MIN, ws.lam_bar)
    end
    return ws.lam_bar
end

@inline _alg3o_dir(dirs::NamedTuple, ::_ACPhase2Workspace, key::Symbol) = dirs[key]
@inline function _alg3o_dir(dirs, ws::_ACPhase2Workspace, key::Symbol)
    idx = findfirst(==(key), collect(keys(ws.vars)))
    return dirs[idx]
end

function _dual_slope(ws::ACStateWorkspaceALADIN, dirs)
    d_th = _alg3o_dir(dirs, ws, :theta)
    d_vm = _alg3o_dir(dirs, ws, :v_mag)
    s = 0.0
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        s += d_th[k] * (zi.y[pair.theta_col_i] - zj.y[pair.theta_col_j])
        s += d_vm[k] * (zi.y[pair.vm_col_i] - zj.y[pair.vm_col_j])
    end
    return s
end

function _dual_slope(ws::ACPowerWorkspaceALADIN, dirs)
    accy = (z, c) -> z.y[c]
    d_p = _alg3o_dir(dirs, ws, :p_flow)
    d_q = _alg3o_dir(dirs, ws, :q_flow)
    d_pd = _alg3o_dir(dirs, ws, :p_dummy)
    d_qd = _alg3o_dir(dirs, ws, :q_dummy)
    s = 0.0
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        s += d_p[k] * (_ac_pf(zi, pair.line_local_i, accy) - _ac_pf(zj, pair.line_local_j, accy))
        s += d_q[k] * (_ac_qf(zi, pair.line_local_i, accy) - _ac_qf(zj, pair.line_local_j, accy))
    end
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        s += d_pd[k] * (_ac_dpf(pair, zi, :i, accy) - _ac_dpf(pair, zj, :j, accy))
        s += d_qd[k] * (_ac_dqf(pair, zi, :i, accy) - _ac_dqf(pair, zj, :j, accy))
    end
    return s
end

function _dual_slope(ws::ACHybridWorkspaceALADIN, dirs)
    accy = (z, c) -> z.y[c]
    d_th = _alg3o_dir(dirs, ws, :theta)
    d_vm = _alg3o_dir(dirs, ws, :v_mag)
    d_p = _alg3o_dir(dirs, ws, :p_flow)
    d_q = _alg3o_dir(dirs, ws, :q_flow)
    d_pi = _alg3o_dir(dirs, ws, :p_inj)
    d_qi = _alg3o_dir(dirs, ws, :q_inj)
    d_pd = _alg3o_dir(dirs, ws, :p_dummy)
    d_qd = _alg3o_dir(dirs, ws, :q_dummy)
    s = 0.0
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        s += d_th[k] * (zi.y[pair.theta_col_i] - zj.y[pair.theta_col_j])
        s += d_vm[k] * (zi.y[pair.vm_col_i] - zj.y[pair.vm_col_j])
    end
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        s += d_p[k] * (_ac_pf(zi, pair.line_local_i, accy) - _ac_pf(zj, pair.line_local_j, accy))
        s += d_q[k] * (_ac_qf(zi, pair.line_local_i, accy) - _ac_qf(zj, pair.line_local_j, accy))
    end
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        s += d_pd[k] * (_ac_dpf(pair, zi, :i, accy) - _ac_dpf(pair, zj, :j, accy))
        s += d_qd[k] * (_ac_dqf(pair, zi, :i, accy) - _ac_dqf(pair, zj, :j, accy))
    end
    @inbounds for (k, pair) in enumerate(ws.split_bus_pairs)
        dp = 0.0; dq = 0.0
        for c in eachindex(pair.zones)
            zone = ws.zones[pair.zones[c]]
            dp += zone.y[pair.p_inj_cols[c]]
            dq += zone.y[pair.q_inj_cols[c]]
        end
        s += d_pi[k] * dp + d_qi[k] * dq
    end
    return s
end

function _dual_residual(ws::ACStateWorkspaceALADIN)
    n = length(ws.bus_pairs)
    r_th = Vector{Float64}(undef, n)
    r_vm = Vector{Float64}(undef, n)
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r_th[k] = zi.y[pair.theta_col_i] - zj.y[pair.theta_col_j]
        r_vm[k] = zi.y[pair.vm_col_i] - zj.y[pair.vm_col_j]
    end
    return NamedTuple{keys(ws.vars)}((theta = r_th, v_mag = r_vm))
end

function _dual_residual(ws::ACPowerWorkspaceALADIN)
    accy = (z, c) -> z.y[c]
    n_e = length(ws.interzonal_line_pairs)
    n_d = length(ws.dummy_line_pairs)
    r_p = Vector{Float64}(undef, n_e)
    r_q = Vector{Float64}(undef, n_e)
    r_pd = Vector{Float64}(undef, n_d)
    r_qd = Vector{Float64}(undef, n_d)
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r_p[k] = _ac_pf(zi, pair.line_local_i, accy) - _ac_pf(zj, pair.line_local_j, accy)
        r_q[k] = _ac_qf(zi, pair.line_local_i, accy) - _ac_qf(zj, pair.line_local_j, accy)
    end
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r_pd[k] = _ac_dpf(pair, zi, :i, accy) - _ac_dpf(pair, zj, :j, accy)
        r_qd[k] = _ac_dqf(pair, zi, :i, accy) - _ac_dqf(pair, zj, :j, accy)
    end
    return NamedTuple{keys(ws.vars)}(
        (p_flow = r_p, q_flow = r_q, p_dummy = r_pd, q_dummy = r_qd),
    )
end

function _dual_residual(ws::ACHybridWorkspaceALADIN)
    accy = (z, c) -> z.y[c]
    n_b = length(ws.bus_pairs)
    n_e = length(ws.interzonal_line_pairs)
    n_d = length(ws.dummy_line_pairs)
    n_s = length(ws.split_bus_pairs)
    r_th = Vector{Float64}(undef, n_b)
    r_vm = Vector{Float64}(undef, n_b)
    r_p = Vector{Float64}(undef, n_e)
    r_q = Vector{Float64}(undef, n_e)
    r_pd = Vector{Float64}(undef, n_d)
    r_qd = Vector{Float64}(undef, n_d)
    r_pi = Vector{Float64}(undef, n_s)
    r_qi = Vector{Float64}(undef, n_s)
    @inbounds for (k, pair) in enumerate(ws.bus_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r_th[k] = zi.y[pair.theta_col_i] - zj.y[pair.theta_col_j]
        r_vm[k] = zi.y[pair.vm_col_i] - zj.y[pair.vm_col_j]
    end
    @inbounds for (k, pair) in enumerate(ws.interzonal_line_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r_p[k] = _ac_pf(zi, pair.line_local_i, accy) - _ac_pf(zj, pair.line_local_j, accy)
        r_q[k] = _ac_qf(zi, pair.line_local_i, accy) - _ac_qf(zj, pair.line_local_j, accy)
    end
    @inbounds for (k, pair) in enumerate(ws.dummy_line_pairs)
        zi = ws.zones[pair.zone_i]; zj = ws.zones[pair.zone_j]
        r_pd[k] = _ac_dpf(pair, zi, :i, accy) - _ac_dpf(pair, zj, :j, accy)
        r_qd[k] = _ac_dqf(pair, zi, :i, accy) - _ac_dqf(pair, zj, :j, accy)
    end
    @inbounds for (k, pair) in enumerate(ws.split_bus_pairs)
        dp = 0.0; dq = 0.0
        for c in eachindex(pair.zones)
            zone = ws.zones[pair.zones[c]]
            dp += zone.y[pair.p_inj_cols[c]]
            dq += zone.y[pair.q_inj_cols[c]]
        end
        r_pi[k] = dp
        r_qi[k] = dq
    end
    return NamedTuple{keys(ws.vars)}((
        theta = r_th, v_mag = r_vm, p_flow = r_p, q_flow = r_q,
        p_inj = r_pi, q_inj = r_qi, p_dummy = r_pd, q_dummy = r_qd,
    ))
end

@inline function _alg3o_record!(ws::_ACPhase2Workspace, branch::Symbol, status_c::Symbol)
    log = get!(ws.hist.meta, :alg3_original,
        Vector{Tuple{Symbol, Float64, Float64, Symbol}}())
    push!(log::Vector{Tuple{Symbol, Float64, Float64, Symbol}},
        (branch, ws.alpha3, ws.lam_bar, status_c))
    return nothing
end

function _alg3_original_step!(ws::_ACPhase2Workspace; gamma::Float64 = 1e-4)
    accx = (z, c) -> z.x_prev[c]
    accy = (z, c) -> z.y[c]
    accyd = (z, c) -> z.y[c] + z.delta[c]

    lam_bar = _alg3_original_lam_bar!(ws)

    # Φ(x) = Σ fᵢ(xᵢ) + λ̄‖Σ Aᵢxᵢ − b‖₁ + κ̄ Σ max{0, (hᵢ(xᵢ))ⱼ},  κ̄ = λ̄
    Phi(acc) = _ac_cost_sum(ws, acc) + lam_bar * _ac_resid_l1(ws, acc) +
        lam_bar * _ac_local_infeas_l1(ws, acc)

    # (6.2) right-hand side: γ ( Σ ρ/2‖yᵢ − xᵢ‖²_Σᵢ + λ̄ ‖Σ Aᵢ yᵢ − b‖₁ )
    rhs = gamma * (_ac_prox(ws) + lam_bar * _ac_resid_l1(ws, accy))
    phi_x = Phi(accx)

    # (a) α₁ = α₂ = α₃ = 1, x⁺ = x + (y − x) + Δy = y + Δy
    phi_a = Phi(accyd)
    if isfinite(phi_a) && phi_x - phi_a >= rhs
        ws.alpha1 = 1.0; ws.alpha2 = 1.0; ws.alpha3 = 1.0
        ws.alg3_c_status = :none
        _alg3o_record!(ws, :a, :none)
        return nothing
    end

    # (b) α₁ = 1, α₂ = α₃ = 0, x⁺ = y
    phi_b = Phi(accy)
    if isfinite(phi_b) && phi_x - phi_b >= rhs
        ws.alpha1 = 1.0; ws.alpha2 = 0.0; ws.alpha3 = 0.0
        ws.alg3_c_status = :none
        _alg3o_record!(ws, :b, :none)
        return nothing
    end

    # (c) α₁ = α₂ = 0, α₃ from the Wolfe ascent search on V_ρ(x, λ + α(λ_QP − λ))
    status = _branch_c_dual_linesearch_wolfe!(ws)
    ws.alg3_c_status = status
    _alg3o_record!(ws, :c, status)
    return nothing
end

# ============================================================
# Armijo backtracking step
# ============================================================

function _line_search_step!(
    ws::_ACPhase2Workspace;
    c1::Float64 = 1e-4,
    max_backtrack::Int = 25,
    alpha_min::Float64 = 1e-3,
)
    ws.line_search || return nothing
    ws.globalization === :alg3 && return _alg3_step!(ws)
    ws.globalization === :alg3_original && return _alg3_original_step!(ws)
    # :full_step — HFD Alg.2 step 5 undamped step (α₁=α₂=1, α₃=alpha3_cfg)
    ws.globalization === :full_step && return nothing

    phi_0 = _ac_merit(ws, 0.0)
    dphi_0 = _ac_dphi_0(ws)

    if !isfinite(dphi_0) || dphi_0 >= 0.0
        return _branch_c_dual_linesearch!(ws)
    end

    alpha = ws.alpha2_init
    @inbounds for _ in 1:max_backtrack
        phi_alpha = _ac_merit(ws, alpha)
        if isfinite(phi_alpha) && phi_alpha <= phi_0 + c1 * alpha * dphi_0
            ws.alpha2 = alpha
            return nothing
        end
        alpha *= 0.5
        if alpha < alpha_min
            break
        end
    end
    return _branch_c_dual_linesearch!(ws)
end
