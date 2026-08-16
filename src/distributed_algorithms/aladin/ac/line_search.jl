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
        # half-A: dummy bus is its TO end -> TO-side flow, matching the consensus /
        # coordinator / residual. See _ac_flow_at_to.
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
        # half-A: derivatives of the TO-side flow, so the step search moves along the same
        # quantity the coupling constrains (see _ac_flow_at_to).
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
#
# Generic point evaluators parametrised by an accessor `acc(zone, col)::Float64`
# giving the value of column `col` at the trial primal point. The three points of
# Algorithm 3 are reached by:
#   x  (current) → acc = (z,c) -> z.x_prev[c]
#   y            → acc = (z,c) -> z.y[c]                (step (b): x⁺ = y)
#   y + Δy       → acc = (z,c) -> z.y[c] + z.delta[c]   (step (a): full step)
# The flow formulas mirror the *_at_trial evaluators exactly; only the value
# source differs.

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

# Σ (ρ/2)‖yᵢ − xᵢ‖²_Σᵢ = value of the proximal quadratic (½ dᵀ(ρΣ)d) at d = y − x.
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
    return s
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
    if ws.lam_bar <= 0.0 || ws.lam_bar < 1.5 * maxlam
        ws.lam_bar = max(10.0 * maxlam, 1.0)
    end
    lam_bar = ws.lam_bar

    rhs = gamma * (_ac_prox(ws) + lam_bar * _ac_resid_l1(ws, accy))
    phi_x = _ac_cost_sum(ws, accx) + lam_bar * _ac_resid_l1(ws, accx)

    phi_a = _ac_cost_sum(ws, accyd) + lam_bar * _ac_resid_l1(ws, accyd)
    if isfinite(phi_a) && phi_x - phi_a >= rhs
        ws.alpha1 = 1.0; ws.alpha2 = 1.0; ws.alpha3 = 1.0
        return nothing
    end

    phi_b = _ac_cost_sum(ws, accy) + lam_bar * _ac_resid_l1(ws, accy)
    if isfinite(phi_b) && phi_x - phi_b >= rhs
        ws.alpha1 = 1.0; ws.alpha2 = 0.0; ws.alpha3 = 0.0
        return nothing
    end

    return _branch_c_dual_linesearch!(ws)
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
