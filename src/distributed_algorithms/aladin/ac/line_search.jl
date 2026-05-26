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
    return g * (tau * vf)^2 -
           tau * vf * vt * (g * cos(delta_theta) + b * sin(delta_theta))
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
    return -(b + b_sh / 2) * (tau * vf)^2 -
           tau * vf * vt * (-b * cos(delta_theta) + g * sin(delta_theta))
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
# Armijo backtracking step
# ============================================================

function _line_search_step!(
    ws::_ACPhase2Workspace;
    c1::Float64 = 1e-4,
    max_backtrack::Int = 25,
    alpha_min::Float64 = 1e-3,
)
    ws.line_search || return nothing

    phi_0 = _ac_merit(ws, 0.0)
    dphi_0 = _ac_dphi_0(ws)

    if !isfinite(dphi_0) || dphi_0 >= 0.0
        ws.alpha2 = alpha_min
        return nothing
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
    ws.alpha2 = alpha_min
    return nothing
end
