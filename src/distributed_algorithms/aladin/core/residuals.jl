# ============================================================
# ALADIN — Residuals (Boyd 2-norm or HFD ∞-norm)
# ============================================================

@inline _final_norm(crit::Symbol, sq::Float64, mx::Float64) =
    crit === :hfd ? mx : sqrt(sq)

# ============================================================
# Shared primal step / lambda step / dual residual
# ============================================================

@inline function _progress_norms(zone::AbstractAladinZoneData)
    s_sq = 0.0
    s_inf = 0.0
    y = zone.y
    x = zone.x_prev
    @inbounds for k in 1:zone.n_vars
        d = abs(y[k] - x[k])
        s_sq += d * d
        s_inf = max(s_inf, d)
    end
    return s_sq, s_inf
end

function _compute_r_step_primal!(ws::AbstractALADINWorkspace)
    s_sq = 0.0
    s_inf = 0.0
    @inbounds for zone in ws.zones
        zsq, zinf = _progress_norms(zone)
        s_sq += zsq
        s_inf = max(s_inf, zinf)
    end
    ws.r_step_primal = _final_norm(ws.stopping_criterion, s_sq, s_inf)
    return nothing
end

@inline function _lambda_step_norms(state::LambdaState)
    s_sq = 0.0
    s_inf = 0.0
    @inbounds for k in eachindex(state.lambda)
        d = abs(state.lambda[k] - state.lambda_prev[k])
        s_sq += d * d
        s_inf = max(s_inf, d)
    end
    return s_sq, s_inf
end

function _compute_lambda_step!(ws::AbstractALADINWorkspace)
    s_sq = 0.0
    s_inf = 0.0
    for state in values(ws.vars)
        ssq, sinf = _lambda_step_norms(state)
        s_sq += ssq
        s_inf = max(s_inf, sinf)
    end
    ws.lambda_step = _final_norm(ws.stopping_criterion, s_sq, s_inf)
    return nothing
end

@inline function _compute_r_dual!(ws::AbstractALADINWorkspace)
    ws.r_dual = ws.mu[] * ws.lambda_step
    return nothing
end

# ============================================================
# Per-block consensus accumulators
# Each returns (sum_of_squares, max_abs)
# ============================================================

@inline function _r_pri_bus_dc(ws)
    sq = 0.0
    mx = 0.0
    @inbounds for pair in ws.bus_pairs
        d = abs(ws.zones[pair.zone_i].y[pair.col_i] -
                ws.zones[pair.zone_j].y[pair.col_j])
        sq += d * d
        mx = max(mx, d)
    end
    return sq, mx
end

@inline function _r_pri_edge_dc(ws)
    sq = 0.0
    mx = 0.0
    @inbounds for pair in ws.interzonal_line_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        li = pair.line_local_i
        lj = pair.line_local_j
        pfi = zi.edge_b[li] * (zi.y[zi.edge_from_cols[li]] - zi.y[zi.edge_to_cols[li]])
        pfj = zj.edge_b[lj] * (zj.y[zj.edge_from_cols[lj]] - zj.y[zj.edge_to_cols[lj]])
        d = abs(pfi - pfj)
        sq += d * d
        mx = max(mx, d)
    end
    return sq, mx
end

@inline function _r_pri_dummy_dc(ws)
    sq = 0.0
    mx = 0.0
    @inbounds for pair in ws.dummy_line_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pfi = pair.b_i * (zi.y[pair.from_col_i] - zi.y[pair.to_col_i])
        pfj = pair.b_j * (zj.y[pair.from_col_j] - zj.y[pair.to_col_j])
        d = abs(pfi - pfj)
        sq += d * d
        mx = max(mx, d)
    end
    return sq, mx
end

@inline function _r_pri_split_dc(ws)
    sq = 0.0
    mx = 0.0
    @inbounds for pair in ws.split_bus_pairs
        d = abs(
            ws.zones[pair.zone_i].y[pair.p_inj_col_i] +
            ws.zones[pair.zone_j].y[pair.p_inj_col_j],
        )
        sq += d * d
        mx = max(mx, d)
    end
    return sq, mx
end

# ============================================================
# DC residuals
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function compute_residuals!(ws::DCStateWorkspaceALADIN)
    sq, mx = _r_pri_bus_dc(ws)
    ws.r_primal = _final_norm(ws.stopping_criterion, sq, mx)
    _compute_r_step_primal!(ws)
    _compute_lambda_step!(ws)
    _compute_r_dual!(ws)
    return nothing
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function compute_residuals!(ws::DCPowerWorkspaceALADIN)
    sq_e, mx_e = _r_pri_edge_dc(ws)
    sq_d, mx_d = _r_pri_dummy_dc(ws)
    ws.r_primal = _final_norm(
        ws.stopping_criterion, sq_e + sq_d, max(mx_e, mx_d),
    )
    _compute_r_step_primal!(ws)
    _compute_lambda_step!(ws)
    _compute_r_dual!(ws)
    return nothing
end

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

function compute_residuals!(ws::DCHybridWorkspaceALADIN)
    sq_b, mx_b = _r_pri_bus_dc(ws)
    sq_e, mx_e = _r_pri_edge_dc(ws)
    sq_s, mx_s = _r_pri_split_dc(ws)
    sq_d, mx_d = _r_pri_dummy_dc(ws)

    crit = ws.stopping_criterion
    ws.r_primal_bus = _final_norm(crit, sq_b, mx_b)
    ws.r_primal_edge = _final_norm(crit, sq_e, mx_e)
    ws.r_primal = _final_norm(
        crit, sq_b + sq_e + sq_s + sq_d, max(mx_b, mx_e, mx_s, mx_d),
    )

    _compute_r_step_primal!(ws)
    ws.r_step_primal_bus = ws.r_step_primal       # same scalar, kept for diagnostics
    ws.r_step_primal_edge = ws.r_step_primal

    _compute_lambda_step!(ws)
    _compute_r_dual!(ws)
    return nothing
end

# ============================================================
# AC consensus accumulators
# ============================================================

@inline function _r_pri_bus_ac(ws)
    sq = 0.0
    mx = 0.0
    @inbounds for pair in ws.bus_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        dth = abs(zi.y[pair.theta_col_i] - zj.y[pair.theta_col_j])
        dvm = abs(zi.y[pair.vm_col_i] - zj.y[pair.vm_col_j])
        sq += dth * dth + dvm * dvm
        mx = max(mx, dth, dvm)
    end
    return sq, mx
end

@inline function _r_pri_edge_pq_ac(ws)
    sq = 0.0
    mx = 0.0
    @inbounds for pair in ws.interzonal_line_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        dpf = abs(_ac_pf_at(zi, pair.line_local_i) - _ac_pf_at(zj, pair.line_local_j))
        dqf = abs(_ac_qf_at(zi, pair.line_local_i) - _ac_qf_at(zj, pair.line_local_j))
        sq += dpf * dpf + dqf * dqf
        mx = max(mx, dpf, dqf)
    end
    return sq, mx
end

@inline function _r_pri_split_ac(ws)
    sq = 0.0
    mx = 0.0
    @inbounds for pair in ws.split_bus_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        dp = abs(zi.y[pair.p_inj_col_i] + zj.y[pair.p_inj_col_j])
        dq = abs(zi.y[pair.q_inj_col_i] + zj.y[pair.q_inj_col_j])
        sq += dp * dp + dq * dq
        mx = max(mx, dp, dq)
    end
    return sq, mx
end

# ------------------------------------------------------------
# AC dummy-line flow helpers
# ------------------------------------------------------------

@inline function _dummy_ac_pf_at(pair::DummyLinePairAC, y::Vector{Float64}, side::Symbol)
    if side === :i
        theta_f = y[pair.from_theta_col_i]
        theta_t = y[pair.to_theta_col_i]
        vf = y[pair.from_vm_col_i]
        vt = y[pair.to_vm_col_i]
        g = pair.g_i
        b = pair.b_i
        tau = pair.inv_tap_i
        tap_angle = pair.tap_angle_i
    else
        theta_f = y[pair.from_theta_col_j]
        theta_t = y[pair.to_theta_col_j]
        vf = y[pair.from_vm_col_j]
        vt = y[pair.to_vm_col_j]
        g = pair.g_j
        b = pair.b_j
        tau = pair.inv_tap_j
        tap_angle = pair.tap_angle_j
    end
    delta_theta = theta_f - theta_t - tap_angle
    return g * (tau * vf)^2 -
           tau * vf * vt * (g * cos(delta_theta) + b * sin(delta_theta))
end

@inline function _dummy_ac_qf_at(pair::DummyLinePairAC, y::Vector{Float64}, side::Symbol)
    if side === :i
        theta_f = y[pair.from_theta_col_i]
        theta_t = y[pair.to_theta_col_i]
        vf = y[pair.from_vm_col_i]
        vt = y[pair.to_vm_col_i]
        g = pair.g_i
        b = pair.b_i
        b_sh = pair.b_sh_i
        tau = pair.inv_tap_i
        tap_angle = pair.tap_angle_i
    else
        theta_f = y[pair.from_theta_col_j]
        theta_t = y[pair.to_theta_col_j]
        vf = y[pair.from_vm_col_j]
        vt = y[pair.to_vm_col_j]
        g = pair.g_j
        b = pair.b_j
        b_sh = pair.b_sh_j
        tau = pair.inv_tap_j
        tap_angle = pair.tap_angle_j
    end
    delta_theta = theta_f - theta_t - tap_angle
    return -(b + b_sh / 2) * (tau * vf)^2 -
           tau * vf * vt * (-b * cos(delta_theta) + g * sin(delta_theta))
end

@inline function _r_pri_dummy_pq_ac(ws)
    sq = 0.0
    mx = 0.0
    @inbounds for pair in ws.dummy_line_pairs
        yi = ws.zones[pair.zone_i].y
        yj = ws.zones[pair.zone_j].y
        dpf = abs(_dummy_ac_pf_at(pair, yi, :i) - _dummy_ac_pf_at(pair, yj, :j))
        dqf = abs(_dummy_ac_qf_at(pair, yi, :i) - _dummy_ac_qf_at(pair, yj, :j))
        sq += dpf * dpf + dqf * dqf
        mx = max(mx, dpf, dqf)
    end
    return sq, mx
end

# ============================================================
# AC residuals
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function compute_residuals!(ws::ACStateWorkspaceALADIN)
    sq, mx = _r_pri_bus_ac(ws)
    ws.r_primal = _final_norm(ws.stopping_criterion, sq, mx)
    _compute_r_step_primal!(ws)
    _compute_lambda_step!(ws)
    _compute_r_dual!(ws)
    return nothing
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function compute_residuals!(ws::ACPowerWorkspaceALADIN)
    sq_e, mx_e = _r_pri_edge_pq_ac(ws)
    sq_d, mx_d = _r_pri_dummy_pq_ac(ws)
    ws.r_primal = _final_norm(
        ws.stopping_criterion, sq_e + sq_d, max(mx_e, mx_d),
    )
    _compute_r_step_primal!(ws)
    _compute_lambda_step!(ws)
    _compute_r_dual!(ws)
    return nothing
end

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

function compute_residuals!(ws::ACHybridWorkspaceALADIN)
    sq_b, mx_b = _r_pri_bus_ac(ws)
    sq_e, mx_e = _r_pri_edge_pq_ac(ws)
    sq_s, mx_s = _r_pri_split_ac(ws)
    sq_d, mx_d = _r_pri_dummy_pq_ac(ws)

    ws.r_primal = _final_norm(
        ws.stopping_criterion,
        sq_b + sq_e + sq_s + sq_d,
        max(mx_b, mx_e, mx_s, mx_d),
    )
    _compute_r_step_primal!(ws)
    _compute_lambda_step!(ws)
    _compute_r_dual!(ws)
    return nothing
end

# ============================================================
# SOCWR consensus accumulators
# ============================================================

@inline function _r_pri_bus_socwr(ws)
    sq = 0.0
    mx = 0.0
    @inbounds for pair in ws.bus_pairs
        d = abs(
            ws.zones[pair.zone_i].y[pair.v2_col_i] -
            ws.zones[pair.zone_j].y[pair.v2_col_j],
        )
        sq += d * d
        mx = max(mx, d)
    end
    return sq, mx
end

@inline function _r_pri_rt_socwr(ws)
    sq = 0.0
    mx = 0.0
    @inbounds for pair in ws.interzonal_line_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        dr = abs(
            zi.y[zi.edge_r_cols[pair.line_local_i]] -
            zj.y[zj.edge_r_cols[pair.line_local_j]],
        )
        dt = abs(
            zi.y[zi.edge_t_cols[pair.line_local_i]] -
            zj.y[zj.edge_t_cols[pair.line_local_j]],
        )
        sq += dr * dr + dt * dt
        mx = max(mx, dr, dt)
    end
    return sq, mx
end

@inline function _r_pri_pq_socwr(ws)
    sq = 0.0
    mx = 0.0
    @inbounds for pair in ws.interzonal_line_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        dpf = abs(
            _socwr_flow_value_at_y(zi, pair.line_local_i, :p) -
            _socwr_flow_value_at_y(zj, pair.line_local_j, :p),
        )
        dqf = abs(
            _socwr_flow_value_at_y(zi, pair.line_local_i, :q) -
            _socwr_flow_value_at_y(zj, pair.line_local_j, :q),
        )
        sq += dpf * dpf + dqf * dqf
        mx = max(mx, dpf, dqf)
    end
    return sq, mx
end

@inline function _r_pri_split_socwr(ws)
    sq = 0.0
    mx = 0.0
    @inbounds for pair in ws.split_bus_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        dp = abs(zi.y[pair.p_inj_col_i] + zj.y[pair.p_inj_col_j])
        dq = abs(zi.y[pair.q_inj_col_i] + zj.y[pair.q_inj_col_j])
        sq += dp * dp + dq * dq
        mx = max(mx, dp, dq)
    end
    return sq, mx
end

# ============================================================
# SOCWR residuals
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function compute_residuals!(ws::SOCWRStateWorkspaceALADIN)
    sq_b, mx_b = _r_pri_bus_socwr(ws)
    sq_rt, mx_rt = _r_pri_rt_socwr(ws)
    ws.r_primal = _final_norm(
        ws.stopping_criterion, sq_b + sq_rt, max(mx_b, mx_rt),
    )
    _compute_r_step_primal!(ws)
    _compute_lambda_step!(ws)
    _compute_r_dual!(ws)
    return nothing
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function compute_residuals!(ws::SOCWRPowerWorkspaceALADIN)
    sq, mx = _r_pri_pq_socwr(ws)
    ws.r_primal = _final_norm(ws.stopping_criterion, sq, mx)
    _compute_r_step_primal!(ws)
    _compute_lambda_step!(ws)
    _compute_r_dual!(ws)
    return nothing
end

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

function compute_residuals!(ws::SOCWRHybridWorkspaceALADIN)
    sq_b, mx_b = _r_pri_bus_socwr(ws)
    sq_rt, mx_rt = _r_pri_rt_socwr(ws)
    sq_pq, mx_pq = _r_pri_pq_socwr(ws)
    sq_s, mx_s = _r_pri_split_socwr(ws)
    ws.r_primal = _final_norm(
        ws.stopping_criterion,
        sq_b + sq_rt + sq_pq + sq_s,
        max(mx_b, mx_rt, mx_pq, mx_s),
    )
    _compute_r_step_primal!(ws)
    _compute_lambda_step!(ws)
    _compute_r_dual!(ws)
    return nothing
end
