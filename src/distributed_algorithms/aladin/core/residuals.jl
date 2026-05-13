function _progress_sq(zone::AbstractAladinZoneData)::Float64
    s = 0.0
    y = zone.y
    x = zone.x_prev
    @inbounds for k in 1:zone.n_vars
        d = y[k] - x[k]
        s += d * d
    end
    return s
end

function _compute_r_progress!(ws::AbstractALADINWorkspace)
    s = 0.0
    @inbounds for zone in ws.zones
        s += _progress_sq(zone)
    end
    ws.r_progress = sqrt(s)
    return nothing
end

@inline function _lambda_step_sq(state::LambdaState)::Float64
    s = 0.0
    @inbounds for k in eachindex(state.lambda)
        d = state.lambda[k] - state.lambda_prev[k]
        s += d * d
    end
    return s
end

function _compute_lambda_step!(ws::AbstractALADINWorkspace)
    s = 0.0
    for state in values(ws.vars)
        s += _lambda_step_sq(state)
    end
    ws.lambda_step = sqrt(s)
    return nothing
end

# ---- State ----

function compute_residuals!(ws::DCStateWorkspaceALADIN)
    r_pri_sq = 0.0
    @inbounds for pair in ws.bus_pairs
        yi = ws.zones[pair.zone_i].y[pair.col_i]
        yj = ws.zones[pair.zone_j].y[pair.col_j]
        r_pri_sq += (yi - yj)^2
    end
    ws.r_primal = sqrt(r_pri_sq)

    _compute_r_progress!(ws)
    _compute_lambda_step!(ws)
    return nothing
end

# ---- Power ----

@inline function _dummy_pf_diff_sq(ws)
    s = 0.0
    @inbounds for pair in ws.dummy_line_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pfi = pair.b_i * (zi.y[pair.from_col_i] - zi.y[pair.to_col_i])
        pfj = pair.b_j * (zj.y[pair.from_col_j] - zj.y[pair.to_col_j])
        s += (pfi - pfj)^2
    end
    return s
end

function compute_residuals!(ws::DCPowerWorkspaceALADIN)
    r_pri_sq = 0.0
    @inbounds for pair in ws.interzonal_line_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        li = pair.line_local_i
        lj = pair.line_local_j
        pfi = zi.edge_b[li] * (zi.y[zi.edge_from_cols[li]] - zi.y[zi.edge_to_cols[li]])
        pfj = zj.edge_b[lj] * (zj.y[zj.edge_from_cols[lj]] - zj.y[zj.edge_to_cols[lj]])
        r_pri_sq += (pfi - pfj)^2
    end
    r_pri_sq += _dummy_pf_diff_sq(ws)
    ws.r_primal = sqrt(r_pri_sq)

    _compute_r_progress!(ws)
    _compute_lambda_step!(ws)
    return nothing
end

# ---- DC Hybrid ----

function compute_residuals!(ws::DCHybridWorkspaceALADIN)
    # Bus consensus block
    r_pri_bus_sq = 0.0
    @inbounds for pair in ws.bus_pairs
        yi = ws.zones[pair.zone_i].y[pair.col_i]
        yj = ws.zones[pair.zone_j].y[pair.col_j]
        r_pri_bus_sq += (yi - yj)^2
    end
    ws.r_primal_bus = sqrt(r_pri_bus_sq)

    # Edge consensus block
    r_pri_edge_sq = 0.0
    @inbounds for pair in ws.interzonal_line_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        li = pair.line_local_i
        lj = pair.line_local_j
        pfi = zi.edge_b[li] * (zi.y[zi.edge_from_cols[li]] - zi.y[zi.edge_to_cols[li]])
        pfj = zj.edge_b[lj] * (zj.y[zj.edge_from_cols[lj]] - zj.y[zj.edge_to_cols[lj]])
        r_pri_edge_sq += (pfi - pfj)^2
    end
    ws.r_primal_edge = sqrt(r_pri_edge_sq)

    r_pri_split_sq = 0.0
    @inbounds for pair in ws.split_bus_pairs
        yi = ws.zones[pair.zone_i].y[pair.p_inj_col_i]
        yj = ws.zones[pair.zone_j].y[pair.p_inj_col_j]
        r_pri_split_sq += (yi + yj)^2
    end

    r_pri_dummy_sq = _dummy_pf_diff_sq(ws)

    ws.r_primal = sqrt(
        r_pri_bus_sq + r_pri_edge_sq + r_pri_split_sq + r_pri_dummy_sq,
    )

    _compute_r_progress!(ws)
    ws.r_progress_bus = ws.r_progress       # same scalar, kept for diagnostics
    ws.r_progress_edge = ws.r_progress

    _compute_lambda_step!(ws)
    return nothing
end

# ------------------------------------------------------------
# AC residuals
# ------------------------------------------------------------

function compute_residuals!(ws::ACStateWorkspaceALADIN)
    r_pri_sq = 0.0
    @inbounds for pair in ws.bus_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        theta_i = zi.y[pair.theta_col_i]
        theta_j = zj.y[pair.theta_col_j]
        vi = zi.y[pair.vm_col_i]
        vj = zj.y[pair.vm_col_j]
        r_pri_sq += (theta_i - theta_j)^2 + (vi - vj)^2
    end
    ws.r_primal = sqrt(r_pri_sq)

    _compute_r_progress!(ws)
    _compute_lambda_step!(ws)
    return nothing
end

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

@inline function _r_pri_dummy_pq_sq_ac(ws)
    s = 0.0
    @inbounds for pair in ws.dummy_line_pairs
        yi = ws.zones[pair.zone_i].y
        yj = ws.zones[pair.zone_j].y
        pfi = _dummy_ac_pf_at(pair, yi, :i)
        pfj = _dummy_ac_pf_at(pair, yj, :j)
        qfi = _dummy_ac_qf_at(pair, yi, :i)
        qfj = _dummy_ac_qf_at(pair, yj, :j)
        s += (pfi - pfj)^2 + (qfi - qfj)^2
    end
    return s
end

function compute_residuals!(ws::ACPowerWorkspaceALADIN)
    r_pri_sq = 0.0
    @inbounds for pair in ws.interzonal_line_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pfi = _ac_pf_at(zi, pair.line_local_i)
        pfj = _ac_pf_at(zj, pair.line_local_j)
        qfi = _ac_qf_at(zi, pair.line_local_i)
        qfj = _ac_qf_at(zj, pair.line_local_j)
        r_pri_sq += (pfi - pfj)^2 + (qfi - qfj)^2
    end
    r_pri_sq += _r_pri_dummy_pq_sq_ac(ws)
    ws.r_primal = sqrt(r_pri_sq)

    _compute_r_progress!(ws)
    _compute_lambda_step!(ws)
    return nothing
end

function compute_residuals!(ws::ACHybridWorkspaceALADIN)
    r_pri_sq = 0.0
    @inbounds for pair in ws.bus_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        theta_i = zi.y[pair.theta_col_i]
        theta_j = zj.y[pair.theta_col_j]
        vi = zi.y[pair.vm_col_i]
        vj = zj.y[pair.vm_col_j]
        r_pri_sq += (theta_i - theta_j)^2 + (vi - vj)^2
    end
    @inbounds for pair in ws.interzonal_line_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pfi = _ac_pf_at(zi, pair.line_local_i)
        pfj = _ac_pf_at(zj, pair.line_local_j)
        qfi = _ac_qf_at(zi, pair.line_local_i)
        qfj = _ac_qf_at(zj, pair.line_local_j)
        r_pri_sq += (pfi - pfj)^2 + (qfi - qfj)^2
    end
    @inbounds for pair in ws.split_bus_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = zi.y[pair.p_inj_col_i]
        pj_val = zj.y[pair.p_inj_col_j]
        qi_val = zi.y[pair.q_inj_col_i]
        qj_val = zj.y[pair.q_inj_col_j]
        r_pri_sq += (pi_val + pj_val)^2 + (qi_val + qj_val)^2
    end
    r_pri_sq += _r_pri_dummy_pq_sq_ac(ws)
    ws.r_primal = sqrt(r_pri_sq)

    _compute_r_progress!(ws)
    _compute_lambda_step!(ws)
    return nothing
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
    @inbounds for pair in ws.interzonal_line_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        ri = zi.y[zi.edge_r_cols[pair.line_local_i]]
        rj = zj.y[zj.edge_r_cols[pair.line_local_j]]
        ti = zi.y[zi.edge_t_cols[pair.line_local_i]]
        tj = zj.y[zj.edge_t_cols[pair.line_local_j]]
        s += (ri - rj)^2 + (ti - tj)^2
    end
    return s
end

@inline function _r_pri_pq_sq_socwr(ws)
    s = 0.0
    @inbounds for pair in ws.interzonal_line_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pfi = _socwr_flow_value_at_y(zi, pair.line_local_i, :p)
        pfj = _socwr_flow_value_at_y(zj, pair.line_local_j, :p)
        qfi = _socwr_flow_value_at_y(zi, pair.line_local_i, :q)
        qfj = _socwr_flow_value_at_y(zj, pair.line_local_j, :q)
        s += (pfi - pfj)^2 + (qfi - qfj)^2
    end
    return s
end

@inline function _r_pri_split_pq_sq_socwr(ws)
    s = 0.0
    @inbounds for pair in ws.split_bus_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = zi.y[pair.p_inj_col_i]
        pj_val = zj.y[pair.p_inj_col_j]
        qi_val = zi.y[pair.q_inj_col_i]
        qj_val = zj.y[pair.q_inj_col_j]
        s += (pi_val + pj_val)^2 + (qi_val + qj_val)^2
    end
    return s
end

function compute_residuals!(ws::SOCWRStateWorkspaceALADIN)
    ws.r_primal = sqrt(_r_pri_bus_sq_socwr(ws) + _r_pri_rt_sq_socwr(ws))
    _compute_r_progress!(ws)
    _compute_lambda_step!(ws)
    return nothing
end

function compute_residuals!(ws::SOCWRPowerWorkspaceALADIN)
    ws.r_primal = sqrt(_r_pri_pq_sq_socwr(ws))
    _compute_r_progress!(ws)
    _compute_lambda_step!(ws)
    return nothing
end

function compute_residuals!(ws::SOCWRHybridWorkspaceALADIN)
    ws.r_primal = sqrt(
        _r_pri_bus_sq_socwr(ws) +
        _r_pri_rt_sq_socwr(ws) +
        _r_pri_pq_sq_socwr(ws) +
        _r_pri_split_pq_sq_socwr(ws),
    )
    _compute_r_progress!(ws)
    _compute_lambda_step!(ws)
    return nothing
end
