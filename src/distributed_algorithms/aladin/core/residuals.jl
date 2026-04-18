function _progress_sq(
    zone::Union{AladinZoneData, AladinZoneDataAC, AladinZoneDataSOCWR},
)::Float64
    s = 0.0
    y = zone.y
    x = zone.x_prev
    @inbounds for k in 1:zone.n_vars
        d = y[k] - x[k]
        s += d * d
    end
    return s
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

    prog_sq = 0.0
    @inbounds for zone in ws.zones
        prog_sq += _progress_sq(zone)
    end
    ws.r_progress = sqrt(prog_sq)

    lam_sq = 0.0
    @inbounds for k in eachindex(ws.vars.theta.lambda)
        d = ws.vars.theta.lambda[k] - ws.vars.theta.lambda_prev[k]
        lam_sq += d * d
    end
    ws.lambda_step = sqrt(lam_sq)
    return nothing
end

# ---- Power ----

@inline function _dummy_pf_diff_sq(ws)
    s = 0.0
    @inbounds for pair in ws.dummy_branch_pairs
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
    @inbounds for pair in ws.edge_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        li = pair.edge_local_i
        lj = pair.edge_local_j
        pfi = zi.edge_b[li] * (zi.y[zi.edge_from_cols[li]] - zi.y[zi.edge_to_cols[li]])
        pfj = zj.edge_b[lj] * (zj.y[zj.edge_from_cols[lj]] - zj.y[zj.edge_to_cols[lj]])
        r_pri_sq += (pfi - pfj)^2
    end
    r_pri_sq += _dummy_pf_diff_sq(ws)
    ws.r_primal = sqrt(r_pri_sq)

    prog_sq = 0.0
    @inbounds for zone in ws.zones
        prog_sq += _progress_sq(zone)
    end
    ws.r_progress = sqrt(prog_sq)

    lam_sq = 0.0
    @inbounds for k in eachindex(ws.vars.p_flow.lambda)
        d = ws.vars.p_flow.lambda[k] - ws.vars.p_flow.lambda_prev[k]
        lam_sq += d * d
    end
    @inbounds for k in eachindex(ws.vars.p_dummy.lambda)
        d = ws.vars.p_dummy.lambda[k] - ws.vars.p_dummy.lambda_prev[k]
        lam_sq += d * d
    end
    ws.lambda_step = sqrt(lam_sq)
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
    @inbounds for pair in ws.edge_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        li = pair.edge_local_i
        lj = pair.edge_local_j
        pfi = zi.edge_b[li] * (zi.y[zi.edge_from_cols[li]] - zi.y[zi.edge_to_cols[li]])
        pfj = zj.edge_b[lj] * (zj.y[zj.edge_from_cols[lj]] - zj.y[zj.edge_to_cols[lj]])
        r_pri_edge_sq += (pfi - pfj)^2
    end
    ws.r_primal_edge = sqrt(r_pri_edge_sq)

    # SplitBus p_inj consensus block: residual is p_inj_i + p_inj_j (target 0).
    r_pri_split_sq = 0.0
    @inbounds for pair in ws.split_bus_pairs
        yi = ws.zones[pair.zone_i].y[pair.pinj_col_i]
        yj = ws.zones[pair.zone_j].y[pair.pinj_col_j]
        r_pri_split_sq += (yi + yj)^2
    end

    # DummyBus half-branch flow consensus (p_i - p_j = 0).
    r_pri_dummy_sq = _dummy_pf_diff_sq(ws)

    ws.r_primal = sqrt(
        r_pri_bus_sq + r_pri_edge_sq + r_pri_split_sq + r_pri_dummy_sq,
    )

    # Progress block (shared across couplings: it is a single per-zone quantity)
    prog_sq = 0.0
    @inbounds for zone in ws.zones
        prog_sq += _progress_sq(zone)
    end
    ws.r_progress = sqrt(prog_sq)
    ws.r_progress_bus = ws.r_progress       # same scalar, kept for diagnostics
    ws.r_progress_edge = ws.r_progress

    lam_sq = 0.0
    @inbounds for k in eachindex(ws.vars.theta.lambda)
        d = ws.vars.theta.lambda[k] - ws.vars.theta.lambda_prev[k]
        lam_sq += d * d
    end
    @inbounds for k in eachindex(ws.vars.p_flow.lambda)
        d = ws.vars.p_flow.lambda[k] - ws.vars.p_flow.lambda_prev[k]
        lam_sq += d * d
    end
    @inbounds for k in eachindex(ws.vars.p_inj.lambda)
        d = ws.vars.p_inj.lambda[k] - ws.vars.p_inj.lambda_prev[k]
        lam_sq += d * d
    end
    @inbounds for k in eachindex(ws.vars.p_dummy.lambda)
        d = ws.vars.p_dummy.lambda[k] - ws.vars.p_dummy.lambda_prev[k]
        lam_sq += d * d
    end
    ws.lambda_step = sqrt(lam_sq)
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
        θi = zi.y[pair.theta_col_i]
        θj = zj.y[pair.theta_col_j]
        vi = zi.y[pair.vm_col_i]
        vj = zj.y[pair.vm_col_j]
        r_pri_sq += (θi - θj)^2 + (vi - vj)^2
    end
    ws.r_primal = sqrt(r_pri_sq)

    prog_sq = 0.0
    @inbounds for zone in ws.zones
        prog_sq += _progress_sq(zone)
    end
    ws.r_progress = sqrt(prog_sq)

    lam_sq = 0.0
    @inbounds for k in eachindex(ws.vars.theta.lambda)
        dt = ws.vars.theta.lambda[k] - ws.vars.theta.lambda_prev[k]
        dv = ws.vars.v_mag.lambda[k] - ws.vars.v_mag.lambda_prev[k]
        lam_sq += dt * dt + dv * dv
    end
    ws.lambda_step = sqrt(lam_sq)
    return nothing
end

@inline function _dummy_ac_pf_at(pair::DummyBranchPairAC, y::Vector{Float64}, side::Symbol)
    if side === :i
        θf = y[pair.from_theta_col_i]
        θt = y[pair.to_theta_col_i]
        vf = y[pair.from_vm_col_i]
        vt = y[pair.to_vm_col_i]
        g = pair.g_i
        b = pair.b_i
        tau = pair.inv_tap_i
        tap_angle = pair.tap_angle_i
    else
        θf = y[pair.from_theta_col_j]
        θt = y[pair.to_theta_col_j]
        vf = y[pair.from_vm_col_j]
        vt = y[pair.to_vm_col_j]
        g = pair.g_j
        b = pair.b_j
        tau = pair.inv_tap_j
        tap_angle = pair.tap_angle_j
    end
    Δθ = θf - θt - tap_angle
    return g * (tau * vf)^2 - tau * vf * vt * (g * cos(Δθ) + b * sin(Δθ))
end

@inline function _dummy_ac_qf_at(pair::DummyBranchPairAC, y::Vector{Float64}, side::Symbol)
    if side === :i
        θf = y[pair.from_theta_col_i]
        θt = y[pair.to_theta_col_i]
        vf = y[pair.from_vm_col_i]
        vt = y[pair.to_vm_col_i]
        g = pair.g_i
        b = pair.b_i
        b_sh = pair.b_sh_i
        tau = pair.inv_tap_i
        tap_angle = pair.tap_angle_i
    else
        θf = y[pair.from_theta_col_j]
        θt = y[pair.to_theta_col_j]
        vf = y[pair.from_vm_col_j]
        vt = y[pair.to_vm_col_j]
        g = pair.g_j
        b = pair.b_j
        b_sh = pair.b_sh_j
        tau = pair.inv_tap_j
        tap_angle = pair.tap_angle_j
    end
    Δθ = θf - θt - tap_angle
    return -(b + b_sh / 2) * (tau * vf)^2 -
           tau * vf * vt * (-b * cos(Δθ) + g * sin(Δθ))
end

@inline function _r_pri_dummy_pq_sq_ac(ws)
    s = 0.0
    @inbounds for pair in ws.dummy_branch_pairs
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
    @inbounds for pair in ws.edge_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pfi = _ac_pf_at(zi, pair.edge_local_i)
        pfj = _ac_pf_at(zj, pair.edge_local_j)
        qfi = _ac_qf_at(zi, pair.edge_local_i)
        qfj = _ac_qf_at(zj, pair.edge_local_j)
        r_pri_sq += (pfi - pfj)^2 + (qfi - qfj)^2
    end
    r_pri_sq += _r_pri_dummy_pq_sq_ac(ws)
    ws.r_primal = sqrt(r_pri_sq)

    prog_sq = 0.0
    @inbounds for zone in ws.zones
        prog_sq += _progress_sq(zone)
    end
    ws.r_progress = sqrt(prog_sq)

    lam_sq = 0.0
    @inbounds for k in eachindex(ws.vars.p_flow.lambda)
        dp = ws.vars.p_flow.lambda[k] - ws.vars.p_flow.lambda_prev[k]
        dq = ws.vars.q_flow.lambda[k] - ws.vars.q_flow.lambda_prev[k]
        lam_sq += dp * dp + dq * dq
    end
    @inbounds for k in eachindex(ws.vars.p_dummy.lambda)
        dp = ws.vars.p_dummy.lambda[k] - ws.vars.p_dummy.lambda_prev[k]
        dq = ws.vars.q_dummy.lambda[k] - ws.vars.q_dummy.lambda_prev[k]
        lam_sq += dp * dp + dq * dq
    end
    ws.lambda_step = sqrt(lam_sq)
    return nothing
end

function compute_residuals!(ws::ACHybridWorkspaceALADIN)
    r_pri_sq = 0.0
    @inbounds for pair in ws.bus_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        θi = zi.y[pair.theta_col_i]
        θj = zj.y[pair.theta_col_j]
        vi = zi.y[pair.vm_col_i]
        vj = zj.y[pair.vm_col_j]
        r_pri_sq += (θi - θj)^2 + (vi - vj)^2
    end
    @inbounds for pair in ws.edge_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pfi = _ac_pf_at(zi, pair.edge_local_i)
        pfj = _ac_pf_at(zj, pair.edge_local_j)
        qfi = _ac_qf_at(zi, pair.edge_local_i)
        qfj = _ac_qf_at(zj, pair.edge_local_j)
        r_pri_sq += (pfi - pfj)^2 + (qfi - qfj)^2
    end
    @inbounds for pair in ws.split_bus_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = zi.y[pair.pinj_col_i]
        pj_val = zj.y[pair.pinj_col_j]
        qi_val = zi.y[pair.qinj_col_i]
        qj_val = zj.y[pair.qinj_col_j]
        r_pri_sq += (pi_val + pj_val)^2 + (qi_val + qj_val)^2
    end
    r_pri_sq += _r_pri_dummy_pq_sq_ac(ws)
    ws.r_primal = sqrt(r_pri_sq)

    prog_sq = 0.0
    @inbounds for zone in ws.zones
        prog_sq += _progress_sq(zone)
    end
    ws.r_progress = sqrt(prog_sq)

    lam_sq = 0.0
    @inbounds for k in eachindex(ws.vars.theta.lambda)
        dt = ws.vars.theta.lambda[k] - ws.vars.theta.lambda_prev[k]
        dv = ws.vars.v_mag.lambda[k] - ws.vars.v_mag.lambda_prev[k]
        lam_sq += dt * dt + dv * dv
    end
    @inbounds for k in eachindex(ws.vars.p_flow.lambda)
        dp = ws.vars.p_flow.lambda[k] - ws.vars.p_flow.lambda_prev[k]
        dq = ws.vars.q_flow.lambda[k] - ws.vars.q_flow.lambda_prev[k]
        lam_sq += dp * dp + dq * dq
    end
    @inbounds for k in eachindex(ws.vars.p_inj.lambda)
        dp = ws.vars.p_inj.lambda[k] - ws.vars.p_inj.lambda_prev[k]
        dq = ws.vars.q_inj.lambda[k] - ws.vars.q_inj.lambda_prev[k]
        lam_sq += dp * dp + dq * dq
    end
    @inbounds for k in eachindex(ws.vars.p_dummy.lambda)
        dp = ws.vars.p_dummy.lambda[k] - ws.vars.p_dummy.lambda_prev[k]
        dq = ws.vars.q_dummy.lambda[k] - ws.vars.q_dummy.lambda_prev[k]
        lam_sq += dp * dp + dq * dq
    end
    ws.lambda_step = sqrt(lam_sq)
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

@inline function _r_pri_split_pq_sq_socwr(ws)
    s = 0.0
    @inbounds for pair in ws.split_bus_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = zi.y[pair.pinj_col_i]
        pj_val = zj.y[pair.pinj_col_j]
        qi_val = zi.y[pair.qinj_col_i]
        qj_val = zj.y[pair.qinj_col_j]
        s += (pi_val + pj_val)^2 + (qi_val + qj_val)^2
    end
    return s
end

@inline function _sum_progress_sq_socwr(ws)
    s = 0.0
    @inbounds for zone in ws.zones
        s += _progress_sq(zone)
    end
    return s
end

function compute_residuals!(ws::SOCWRStateWorkspaceALADIN)
    ws.r_primal = sqrt(_r_pri_bus_sq_socwr(ws) + _r_pri_rt_sq_socwr(ws))
    ws.r_progress = sqrt(_sum_progress_sq_socwr(ws))

    lam_sq = 0.0
    @inbounds for k in eachindex(ws.vars.v2.lambda)
        d = ws.vars.v2.lambda[k] - ws.vars.v2.lambda_prev[k]
        lam_sq += d * d
    end
    @inbounds for k in eachindex(ws.vars.r.lambda)
        dr = ws.vars.r.lambda[k] - ws.vars.r.lambda_prev[k]
        dt = ws.vars.t.lambda[k] - ws.vars.t.lambda_prev[k]
        lam_sq += dr * dr + dt * dt
    end
    ws.lambda_step = sqrt(lam_sq)
    return nothing
end

function compute_residuals!(ws::SOCWRPowerWorkspaceALADIN)
    ws.r_primal = sqrt(_r_pri_pq_sq_socwr(ws))
    ws.r_progress = sqrt(_sum_progress_sq_socwr(ws))

    lam_sq = 0.0
    @inbounds for k in eachindex(ws.vars.p_flow.lambda)
        dp = ws.vars.p_flow.lambda[k] - ws.vars.p_flow.lambda_prev[k]
        dq = ws.vars.q_flow.lambda[k] - ws.vars.q_flow.lambda_prev[k]
        lam_sq += dp * dp + dq * dq
    end
    ws.lambda_step = sqrt(lam_sq)
    return nothing
end

function compute_residuals!(ws::SOCWRHybridWorkspaceALADIN)
    ws.r_primal = sqrt(
        _r_pri_bus_sq_socwr(ws) +
        _r_pri_rt_sq_socwr(ws) +
        _r_pri_pq_sq_socwr(ws) +
        _r_pri_split_pq_sq_socwr(ws),
    )
    ws.r_progress = sqrt(_sum_progress_sq_socwr(ws))

    lam_sq = 0.0
    @inbounds for k in eachindex(ws.vars.v2.lambda)
        d = ws.vars.v2.lambda[k] - ws.vars.v2.lambda_prev[k]
        lam_sq += d * d
    end
    @inbounds for k in eachindex(ws.vars.r.lambda)
        dr = ws.vars.r.lambda[k] - ws.vars.r.lambda_prev[k]
        dt = ws.vars.t.lambda[k] - ws.vars.t.lambda_prev[k]
        dp = ws.vars.p_flow.lambda[k] - ws.vars.p_flow.lambda_prev[k]
        dq = ws.vars.q_flow.lambda[k] - ws.vars.q_flow.lambda_prev[k]
        lam_sq += dr * dr + dt * dt + dp * dp + dq * dq
    end
    @inbounds for k in eachindex(ws.vars.p_inj.lambda)
        dp = ws.vars.p_inj.lambda[k] - ws.vars.p_inj.lambda_prev[k]
        dq = ws.vars.q_inj.lambda[k] - ws.vars.q_inj.lambda_prev[k]
        lam_sq += dp * dp + dq * dq
    end
    ws.lambda_step = sqrt(lam_sq)
    return nothing
end
