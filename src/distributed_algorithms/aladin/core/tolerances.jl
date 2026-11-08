# ============================================================
# ALADIN — Convergence tolerances (Boyd-style and HFD-style)
# ============================================================

@inline function _refresh_eps_step!(ws::AbstractALADINWorkspace, eps_abs::Float64)
    ws.eps_step_primal == 0.0 &&
        (ws.eps_step_primal = _structural_nprox(ws) * eps_abs)
    return nothing
end

@inline function _norm_sq(v::Vector{Float64})::Float64
    s = 0.0
    @inbounds for x in v
        s += x * x
    end
    return s
end

# ============================================================
# Boyd-style scaling helpers
# ============================================================

function _progress_tol(
    ws::AbstractALADINWorkspace,
    eps_abs::Float64,
    eps_rel::Float64,
)::Float64
    n = 0
    ny_sq = 0.0
    nx_sq = 0.0
    for zone in ws.zones
        n += zone.n_vars
        ny_sq += _norm_sq(zone.y)
        nx_sq += _norm_sq(zone.x_prev)
    end
    return sqrt(max(n, 1)) * eps_abs + eps_rel * max(sqrt(ny_sq), sqrt(nx_sq))
end

function _dual_tol(
    ws::AbstractALADINWorkspace,
    eps_abs::Float64,
    eps_rel::Float64,
)::Float64
    n = 0
    norm_max = 0.0
    for state in values(ws.vars)
        n += length(state.lambda)
        nl = sqrt(_norm_sq(state.lambda))
        norm_max = max(norm_max, nl)
    end
    return sqrt(max(n, 1)) * eps_abs + eps_rel * norm_max
end

# ============================================================
# Boyd / HFD branch — final ε assignment
# ============================================================

@inline function _finalize_eps!(
    ws::AbstractALADINWorkspace,
    n_pairs::Int,
    norm_i_sq::Float64,
    norm_j_sq::Float64,
    eps_abs::Float64,
    eps_rel::Float64,
)
    if ws.stopping_criterion === :hfd
        ws.eps_primal = eps_abs
        _refresh_eps_step!(ws, eps_abs)
        ws.eps_dual = eps_abs
    else
        ws.eps_primal =
            sqrt(max(n_pairs, 1)) * eps_abs +
            eps_rel * max(sqrt(norm_i_sq), sqrt(norm_j_sq))
        _refresh_eps_step!(ws, eps_abs)
        ws.eps_dual = _dual_tol(ws, eps_abs, eps_rel)
    end
    return nothing
end

# ============================================================
# DC tolerances
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function compute_tolerances!(
    ws::DCStateWorkspaceALADIN,
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = length(ws.bus_pairs)
    norm_i = 0.0
    norm_j = 0.0
    @inbounds for pair in ws.bus_pairs
        yi = ws.zones[pair.zone_i].y[pair.col_i]
        yj = ws.zones[pair.zone_j].y[pair.col_j]
        norm_i += yi * yi
        norm_j += yj * yj
    end
    _finalize_eps!(ws, n, norm_i, norm_j, eps_abs, eps_rel)
    return nothing
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

@inline function _dummy_pf_norms(ws)
    ni = 0.0
    nj = 0.0
    @inbounds for pair in ws.dummy_line_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pfi = pair.b_i * (zi.y[pair.from_col_i] - zi.y[pair.to_col_i])
        pfj = pair.b_j * (zj.y[pair.from_col_j] - zj.y[pair.to_col_j])
        ni += pfi * pfi
        nj += pfj * pfj
    end
    return ni, nj
end

function compute_tolerances!(
    ws::DCPowerWorkspaceALADIN,
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = length(ws.interzonal_line_pairs) + length(ws.dummy_line_pairs)
    norm_i = 0.0
    norm_j = 0.0
    @inbounds for pair in ws.interzonal_line_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        li = pair.line_local_i
        lj = pair.line_local_j
        pfi = zi.edge_b[li] * (zi.y[zi.edge_from_cols[li]] - zi.y[zi.edge_to_cols[li]])
        pfj = zj.edge_b[lj] * (zj.y[zj.edge_from_cols[lj]] - zj.y[zj.edge_to_cols[lj]])
        norm_i += pfi * pfi
        norm_j += pfj * pfj
    end
    nd_i, nd_j = _dummy_pf_norms(ws)
    norm_i += nd_i
    norm_j += nd_j
    _finalize_eps!(ws, n, norm_i, norm_j, eps_abs, eps_rel)
    return nothing
end

# ------------------------------------------------------------
# Hybrid (Pythagoras of per-block ε for Boyd path)
# ------------------------------------------------------------

function compute_tolerances!(
    ws::DCHybridWorkspaceALADIN,
    eps_abs::Float64,
    eps_rel::Float64,
)
    if ws.stopping_criterion === :hfd
        ws.eps_primal = eps_abs
        _refresh_eps_step!(ws, eps_abs)
        ws.eps_dual = eps_abs
        return nothing
    end

    # Bus block
    n_bus = length(ws.bus_pairs)
    norm_i_bus = 0.0
    norm_j_bus = 0.0
    @inbounds for pair in ws.bus_pairs
        yi = ws.zones[pair.zone_i].y[pair.col_i]
        yj = ws.zones[pair.zone_j].y[pair.col_j]
        norm_i_bus += yi * yi
        norm_j_bus += yj * yj
    end
    eps_pri_bus =
        sqrt(max(n_bus, 1)) * eps_abs +
        eps_rel * max(sqrt(norm_i_bus), sqrt(norm_j_bus))

    # Edge block
    n_edge = length(ws.interzonal_line_pairs)
    norm_i_edge = 0.0
    norm_j_edge = 0.0
    @inbounds for pair in ws.interzonal_line_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        li = pair.line_local_i
        lj = pair.line_local_j
        pfi = zi.edge_b[li] * (zi.y[zi.edge_from_cols[li]] - zi.y[zi.edge_to_cols[li]])
        pfj = zj.edge_b[lj] * (zj.y[zj.edge_from_cols[lj]] - zj.y[zj.edge_to_cols[lj]])
        norm_i_edge += pfi * pfi
        norm_j_edge += pfj * pfj
    end
    eps_pri_edge =
        sqrt(max(n_edge, 1)) * eps_abs +
        eps_rel * max(sqrt(norm_i_edge), sqrt(norm_j_edge))

    # SplitBus block (first copy scales the "i" side, the rest the "j" side)
    n_split = length(ws.split_bus_pairs)
    norm_i_split = 0.0
    norm_j_split = 0.0
    @inbounds for pair in ws.split_bus_pairs
        for c in eachindex(pair.zones)
            yc = ws.zones[pair.zones[c]].y[pair.p_inj_cols[c]]
            if c == 1
                norm_i_split += yc * yc
            else
                norm_j_split += yc * yc
            end
        end
    end
    eps_pri_split =
        sqrt(max(n_split, 1)) * eps_abs +
        eps_rel * max(sqrt(norm_i_split), sqrt(norm_j_split))

    # DummyBus block
    n_dummy = length(ws.dummy_line_pairs)
    nd_i, nd_j = _dummy_pf_norms(ws)
    eps_pri_dummy =
        sqrt(max(n_dummy, 1)) * eps_abs +
        eps_rel * max(sqrt(nd_i), sqrt(nd_j))

    ws.eps_primal = sqrt(
        eps_pri_bus^2 + eps_pri_edge^2 + eps_pri_split^2 + eps_pri_dummy^2,
    )
    _refresh_eps_step!(ws, eps_abs)
    ws.eps_dual = _dual_tol(ws, eps_abs, eps_rel)
    return nothing
end

# ============================================================
# AC tolerances
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function compute_tolerances!(
    ws::ACStateWorkspaceALADIN,
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = 2 * length(ws.bus_pairs)
    norm_i = 0.0
    norm_j = 0.0
    @inbounds for pair in ws.bus_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        theta_i = zi.y[pair.theta_col_i]
        theta_j = zj.y[pair.theta_col_j]
        vi = zi.y[pair.vm_col_i]
        vj = zj.y[pair.vm_col_j]
        norm_i += theta_i * theta_i + vi * vi
        norm_j += theta_j * theta_j + vj * vj
    end
    _finalize_eps!(ws, n, norm_i, norm_j, eps_abs, eps_rel)
    return nothing
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

@inline function _norms_dummy_pq_ac(ws)
    ni = 0.0
    nj = 0.0
    @inbounds for pair in ws.dummy_line_pairs
        yi = ws.zones[pair.zone_i].y
        yj = ws.zones[pair.zone_j].y
        pfi = _dummy_ac_pf_at(pair, yi, :i)
        pfj = _dummy_ac_pf_at(pair, yj, :j)
        qfi = _dummy_ac_qf_at(pair, yi, :i)
        qfj = _dummy_ac_qf_at(pair, yj, :j)
        ni += pfi * pfi + qfi * qfi
        nj += pfj * pfj + qfj * qfj
    end
    return ni, nj
end

function compute_tolerances!(
    ws::ACPowerWorkspaceALADIN,
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = 2 * length(ws.interzonal_line_pairs) + 2 * length(ws.dummy_line_pairs)
    norm_i = 0.0
    norm_j = 0.0
    @inbounds for pair in ws.interzonal_line_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pfi = _ac_pf_at(zi, pair.line_local_i)
        pfj = _ac_pf_at(zj, pair.line_local_j)
        qfi = _ac_qf_at(zi, pair.line_local_i)
        qfj = _ac_qf_at(zj, pair.line_local_j)
        norm_i += pfi * pfi + qfi * qfi
        norm_j += pfj * pfj + qfj * qfj
    end
    nd_i, nd_j = _norms_dummy_pq_ac(ws)
    norm_i += nd_i
    norm_j += nd_j
    _finalize_eps!(ws, n, norm_i, norm_j, eps_abs, eps_rel)
    return nothing
end

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

function compute_tolerances!(
    ws::ACHybridWorkspaceALADIN,
    eps_abs::Float64,
    eps_rel::Float64,
)
    n =
        2 * length(ws.bus_pairs) + 2 * length(ws.interzonal_line_pairs) +
        2 * length(ws.split_bus_pairs) + 2 * length(ws.dummy_line_pairs)
    norm_i = 0.0
    norm_j = 0.0
    @inbounds for pair in ws.bus_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        theta_i = zi.y[pair.theta_col_i]
        theta_j = zj.y[pair.theta_col_j]
        vi = zi.y[pair.vm_col_i]
        vj = zj.y[pair.vm_col_j]
        norm_i += theta_i * theta_i + vi * vi
        norm_j += theta_j * theta_j + vj * vj
    end
    @inbounds for pair in ws.interzonal_line_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pfi = _ac_pf_at(zi, pair.line_local_i)
        pfj = _ac_pf_at(zj, pair.line_local_j)
        qfi = _ac_qf_at(zi, pair.line_local_i)
        qfj = _ac_qf_at(zj, pair.line_local_j)
        norm_i += pfi * pfi + qfi * qfi
        norm_j += pfj * pfj + qfj * qfj
    end
    @inbounds for pair in ws.split_bus_pairs
        for c in eachindex(pair.zones)
            zone = ws.zones[pair.zones[c]]
            pv = zone.y[pair.p_inj_cols[c]]
            qv = zone.y[pair.q_inj_cols[c]]
            if c == 1
                norm_i += pv * pv + qv * qv
            else
                norm_j += pv * pv + qv * qv
            end
        end
    end
    nd_i, nd_j = _norms_dummy_pq_ac(ws)
    norm_i += nd_i
    norm_j += nd_j
    _finalize_eps!(ws, n, norm_i, norm_j, eps_abs, eps_rel)
    return nothing
end

# ============================================================
# SOCWR tolerances
# ============================================================

@inline function _norms_bus_v2_socwr(ws)
    ni = 0.0
    nj = 0.0
    @inbounds for pair in ws.bus_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        vi = zi.y[pair.v2_col_i]
        vj = zj.y[pair.v2_col_j]
        ni += vi * vi
        nj += vj * vj
    end
    return ni, nj
end

@inline function _norms_edge_rt_socwr(ws)
    ni = 0.0
    nj = 0.0
    @inbounds for pair in ws.interzonal_line_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        ri = zi.y[zi.edge_r_cols[pair.line_local_i]]
        rj = zj.y[zj.edge_r_cols[pair.line_local_j]]
        ti = zi.y[zi.edge_t_cols[pair.line_local_i]]
        tj = zj.y[zj.edge_t_cols[pair.line_local_j]]
        ni += ri * ri + ti * ti
        nj += rj * rj + tj * tj
    end
    return ni, nj
end

@inline function _norms_split_pq_socwr(ws)
    ni = 0.0
    nj = 0.0
    @inbounds for pair in ws.split_bus_pairs
        for c in eachindex(pair.zones)
            zone = ws.zones[pair.zones[c]]
            pv = zone.y[pair.p_inj_cols[c]]
            qv = zone.y[pair.q_inj_cols[c]]
            if c == 1
                ni += pv * pv + qv * qv
            else
                nj += pv * pv + qv * qv
            end
        end
    end
    return ni, nj
end

@inline function _norms_edge_pq_socwr(ws)
    ni = 0.0
    nj = 0.0
    @inbounds for pair in ws.interzonal_line_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = _socwr_flow_value_at_y(zi, pair.line_local_i, :p)
        pj_val = _socwr_flow_value_at_y(zj, pair.line_local_j, :p)
        qi_val = _socwr_flow_value_at_y(zi, pair.line_local_i, :q)
        qj_val = _socwr_flow_value_at_y(zj, pair.line_local_j, :q)
        ni += pi_val * pi_val + qi_val * qi_val
        nj += pj_val * pj_val + qj_val * qj_val
    end
    return ni, nj
end

@inline function _norms_dummy_pq_socwr(ws)
    ni = 0.0
    nj = 0.0
    @inbounds for pair in ws.dummy_line_pairs
        yi = ws.zones[pair.zone_i].y
        yj = ws.zones[pair.zone_j].y
        pfi = _socwr_dummy_flow_value_at_y(pair, yi, :i, :p)
        pfj = _socwr_dummy_flow_value_at_y(pair, yj, :j, :p)
        qfi = _socwr_dummy_flow_value_at_y(pair, yi, :i, :q)
        qfj = _socwr_dummy_flow_value_at_y(pair, yj, :j, :q)
        ni += pfi * pfi + qfi * qfi
        nj += pfj * pfj + qfj * qfj
    end
    return ni, nj
end

function compute_tolerances!(
    ws::SOCWRStateWorkspaceALADIN,
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = length(ws.bus_pairs) + 2 * length(ws.interzonal_line_pairs)
    ni_b, nj_b = _norms_bus_v2_socwr(ws)
    ni_rt, nj_rt = _norms_edge_rt_socwr(ws)
    ni = ni_b + ni_rt
    nj = nj_b + nj_rt
    _finalize_eps!(ws, n, ni, nj, eps_abs, eps_rel)
    return nothing
end

function compute_tolerances!(
    ws::SOCWRPowerWorkspaceALADIN,
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = 2 * length(ws.interzonal_line_pairs) + 2 * length(ws.dummy_line_pairs)
    ni_e, nj_e = _norms_edge_pq_socwr(ws)
    ni_d, nj_d = _norms_dummy_pq_socwr(ws)
    ni = ni_e + ni_d
    nj = nj_e + nj_d
    _finalize_eps!(ws, n, ni, nj, eps_abs, eps_rel)
    return nothing
end

function compute_tolerances!(
    ws::SOCWRHybridWorkspaceALADIN,
    eps_abs::Float64,
    eps_rel::Float64,
)
    n =
        length(ws.bus_pairs) + 4 * length(ws.interzonal_line_pairs) +
        2 * length(ws.split_bus_pairs) + 2 * length(ws.dummy_line_pairs)
    ni_b, nj_b = _norms_bus_v2_socwr(ws)
    ni_rt, nj_rt = _norms_edge_rt_socwr(ws)
    ni_pq, nj_pq = _norms_edge_pq_socwr(ws)
    ni_sp, nj_sp = _norms_split_pq_socwr(ws)
    ni_d, nj_d = _norms_dummy_pq_socwr(ws)
    ni = ni_b + ni_rt + ni_pq + ni_sp + ni_d
    nj = nj_b + nj_rt + nj_pq + nj_sp + nj_d
    _finalize_eps!(ws, n, ni, nj, eps_abs, eps_rel)
    return nothing
end
