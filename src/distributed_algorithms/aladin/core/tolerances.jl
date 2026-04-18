# ============================================================
# ALADIN — Convergence tolerances (Boyd-style, HFD-aligned)
# ============================================================

@inline function _norm_sq(v::Vector{Float64})::Float64
    s = 0.0
    @inbounds for x in v
        s += x * x
    end
    return s
end

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

# ---- State ----

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
    ws.eps_primal =
        sqrt(max(n, 1)) * eps_abs +
        eps_rel * max(sqrt(norm_i), sqrt(norm_j))
    ws.eps_progress = _progress_tol(ws, eps_abs, eps_rel)
    return nothing
end

# ---- Power ----

@inline function _dummy_pf_norms(ws)
    ni = 0.0
    nj = 0.0
    @inbounds for pair in ws.dummy_branch_pairs
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
    n = length(ws.edge_pairs) + length(ws.dummy_branch_pairs)
    norm_i = 0.0
    norm_j = 0.0
    @inbounds for pair in ws.edge_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        li = pair.edge_local_i
        lj = pair.edge_local_j
        pfi = zi.edge_b[li] * (zi.y[zi.edge_from_cols[li]] - zi.y[zi.edge_to_cols[li]])
        pfj = zj.edge_b[lj] * (zj.y[zj.edge_from_cols[lj]] - zj.y[zj.edge_to_cols[lj]])
        norm_i += pfi * pfi
        norm_j += pfj * pfj
    end
    nd_i, nd_j = _dummy_pf_norms(ws)
    norm_i += nd_i
    norm_j += nd_j
    ws.eps_primal =
        sqrt(max(n, 1)) * eps_abs +
        eps_rel * max(sqrt(norm_i), sqrt(norm_j))
    ws.eps_progress = _progress_tol(ws, eps_abs, eps_rel)
    return nothing
end

# ---- Hybrid ----

function compute_tolerances!(
    ws::DCHybridWorkspaceALADIN,
    eps_abs::Float64,
    eps_rel::Float64,
)
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
    n_edge = length(ws.edge_pairs)
    norm_i_edge = 0.0
    norm_j_edge = 0.0
    @inbounds for pair in ws.edge_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        li = pair.edge_local_i
        lj = pair.edge_local_j
        pfi = zi.edge_b[li] * (zi.y[zi.edge_from_cols[li]] - zi.y[zi.edge_to_cols[li]])
        pfj = zj.edge_b[lj] * (zj.y[zj.edge_from_cols[lj]] - zj.y[zj.edge_to_cols[lj]])
        norm_i_edge += pfi * pfi
        norm_j_edge += pfj * pfj
    end
    eps_pri_edge =
        sqrt(max(n_edge, 1)) * eps_abs +
        eps_rel * max(sqrt(norm_i_edge), sqrt(norm_j_edge))

    # SplitBus block: norm on p_inj values (sign-flipped consensus)
    n_split = length(ws.split_bus_pairs)
    norm_i_split = 0.0
    norm_j_split = 0.0
    @inbounds for pair in ws.split_bus_pairs
        yi = ws.zones[pair.zone_i].y[pair.pinj_col_i]
        yj = ws.zones[pair.zone_j].y[pair.pinj_col_j]
        norm_i_split += yi * yi
        norm_j_split += yj * yj
    end
    eps_pri_split =
        sqrt(max(n_split, 1)) * eps_abs +
        eps_rel * max(sqrt(norm_i_split), sqrt(norm_j_split))

    # DummyBus block: p_flow norms on half-branches
    n_dummy = length(ws.dummy_branch_pairs)
    nd_i, nd_j = _dummy_pf_norms(ws)
    eps_pri_dummy =
        sqrt(max(n_dummy, 1)) * eps_abs +
        eps_rel * max(sqrt(nd_i), sqrt(nd_j))

    ws.eps_primal = sqrt(
        eps_pri_bus^2 + eps_pri_edge^2 + eps_pri_split^2 + eps_pri_dummy^2,
    )
    ws.eps_progress = _progress_tol(ws, eps_abs, eps_rel)
    return nothing
end

# ------------------------------------------------------------
# AC tolerances
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
        θi = zi.y[pair.theta_col_i]
        θj = zj.y[pair.theta_col_j]
        vi = zi.y[pair.vm_col_i]
        vj = zj.y[pair.vm_col_j]
        norm_i += θi * θi + vi * vi
        norm_j += θj * θj + vj * vj
    end
    ws.eps_primal =
        sqrt(max(n, 1)) * eps_abs +
        eps_rel * max(sqrt(norm_i), sqrt(norm_j))
    ws.eps_progress = _progress_tol(ws, eps_abs, eps_rel)
    return nothing
end

function compute_tolerances!(
    ws::ACPowerWorkspaceALADIN,
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = 2 * length(ws.edge_pairs)
    norm_i = 0.0
    norm_j = 0.0
    @inbounds for pair in ws.edge_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pfi = _ac_pf_at(zi, pair.edge_local_i)
        pfj = _ac_pf_at(zj, pair.edge_local_j)
        qfi = _ac_qf_at(zi, pair.edge_local_i)
        qfj = _ac_qf_at(zj, pair.edge_local_j)
        norm_i += pfi * pfi + qfi * qfi
        norm_j += pfj * pfj + qfj * qfj
    end
    ws.eps_primal =
        sqrt(max(n, 1)) * eps_abs +
        eps_rel * max(sqrt(norm_i), sqrt(norm_j))
    ws.eps_progress = _progress_tol(ws, eps_abs, eps_rel)
    return nothing
end

function compute_tolerances!(
    ws::ACHybridWorkspaceALADIN,
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = 2 * length(ws.bus_pairs) + 2 * length(ws.edge_pairs) +
        2 * length(ws.split_bus_pairs)
    norm_i = 0.0
    norm_j = 0.0
    @inbounds for pair in ws.bus_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        θi = zi.y[pair.theta_col_i]
        θj = zj.y[pair.theta_col_j]
        vi = zi.y[pair.vm_col_i]
        vj = zj.y[pair.vm_col_j]
        norm_i += θi * θi + vi * vi
        norm_j += θj * θj + vj * vj
    end
    @inbounds for pair in ws.edge_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pfi = _ac_pf_at(zi, pair.edge_local_i)
        pfj = _ac_pf_at(zj, pair.edge_local_j)
        qfi = _ac_qf_at(zi, pair.edge_local_i)
        qfj = _ac_qf_at(zj, pair.edge_local_j)
        norm_i += pfi * pfi + qfi * qfi
        norm_j += pfj * pfj + qfj * qfj
    end
    @inbounds for pair in ws.split_bus_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = zi.y[pair.pinj_col_i]
        pj_val = zj.y[pair.pinj_col_j]
        qi_val = zi.y[pair.qinj_col_i]
        qj_val = zj.y[pair.qinj_col_j]
        norm_i += pi_val * pi_val + qi_val * qi_val
        norm_j += pj_val * pj_val + qj_val * qj_val
    end
    ws.eps_primal =
        sqrt(max(n, 1)) * eps_abs +
        eps_rel * max(sqrt(norm_i), sqrt(norm_j))
    ws.eps_progress = _progress_tol(ws, eps_abs, eps_rel)
    return nothing
end

# ------------------------------------------------------------
# SOCWR tolerances
# ------------------------------------------------------------

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
    @inbounds for pair in ws.edge_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        ri = zi.y[zi.edge_r_cols[pair.edge_local_i]]
        rj = zj.y[zj.edge_r_cols[pair.edge_local_j]]
        ti = zi.y[zi.edge_t_cols[pair.edge_local_i]]
        tj = zj.y[zj.edge_t_cols[pair.edge_local_j]]
        ni += ri * ri + ti * ti
        nj += rj * rj + tj * tj
    end
    return ni, nj
end

@inline function _norms_split_pq_socwr(ws)
    ni = 0.0
    nj = 0.0
    @inbounds for pair in ws.split_bus_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = zi.y[pair.pinj_col_i]
        pj_val = zj.y[pair.pinj_col_j]
        qi_val = zi.y[pair.qinj_col_i]
        qj_val = zj.y[pair.qinj_col_j]
        ni += pi_val * pi_val + qi_val * qi_val
        nj += pj_val * pj_val + qj_val * qj_val
    end
    return ni, nj
end

@inline function _norms_edge_pq_socwr(ws)
    ni = 0.0
    nj = 0.0
    @inbounds for pair in ws.edge_pairs
        zi = ws.zones[pair.zone_i]
        zj = ws.zones[pair.zone_j]
        pi_val = _socwr_flow_value_at_y(zi, pair.edge_local_i, :p)
        pj_val = _socwr_flow_value_at_y(zj, pair.edge_local_j, :p)
        qi_val = _socwr_flow_value_at_y(zi, pair.edge_local_i, :q)
        qj_val = _socwr_flow_value_at_y(zj, pair.edge_local_j, :q)
        ni += pi_val * pi_val + qi_val * qi_val
        nj += pj_val * pj_val + qj_val * qj_val
    end
    return ni, nj
end

function compute_tolerances!(
    ws::SOCWRStateWorkspaceALADIN,
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = length(ws.bus_pairs) + 2 * length(ws.edge_pairs)
    ni_b, nj_b = _norms_bus_v2_socwr(ws)
    ni_rt, nj_rt = _norms_edge_rt_socwr(ws)
    ni = ni_b + ni_rt
    nj = nj_b + nj_rt
    ws.eps_primal =
        sqrt(max(n, 1)) * eps_abs +
        eps_rel * max(sqrt(ni), sqrt(nj))
    ws.eps_progress = _progress_tol(ws, eps_abs, eps_rel)
    return nothing
end

function compute_tolerances!(
    ws::SOCWRPowerWorkspaceALADIN,
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = 2 * length(ws.edge_pairs)
    ni, nj = _norms_edge_pq_socwr(ws)
    ws.eps_primal =
        sqrt(max(n, 1)) * eps_abs +
        eps_rel * max(sqrt(ni), sqrt(nj))
    ws.eps_progress = _progress_tol(ws, eps_abs, eps_rel)
    return nothing
end

function compute_tolerances!(
    ws::SOCWRHybridWorkspaceALADIN,
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = length(ws.bus_pairs) + 4 * length(ws.edge_pairs) +
        2 * length(ws.split_bus_pairs)
    ni_b, nj_b = _norms_bus_v2_socwr(ws)
    ni_rt, nj_rt = _norms_edge_rt_socwr(ws)
    ni_pq, nj_pq = _norms_edge_pq_socwr(ws)
    ni_sp, nj_sp = _norms_split_pq_socwr(ws)
    ni = ni_b + ni_rt + ni_pq + ni_sp
    nj = nj_b + nj_rt + nj_pq + nj_sp
    ws.eps_primal =
        sqrt(max(n, 1)) * eps_abs +
        eps_rel * max(sqrt(ni), sqrt(nj))
    ws.eps_progress = _progress_tol(ws, eps_abs, eps_rel)
    return nothing
end
