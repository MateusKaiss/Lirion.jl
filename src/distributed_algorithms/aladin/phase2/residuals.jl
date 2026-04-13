# ============================================================
# ALADIN — Residuals (canonical HFD 2016 Theorem 4.1)
# ------------------------------------------------------------
# Primal residual (consensus feasibility):
#     r_primal = ‖∑_i A_i y_i‖₂
#   For our pairwise formulation this equals
#     r_primal² = ∑_pairs (coupling_i - coupling_j)²
#
# Progress residual (Phase-1 KKT proxy):
#     r_progress = √(∑_i ‖y_i - x_iᵏ‖²)
#   This is the quantity that appears in HFD 2016 Theorem 4.1 as the
#   sufficient condition for KKT convergence: when y_i stabilises at x_iᵏ
#   across all zones, every local subproblem is at its KKT point with the
#   current Lagrange multipliers, and the consensus residual captures the
#   remaining feasibility gap.
#
# Diagnostic (not part of the stopping rule):
#     lambda_step = ‖λᵏ⁺¹ - λᵏ‖₂
#   Tracked to monitor how the outer dual is moving.
# ============================================================

function _progress_sq(zone::AladinZoneData)::Float64
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
    @inbounds for k in eachindex(ws.lambda_bus)
        d = ws.lambda_bus[k] - ws.lambda_bus_prev[k]
        lam_sq += d * d
    end
    ws.lambda_step = sqrt(lam_sq)
    return nothing
end

# ---- Power ----

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
    ws.r_primal = sqrt(r_pri_sq)

    prog_sq = 0.0
    @inbounds for zone in ws.zones
        prog_sq += _progress_sq(zone)
    end
    ws.r_progress = sqrt(prog_sq)

    lam_sq = 0.0
    @inbounds for k in eachindex(ws.lambda_edge)
        d = ws.lambda_edge[k] - ws.lambda_edge_prev[k]
        lam_sq += d * d
    end
    ws.lambda_step = sqrt(lam_sq)
    return nothing
end

# ---- Hybrid ----

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

    ws.r_primal = hypot(ws.r_primal_bus, ws.r_primal_edge)

    # Progress block (shared across couplings: it is a single per-zone quantity)
    prog_sq = 0.0
    @inbounds for zone in ws.zones
        prog_sq += _progress_sq(zone)
    end
    ws.r_progress = sqrt(prog_sq)
    ws.r_progress_bus = ws.r_progress       # same scalar, kept for diagnostics
    ws.r_progress_edge = ws.r_progress

    lam_sq = 0.0
    @inbounds for k in eachindex(ws.lambda_bus)
        d = ws.lambda_bus[k] - ws.lambda_bus_prev[k]
        lam_sq += d * d
    end
    @inbounds for k in eachindex(ws.lambda_edge)
        d = ws.lambda_edge[k] - ws.lambda_edge_prev[k]
        lam_sq += d * d
    end
    ws.lambda_step = sqrt(lam_sq)
    return nothing
end
