# ============================================================
# ALADIN — Convergence tolerances (Boyd-style, HFD-aligned)
# ------------------------------------------------------------
# Primal (consensus feasibility):
#     eps_primal = √n_pri · eps_abs
#                  + eps_rel · max(‖coupling_i‖₂, ‖coupling_j‖₂)
#
# Progress (Phase-1 KKT proxy, HFD 2016 Thm 4.1):
#     eps_progress = √n_prog · eps_abs
#                    + eps_rel · max(‖y‖₂, ‖x_prev‖₂)
# where n_prog = ∑_i n_vars_i (total primal dimension across zones).
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

function compute_tolerances!(
    ws::DCPowerWorkspaceALADIN,
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = length(ws.edge_pairs)
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

    ws.eps_primal = hypot(eps_pri_bus, eps_pri_edge)
    ws.eps_progress = _progress_tol(ws, eps_abs, eps_rel)
    return nothing
end
