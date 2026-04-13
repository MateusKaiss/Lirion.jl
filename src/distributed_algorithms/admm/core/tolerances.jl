# ============================================================
# ADMM Tolerances (Boyd et al.)
# ------------------------------------------------------------
#   eps_pri  = sqrt(n) * eps_abs + eps_rel * max(‖Ax‖, ‖Bz‖)
#   eps_dual = sqrt(n) * eps_abs + eps_rel * ‖y‖
#
# All quantities are summed over expanded (zone, var) pairs to match
# the expanded form of the residuals (Boyd canonical consensus form).
# ============================================================

# ------------------------------------------------------------
# Bus tolerance helpers
# ------------------------------------------------------------

function _eps_primal_bus(
    partitions::Partitions,
    consensus_vecs,
    snapshot_fns,
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = 0
    norm_x_sq = 0.0
    norm_z_sq = 0.0
    for part in partitions
        bufs = map(fn -> fn(part), snapshot_fns)
        n += length(snapshot_fns) * length(part.zone_shared_bus_ids)
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            li = part.shared_bus_local_idx[i]
            cj = part.shared_bus_consensus_idx[i]
            for k in eachindex(bufs)
                norm_x_sq += bufs[k][li]^2
                norm_z_sq += consensus_vecs[k][cj]^2
            end
        end
    end
    n = max(n, 1)

    return sqrt(n) * eps_abs + eps_rel * max(sqrt(norm_x_sq), sqrt(norm_z_sq))
end

function _eps_dual_bus(
    duals_list,       # tuple of Matrix{Float64} (n_shared × n_zones)
    partitions::Partitions,
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = 0
    norm_y_sq = 0.0
    for part in partitions
        p = part.position
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            j = part.shared_bus_consensus_idx[i]
            for y in duals_list
                n += 1
                norm_y_sq += y[j, p]^2
            end
        end
    end
    n = max(n, 1)
    return sqrt(n) * eps_abs + eps_rel * sqrt(norm_y_sq)
end

# ------------------------------------------------------------
# Branch tolerance helpers
# ------------------------------------------------------------

function _eps_primal_branch(
    partitions::Partitions,
    consensus_vecs,            # tuple of Vector{Float64}, length n_inter
    interzonal::Vector{Int},
    n_vars_per_edge::Int,
    eps_abs::Float64,
    eps_rel::Float64,
)
    dim = 0
    norm_z_sq = 0.0
    for part in partitions
        @inbounds for (k, e) in enumerate(interzonal)
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            dim += n_vars_per_edge
            for z in consensus_vecs
                norm_z_sq += z[k]^2
            end
        end
    end

    return sqrt(max(1, dim)) * eps_abs + eps_rel * max(sqrt(norm_z_sq), 1e-12)
end

function _eps_dual_branch(
    duals_list,         # tuple of Matrix{Float64} (n_inter × n_zones)
    partitions::Partitions,
    interzonal::Vector{Int},
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = 0
    norm_y_sq = 0.0
    for part in partitions
        p = part.position
        @inbounds for (k, e) in enumerate(interzonal)
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            for y in duals_list
                n += 1
                norm_y_sq += y[k, p]^2
            end
        end
    end
    n = max(n, 1)
    return sqrt(n) * eps_abs + eps_rel * sqrt(norm_y_sq)
end

# ############################################################
#  DISPATCH: compute_tolerances!(ws, eps_abs, eps_rel)
# ############################################################

# ---- DC ----

function compute_tolerances!(ws::DCStateWorkspace, eps_abs::Float64, eps_rel::Float64)
    eps_r = _eps_primal_bus(
        ws.partitions,
        (ws.theta_consensus,),
        (snapshot_theta!,),
        eps_abs,
        eps_rel,
    )
    eps_s = _eps_dual_bus(
        (ws.duals_theta,),
        ws.partitions,
        eps_abs,
        eps_rel,
    )
    ws.eps_primal = eps_r
    ws.eps_dual = eps_s
    return nothing
end

function compute_tolerances!(ws::DCPowerWorkspace, eps_abs::Float64, eps_rel::Float64)
    eps_r = _eps_primal_branch(
        ws.partitions,
        (ws.p_flow_consensus,),
        ws.interzonal,
        1,
        eps_abs,
        eps_rel,
    )
    eps_s = _eps_dual_branch(
        (ws.duals_p_flow,),
        ws.partitions,
        ws.interzonal,
        eps_abs,
        eps_rel,
    )
    ws.eps_primal = eps_r
    ws.eps_dual = eps_s
    return nothing
end

function compute_tolerances!(ws::DCHybridWorkspace, eps_abs::Float64, eps_rel::Float64)
    eps_r = hypot(
        _eps_primal_bus(
            ws.partitions,
            (ws.theta_consensus,),
            (snapshot_theta!,),
            eps_abs,
            eps_rel,
        ),
        _eps_primal_branch(
            ws.partitions,
            (ws.p_flow_consensus,),
            ws.interzonal,
            1,
            eps_abs,
            eps_rel,
        ),
    )
    eps_s = hypot(
        _eps_dual_bus(
            (ws.duals_theta,),
            ws.partitions,
            eps_abs,
            eps_rel,
        ),
        _eps_dual_branch(
            (ws.duals_p_flow,),
            ws.partitions,
            ws.interzonal,
            eps_abs,
            eps_rel,
        ),
    )
    ws.eps_primal = eps_r
    ws.eps_dual = eps_s
    return nothing
end

# ---- AC ----

function compute_tolerances!(ws::ACStateWorkspace, eps_abs::Float64, eps_rel::Float64)
    eps_r = _eps_primal_bus(
        ws.partitions,
        (ws.theta_consensus, ws.v_mag_consensus),
        (snapshot_theta!, snapshot_v_mag!),
        eps_abs,
        eps_rel,
    )
    eps_s = _eps_dual_bus(
        (ws.duals_theta, ws.duals_v_mag),
        ws.partitions,
        eps_abs,
        eps_rel,
    )
    ws.eps_primal = eps_r
    ws.eps_dual = eps_s
    return nothing
end

function compute_tolerances!(ws::ACPowerWorkspace, eps_abs::Float64, eps_rel::Float64)
    eps_r = _eps_primal_branch(
        ws.partitions,
        (ws.p_flow_consensus, ws.q_flow_consensus),
        ws.interzonal,
        2,
        eps_abs,
        eps_rel,
    )
    eps_s = _eps_dual_branch(
        (ws.duals_p_flow, ws.duals_q_flow),
        ws.partitions,
        ws.interzonal,
        eps_abs,
        eps_rel,
    )
    ws.eps_primal = eps_r
    ws.eps_dual = eps_s
    return nothing
end

function compute_tolerances!(ws::ACHybridWorkspace, eps_abs::Float64, eps_rel::Float64)
    eps_r = hypot(
        _eps_primal_bus(
            ws.partitions,
            (ws.theta_consensus, ws.v_mag_consensus),
            (snapshot_theta!, snapshot_v_mag!),
            eps_abs,
            eps_rel,
        ),
        _eps_primal_branch(
            ws.partitions,
            (ws.p_flow_consensus, ws.q_flow_consensus),
            ws.interzonal,
            2,
            eps_abs,
            eps_rel,
        ),
    )
    eps_s = hypot(
        _eps_dual_bus(
            (ws.duals_theta, ws.duals_v_mag),
            ws.partitions,
            eps_abs,
            eps_rel,
        ),
        _eps_dual_branch(
            (ws.duals_p_flow, ws.duals_q_flow),
            ws.partitions,
            ws.interzonal,
            eps_abs,
            eps_rel,
        ),
    )
    ws.eps_primal = eps_r
    ws.eps_dual = eps_s
    return nothing
end

# ---- SOCWR ----

function compute_tolerances!(ws::SOCWRStateWorkspace, eps_abs::Float64, eps_rel::Float64)
    eps_r = hypot(
        _eps_primal_bus(
            ws.partitions,
            (ws.v2_consensus,),
            (snapshot_v2!,),
            eps_abs,
            eps_rel,
        ),
        _eps_primal_branch(
            ws.partitions,
            (ws.r_consensus, ws.t_consensus),
            ws.interzonal,
            2,
            eps_abs,
            eps_rel,
        ),
    )
    eps_s = hypot(
        _eps_dual_bus(
            (ws.duals_v2,),
            ws.partitions,
            eps_abs,
            eps_rel,
        ),
        _eps_dual_branch(
            (ws.duals_r, ws.duals_t),
            ws.partitions,
            ws.interzonal,
            eps_abs,
            eps_rel,
        ),
    )
    ws.eps_primal = eps_r
    ws.eps_dual = eps_s
    return nothing
end

function compute_tolerances!(ws::SOCWRPowerWorkspace, eps_abs::Float64, eps_rel::Float64)
    eps_r = _eps_primal_branch(
        ws.partitions,
        (ws.p_flow_consensus, ws.q_flow_consensus),
        ws.interzonal,
        2,
        eps_abs,
        eps_rel,
    )
    eps_s = _eps_dual_branch(
        (ws.duals_p, ws.duals_q),
        ws.partitions,
        ws.interzonal,
        eps_abs,
        eps_rel,
    )
    ws.eps_primal = eps_r
    ws.eps_dual = eps_s
    return nothing
end

function compute_tolerances!(ws::SOCWRHybridWorkspace, eps_abs::Float64, eps_rel::Float64)
    eps_r = sqrt(
        _eps_primal_bus(
            ws.partitions,
            (ws.v2_consensus,),
            (snapshot_v2!,),
            eps_abs,
            eps_rel,
        )^2 +
        _eps_primal_branch(
            ws.partitions,
            (ws.p_flow_consensus, ws.q_flow_consensus),
            ws.interzonal,
            2,
            eps_abs,
            eps_rel,
        )^2 +
        _eps_primal_branch(
            ws.partitions,
            (ws.r_consensus, ws.t_consensus),
            ws.interzonal,
            2,
            eps_abs,
            eps_rel,
        )^2,
    )
    eps_s = sqrt(
        _eps_dual_bus(
            (ws.duals_v2,),
            ws.partitions,
            eps_abs,
            eps_rel,
        )^2 +
        _eps_dual_branch(
            (ws.duals_p, ws.duals_q),
            ws.partitions,
            ws.interzonal,
            eps_abs,
            eps_rel,
        )^2 +
        _eps_dual_branch(
            (ws.duals_r, ws.duals_t),
            ws.partitions,
            ws.interzonal,
            eps_abs,
            eps_rel,
        )^2,
    )
    ws.eps_primal = eps_r
    ws.eps_dual = eps_s
    return nothing
end
