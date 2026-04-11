# ============================================================
# ADMM Tolerances (Boyd et al.)
# ------------------------------------------------------------
#   eps_pri  = sqrt(n) * eps_abs + eps_rel * max(|x|, |z|)
#   eps_dual = sqrt(n) * eps_abs + eps_rel * |y|
# ============================================================

# ------------------------------------------------------------
# Bus tolerance helpers
# ------------------------------------------------------------

function _eps_primal_bus(
    partitions::Partitions,
    consensus_vecs,
    consensus_bus_idx_map::Dict{Int, Int},
    snapshot_fns,
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = 0
    norm_x_sq = 0.0
    norm_z_sq = 0.0
    for part in partitions
        bufs = [fn(part) for fn in snapshot_fns]
        pos = part.bus_local_index
        n += length(snapshot_fns) * length(part.zone_shared_bus_ids)
        @inbounds for b in part.zone_shared_bus_ids
            j = consensus_bus_idx_map[b]
            for (k, buf) in enumerate(bufs)
                norm_x_sq += buf[pos[b]]^2
                norm_z_sq += consensus_vecs[k][j]^2
            end
        end
    end
    n = max(n, 1)

    return sqrt(n) * eps_abs + eps_rel * max(sqrt(norm_x_sq), sqrt(norm_z_sq))
end

function _eps_dual_bus(
    duals_list,
    partitions::Partitions,
    zone_ids::Vector{Int},
    consensus_bus_idx_map::Dict{Int, Int},
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = 0
    norm_y_sq = 0.0
    @inbounds for i in eachindex(zone_ids)
        zid = zone_ids[i]
        part = partitions.parts[i]
        for duals in duals_list
            y = duals[zid]
            for b in part.zone_shared_bus_ids
                j = consensus_bus_idx_map[b]
                n += 1
                norm_y_sq += y[j]^2
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
    consensus_dicts,
    interzonal::Vector{Int},
    n_vars_per_edge::Int,
    eps_abs::Float64,
    eps_rel::Float64,
)
    dim = 0
    norm_z_sq = 0.0
    for part in partitions
        @inbounds for e in interzonal
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            dim += n_vars_per_edge
            for d in consensus_dicts
                norm_z_sq += get(d, e, 0.0)^2
            end
        end
    end

    return sqrt(max(1, dim)) * eps_abs + eps_rel * max(sqrt(norm_z_sq), 1e-12)
end

function _eps_dual_branch(
    duals_list,
    zone_ids::Vector{Int},
    interzonal::Vector{Int},
    eps_abs::Float64,
    eps_rel::Float64,
)
    n = 0
    norm_y_sq = 0.0
    @inbounds for zid in zone_ids
        for duals in duals_list
            y = duals[zid]
            for e in interzonal
                haskey(y, e) || continue
                n += 1
                norm_y_sq += get(y, e, 0.0)^2
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
        ws.consensus_bus_idx_map,
        (snapshot_theta!,),
        eps_abs,
        eps_rel,
    )
    eps_s = _eps_dual_bus(
        (ws.duals_theta,),
        ws.partitions,
        ws.id,
        ws.consensus_bus_idx_map,
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
        ws.id,
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
            ws.consensus_bus_idx_map,
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
            ws.id,
            ws.consensus_bus_idx_map,
            eps_abs,
            eps_rel,
        ),
        _eps_dual_branch((ws.duals_p_flow,), ws.id, ws.interzonal, eps_abs, eps_rel),
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
        ws.consensus_bus_idx_map,
        (snapshot_theta!, snapshot_v_mag!),
        eps_abs,
        eps_rel,
    )
    eps_s = _eps_dual_bus(
        (ws.duals_theta, ws.duals_v_mag),
        ws.partitions,
        ws.id,
        ws.consensus_bus_idx_map,
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
        ws.id,
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
            ws.consensus_bus_idx_map,
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
            ws.id,
            ws.consensus_bus_idx_map,
            eps_abs,
            eps_rel,
        ),
        _eps_dual_branch(
            (ws.duals_p_flow, ws.duals_q_flow),
            ws.id,
            ws.interzonal,
            eps_abs,
            eps_rel,
        ),
    )
    ws.eps_primal = eps_r
    ws.eps_dual = eps_s
    return nothing
end

# ---- SOCP ----

function compute_tolerances!(ws::SOCPStateWorkspace, eps_abs::Float64, eps_rel::Float64)
    eps_r = hypot(
        _eps_primal_bus(
            ws.partitions,
            (ws.v2_consensus,),
            ws.consensus_bus_idx_map,
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
            ws.id,
            ws.consensus_bus_idx_map,
            eps_abs,
            eps_rel,
        ),
        _eps_dual_branch(
            (ws.duals_r, ws.duals_t),
            ws.id,
            ws.interzonal,
            eps_abs,
            eps_rel,
        ),
    )
    ws.eps_primal = eps_r
    ws.eps_dual = eps_s
    return nothing
end

function compute_tolerances!(ws::SOCPPowerWorkspace, eps_abs::Float64, eps_rel::Float64)
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
        ws.id,
        ws.interzonal,
        eps_abs,
        eps_rel,
    )
    ws.eps_primal = eps_r
    ws.eps_dual = eps_s
    return nothing
end

function compute_tolerances!(ws::SOCPHybridWorkspace, eps_abs::Float64, eps_rel::Float64)
    eps_r = sqrt(
        _eps_primal_bus(
            ws.partitions,
            (ws.v2_consensus,),
            ws.consensus_bus_idx_map,
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
            ws.id,
            ws.consensus_bus_idx_map,
            eps_abs,
            eps_rel,
        )^2 +
        _eps_dual_branch(
            (ws.duals_p, ws.duals_q),
            ws.id,
            ws.interzonal,
            eps_abs,
            eps_rel,
        )^2 +
        _eps_dual_branch(
            (ws.duals_r, ws.duals_t),
            ws.id,
            ws.interzonal,
            eps_abs,
            eps_rel,
        )^2,
    )
    ws.eps_primal = eps_r
    ws.eps_dual = eps_s
    return nothing
end
