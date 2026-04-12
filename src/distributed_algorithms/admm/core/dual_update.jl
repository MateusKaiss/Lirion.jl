# ============================================================
# ADMM Dual Update — unified for DC · AC · SOCWR
# ============================================================

# ------------------------------------------------------------
# Bus duals: single variable, linear (DC theta, SOCWR V²)
# ------------------------------------------------------------
function _update_bus_duals_single!(
    duals::Matrix{Float64},
    partitions::Partitions,
    consensus::Vector{Float64},
    rho::Float64,
    snapshot_fn,    # (part) -> Vector{Float64}: shared-bus values
)
    for part in partitions
        x_local = snapshot_fn(part)
        p = part.position
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            j = part.shared_bus_consensus_idx[i]
            duals[j, p] += rho * (x_local[i] - consensus[j])
        end
    end
    return nothing
end

# ------------------------------------------------------------
# Bus duals: AC (angle-wrapped theta + linear |V|)
# ------------------------------------------------------------
function _update_bus_duals_ac!(
    duals_theta::Matrix{Float64},
    duals_v_mag::Matrix{Float64},
    partitions::Partitions{PartitionAC},
    theta_consensus::Vector{Float64},
    v_mag_consensus::Vector{Float64},
    rho_theta::Float64,
    rho_v_mag::Float64,
)
    for part in partitions
        theta_loc = snapshot_theta_shared!(part)
        v_mag_loc = snapshot_v_mag_shared!(part)
        p = part.position
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            j = part.shared_bus_consensus_idx[i]
            duals_theta[j, p] +=
                rho_theta * _wrap_angle_diff(theta_loc[i] - theta_consensus[j])
            duals_v_mag[j, p] += rho_v_mag * (v_mag_loc[i] - v_mag_consensus[j])
        end
    end
    return nothing
end

# ------------------------------------------------------------
# Branch duals: single variable (DC P)
# ------------------------------------------------------------
function _update_branch_duals_single!(
    duals::Matrix{Float64},
    partitions::Partitions,
    consensus::Vector{Float64},
    rho::Float64,
    interzonal::Vector{Int},
    extract_fn,    # (part, l) -> Float64
)
    for part in partitions
        p = part.position
        @inbounds for (k, e) in enumerate(interzonal)
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            val = extract_fn(part, l)
            duals[k, p] += rho * (val - consensus[k])
        end
    end
    return nothing
end

# ------------------------------------------------------------
# Branch duals: paired variables (SOCWR PQ, SOCWR RT)
# ------------------------------------------------------------
function _update_branch_duals_paired!(
    duals_a::Matrix{Float64},
    duals_b::Matrix{Float64},
    partitions::Partitions,
    consensus_a::Vector{Float64},
    consensus_b::Vector{Float64},
    rho_a::Float64,
    rho_b::Float64,
    interzonal::Vector{Int},
    extract_fn,    # (part, l) -> (Float64, Float64)
)
    for part in partitions
        p = part.position
        @inbounds for (k, e) in enumerate(interzonal)
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            a, b = extract_fn(part, l)
            duals_a[k, p] += rho_a * (a - consensus_a[k])
            duals_b[k, p] += rho_b * (b - consensus_b[k])
        end
    end
    return nothing
end

# ------------------------------------------------------------
# Branch duals: AC P/Q (with isfinite checks for NLP)
# ------------------------------------------------------------
# Specialized variant because AC flow expressions can return NaN; a bad
# P must NOT skip a good Q update (and vice versa).
function _update_branch_duals_ac!(
    duals_p::Matrix{Float64},
    duals_q::Matrix{Float64},
    partitions::Partitions{PartitionAC},
    p_consensus::Vector{Float64},
    q_consensus::Vector{Float64},
    rho_p::Float64,
    rho_q::Float64,
    interzonal::Vector{Int},
)
    for part in partitions
        p = part.position
        @inbounds for (k, e) in enumerate(interzonal)
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            p_expr, q_expr = _pq_flow_from_exprs(part, l)
            p_val = JuMP.value(p_expr)
            q_val = JuMP.value(q_expr)
            isfinite(p_val) &&
                (duals_p[k, p] += rho_p * (p_val - p_consensus[k]))
            isfinite(q_val) &&
                (duals_q[k, p] += rho_q * (q_val - q_consensus[k]))
        end
    end
    return nothing
end

# ============================================================
# Dispatch: dual_update!(ws)
# ============================================================

# ---- DC ----

function dual_update!(ws::DCStateWorkspace)
    _update_bus_duals_single!(
        ws.duals_theta,
        ws.partitions,
        ws.theta_consensus,
        ws.rho_theta[],
        snapshot_theta_shared!,
    )
    return nothing
end

function dual_update!(ws::DCPowerWorkspace)
    _update_branch_duals_single!(
        ws.duals_p_flow,
        ws.partitions,
        ws.p_flow_consensus,
        ws.rho_p[],
        ws.interzonal,
        _dc_flow_value,
    )
    return nothing
end

function dual_update!(ws::DCHybridWorkspace)
    _update_bus_duals_single!(
        ws.duals_theta,
        ws.partitions,
        ws.theta_consensus,
        ws.rho_theta[],
        snapshot_theta_shared!,
    )
    _update_branch_duals_single!(
        ws.duals_p_flow,
        ws.partitions,
        ws.p_flow_consensus,
        ws.rho_p[],
        ws.interzonal,
        _dc_flow_value,
    )
    return nothing
end

# ---- AC ----

function dual_update!(ws::ACStateWorkspace)
    _update_bus_duals_ac!(
        ws.duals_theta,
        ws.duals_v_mag,
        ws.partitions,
        ws.theta_consensus,
        ws.v_mag_consensus,
        ws.rho_theta[],
        ws.rho_v[],
    )
    return nothing
end

function dual_update!(ws::ACPowerWorkspace)
    _update_branch_duals_ac!(
        ws.duals_p_flow,
        ws.duals_q_flow,
        ws.partitions,
        ws.p_flow_consensus,
        ws.q_flow_consensus,
        ws.rho_p[],
        ws.rho_q[],
        ws.interzonal,
    )
    return nothing
end

function dual_update!(ws::ACHybridWorkspace)
    _update_bus_duals_ac!(
        ws.duals_theta,
        ws.duals_v_mag,
        ws.partitions,
        ws.theta_consensus,
        ws.v_mag_consensus,
        ws.rho_theta[],
        ws.rho_v[],
    )
    _update_branch_duals_ac!(
        ws.duals_p_flow,
        ws.duals_q_flow,
        ws.partitions,
        ws.p_flow_consensus,
        ws.q_flow_consensus,
        ws.rho_p[],
        ws.rho_q[],
        ws.interzonal,
    )
    return nothing
end

# ---- SOCWR ----

function dual_update!(ws::SOCWRStateWorkspace)
    _update_bus_duals_single!(
        ws.duals_v2,
        ws.partitions,
        ws.v2_consensus,
        ws.rho_v2[],
        snapshot_v2_shared!,
    )
    _update_branch_duals_paired!(
        ws.duals_r,
        ws.duals_t,
        ws.partitions,
        ws.r_consensus,
        ws.t_consensus,
        ws.rho_r[],
        ws.rho_t[],
        ws.interzonal,
        _socwr_rt_values,
    )
    return nothing
end

function dual_update!(ws::SOCWRPowerWorkspace)
    _update_branch_duals_paired!(
        ws.duals_p,
        ws.duals_q,
        ws.partitions,
        ws.p_flow_consensus,
        ws.q_flow_consensus,
        ws.rho_p[],
        ws.rho_q[],
        ws.interzonal,
        _socwr_pq_values,
    )
    return nothing
end

function dual_update!(ws::SOCWRHybridWorkspace)
    _update_bus_duals_single!(
        ws.duals_v2,
        ws.partitions,
        ws.v2_consensus,
        ws.rho_v2[],
        snapshot_v2_shared!,
    )
    _update_branch_duals_paired!(
        ws.duals_p,
        ws.duals_q,
        ws.partitions,
        ws.p_flow_consensus,
        ws.q_flow_consensus,
        ws.rho_p[],
        ws.rho_q[],
        ws.interzonal,
        _socwr_pq_values,
    )
    _update_branch_duals_paired!(
        ws.duals_r,
        ws.duals_t,
        ws.partitions,
        ws.r_consensus,
        ws.t_consensus,
        ws.rho_r[],
        ws.rho_t[],
        ws.interzonal,
        _socwr_rt_values,
    )
    return nothing
end
