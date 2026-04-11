# ============================================================
# ADMM Dual Update — unified for DC · AC · SOCP
# ============================================================
# duals[j] += rho * (x_local[i] - consensus[j])

# ------------------------------------------------------------
# Bus duals: single variable, linear (DC theta, SOCP V²)
# ------------------------------------------------------------
function _update_bus_duals_single!(
    duals::Dict{Int, Vector{Float64}},
    partitions::Partitions,
    consensus::Vector{Float64},
    rho::Float64,
    consensus_bus_idx_map::Dict{Int, Int},
    snapshot_fn,    # (part) -> Vector{Float64}: shared bus values
)
    for part in partitions
        x_local = snapshot_fn(part)
        zone_duals = duals[part.id]
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            j = consensus_bus_idx_map[part.zone_shared_bus_ids[i]]
            zone_duals[j] += rho * (x_local[i] - consensus[j])
        end
    end
    return nothing
end

# ------------------------------------------------------------
# Bus duals: AC (angle-wrapped theta + linear |V|)
# ------------------------------------------------------------
function _update_bus_duals_ac!(
    duals_theta::Dict{Int, Vector{Float64}},
    duals_v_mag::Dict{Int, Vector{Float64}},
    partitions::Partitions{PartitionAC},
    theta_consensus::Vector{Float64},
    v_mag_consensus::Vector{Float64},
    rho_theta::Float64,
    rho_v_mag::Float64,
    consensus_bus_idx_map::Dict{Int, Int},
)
    for part in partitions
        theta_loc = snapshot_theta_shared!(part)
        v_mag_loc = snapshot_v_mag_shared!(part)
        theta_zone_duals = duals_theta[part.id]
        v_mag_zone_duals = duals_v_mag[part.id]
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            j = consensus_bus_idx_map[part.zone_shared_bus_ids[i]]
            theta_zone_duals[j] +=
                rho_theta * _wrap_angle_diff(theta_loc[i] - theta_consensus[j])
            v_mag_zone_duals[j] += rho_v_mag * (v_mag_loc[i] - v_mag_consensus[j])
        end
    end
    return nothing
end

# ------------------------------------------------------------
# Branch duals: single variable (DC P)
# ------------------------------------------------------------
function _update_branch_duals_single!(
    duals::Dict{Int, Dict{Int, Float64}},
    partitions::Partitions,
    consensus::Dict{Int, Float64},
    rho::Float64,
    extract_fn,    # (part, l) -> Float64
)
    for part in partitions
        zone_duals = get!(duals, part.id, Dict{Int, Float64}())
        @inbounds for e in part.zone_edges
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            val = extract_fn(part, l)
            zone_duals[e] =
                get(zone_duals, e, 0.0) + rho * (val - get(consensus, e, 0.0))
        end
    end
    return nothing
end

# ------------------------------------------------------------
# Branch duals: paired variables (SOCP PQ, SOCP RT)
# ------------------------------------------------------------
function _update_branch_duals_paired!(
    duals_a::Dict{Int, Dict{Int, Float64}},
    duals_b::Dict{Int, Dict{Int, Float64}},
    partitions::Partitions,
    consensus_a::Dict{Int, Float64},
    consensus_b::Dict{Int, Float64},
    rho_a::Float64,
    rho_b::Float64,
    extract_fn,    # (part, l) -> (Float64, Float64)
)
    for part in partitions
        a_zone_duals = get!(duals_a, part.id, Dict{Int, Float64}())
        b_zone_duals = get!(duals_b, part.id, Dict{Int, Float64}())
        @inbounds for e in part.zone_edges
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            a, b = extract_fn(part, l)
            a_zone_duals[e] =
                get(a_zone_duals, e, 0.0) + rho_a * (a - get(consensus_a, e, 0.0))
            b_zone_duals[e] =
                get(b_zone_duals, e, 0.0) + rho_b * (b - get(consensus_b, e, 0.0))
        end
    end
    return nothing
end

# ------------------------------------------------------------
# Branch duals: AC P/Q (with isfinite checks for NLP)
# ------------------------------------------------------------
# Specialized variant (not using _update_branch_duals_paired!) because AC flow
# expressions can return NaN (NLP convergence quirks), and we want a bad P
# value to NOT skip a good Q update (and vice versa).
function _update_branch_duals_ac!(
    duals_p::Dict{Int, Dict{Int, Float64}},
    duals_q::Dict{Int, Dict{Int, Float64}},
    partitions::Partitions{PartitionAC},
    p_consensus::Dict{Int, Float64},
    q_consensus::Dict{Int, Float64},
    rho_p::Float64,
    rho_q::Float64,
)
    for part in partitions
        p_zone_duals = get!(duals_p, part.id, Dict{Int, Float64}())
        q_zone_duals = get!(duals_q, part.id, Dict{Int, Float64}())
        @inbounds for e in part.zone_edges
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            p, q = _pq_flow_from_exprs(part, l)
            p_val = JuMP.value(p)
            q_val = JuMP.value(q)
            isfinite(p_val) && (
                p_zone_duals[e] =
                    get(p_zone_duals, e, 0.0) + rho_p * (p_val - get(p_consensus, e, 0.0))
            )
            isfinite(q_val) && (
                q_zone_duals[e] =
                    get(q_zone_duals, e, 0.0) + rho_q * (q_val - get(q_consensus, e, 0.0))
            )
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
        ws.consensus_bus_idx_map,
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
        ws.consensus_bus_idx_map,
        snapshot_theta_shared!,
    )
    _update_branch_duals_single!(
        ws.duals_p_flow,
        ws.partitions,
        ws.p_flow_consensus,
        ws.rho_p[],
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
        ws.consensus_bus_idx_map,
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
        ws.consensus_bus_idx_map,
    )
    _update_branch_duals_ac!(
        ws.duals_p_flow,
        ws.duals_q_flow,
        ws.partitions,
        ws.p_flow_consensus,
        ws.q_flow_consensus,
        ws.rho_p[],
        ws.rho_q[],
    )
    return nothing
end

# ---- SOCP ----

function dual_update!(ws::SOCPStateWorkspace)
    _update_bus_duals_single!(
        ws.duals_v2,
        ws.partitions,
        ws.v2_consensus,
        ws.rho_v2[],
        ws.consensus_bus_idx_map,
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
        _socp_rt_values,
    )
    return nothing
end

function dual_update!(ws::SOCPPowerWorkspace)
    _update_branch_duals_paired!(
        ws.duals_p,
        ws.duals_q,
        ws.partitions,
        ws.p_flow_consensus,
        ws.q_flow_consensus,
        ws.rho_p[],
        ws.rho_q[],
        _socp_pq_values,
    )
    return nothing
end

function dual_update!(ws::SOCPHybridWorkspace)
    _update_bus_duals_single!(
        ws.duals_v2,
        ws.partitions,
        ws.v2_consensus,
        ws.rho_v2[],
        ws.consensus_bus_idx_map,
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
        _socp_pq_values,
    )
    _update_branch_duals_paired!(
        ws.duals_r,
        ws.duals_t,
        ws.partitions,
        ws.r_consensus,
        ws.t_consensus,
        ws.rho_r[],
        ws.rho_t[],
        _socp_rt_values,
    )
    return nothing
end
