# ============================================================
# ADMM Residuals & Tolerances — unified for DC · AC · SOCP
# ============================================================

# ############################################################
#  PRIMAL RESIDUALS
# ############################################################

# ------------------------------------------------------------
# Bus primal residual: single variable (DC theta, SOCP V²)
# ------------------------------------------------------------
function _primal_residual_bus_single(
    partitions::Partitions,
    consensus::Vector{Float64},
    consensus_bus_idx_map::Dict{Int, Int},
    snapshot_fn,    # (part) -> Vector{Float64}: shared bus values
)
    acc = 0.0
    for part in partitions
        buf = snapshot_fn(part)
        pos = part.bus_local_index
        @inbounds for b in part.zone_shared_bus_ids
            acc += (buf[pos[b]] - consensus[consensus_bus_idx_map[b]])^2
        end
    end
    return sqrt(acc)
end

# ------------------------------------------------------------
# Bus primal residual: AC components (separate theta / v_mag)
# ------------------------------------------------------------
function _primal_residual_bus_ac_components(
    partitions::Partitions{PartitionAC},
    theta_consensus::Vector{Float64},
    v_mag_consensus::Vector{Float64},
    consensus_bus_idx_map::Dict{Int, Int},
)
    r_theta = 0.0
    r_v_mag = 0.0
    for part in partitions
        theta_buf = snapshot_theta!(part)
        v_mag_buf = snapshot_v_mag!(part)
        pos = part.bus_local_index
        @inbounds for b in part.zone_shared_bus_ids
            j = consensus_bus_idx_map[b]
            r_theta += _wrap_angle_diff(theta_buf[pos[b]] - theta_consensus[j])^2
            r_v_mag += (v_mag_buf[pos[b]] - v_mag_consensus[j])^2
        end
    end
    return sqrt(r_theta), sqrt(r_v_mag)
end

# ------------------------------------------------------------
# Branch primal residual: single variable (DC P)
# ------------------------------------------------------------
function _primal_residual_branch_single(
    partitions::Partitions,
    consensus::Dict{Int, Float64},
    interzonal::Vector{Int},
    extract_fn,    # (part, l) -> Float64
)
    acc = 0.0
    for part in partitions
        @inbounds for e in interzonal
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            acc += (extract_fn(part, l) - get(consensus, e, 0.0))^2
        end
    end
    return sqrt(acc)
end

# ------------------------------------------------------------
# Branch primal residual: paired components (SOCP PQ, SOCP RT)
# ------------------------------------------------------------
function _primal_residual_branch_paired_components(
    partitions::Partitions,
    consensus_a::Dict{Int, Float64},
    consensus_b::Dict{Int, Float64},
    interzonal::Vector{Int},
    extract_fn,    # (part, l) -> (Float64, Float64)
)
    acc_a = 0.0
    acc_b = 0.0
    for part in partitions
        @inbounds for e in interzonal
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            a, b = extract_fn(part, l)
            acc_a += (a - get(consensus_a, e, 0.0))^2
            acc_b += (b - get(consensus_b, e, 0.0))^2
        end
    end
    return sqrt(acc_a), sqrt(acc_b)
end

# ------------------------------------------------------------
# Branch primal residual: AC P/Q components (with isfinite checks for NLP)
# ------------------------------------------------------------
# AC flow expressions can return NaN (NLP convergence quirks), so bad values
# are skipped instead of poisoning the accumulator with NaN.
function _primal_residual_branch_ac_components(
    partitions::Partitions{PartitionAC},
    p_consensus::Dict{Int, Float64},
    q_consensus::Dict{Int, Float64},
    interzonal::Vector{Int},
)
    acc_p = 0.0
    acc_q = 0.0
    for part in partitions
        @inbounds for e in interzonal
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            p, q = _pq_flow_from_exprs(part, l; t = 1)
            p_val = JuMP.value(p)
            q_val = JuMP.value(q)
            isfinite(p_val) && (acc_p += (p_val - get(p_consensus, e, 0.0))^2)
            isfinite(q_val) && (acc_q += (q_val - get(q_consensus, e, 0.0))^2)
        end
    end
    return sqrt(acc_p), sqrt(acc_q)
end

# ############################################################
#  DUAL RESIDUALS (Boyd canonical — expanded consensus form)
# ############################################################
#
# In consensus ADMM (Boyd §5.2), each agent has its own dual y_i for each
# shared variable b. The canonical dual residual stacks N copies of
# (z^{k+1} - z^k):
#
#   ||s||² = Σ_zone Σ_{b ∈ zone} (ρ·(z^b_k - z^b_{k-1}))²
#          = ρ² · Σ_b m_b · (z^b_k - z^b_{k-1})²
#
# where m_b is the multiplicity (# zones sharing b). We iterate over
# (zone, var) pairs so m_b arises naturally and matches ε_dual which
# already sums over all pairs.

# ------------------------------------------------------------
# Bus dual residual: single vector (DC theta, SOCP V²)
# ------------------------------------------------------------
function _dual_residual_bus_single(
    partitions::Partitions,
    consensus::Vector{Float64},
    consensus_prev::Vector{Float64},
    consensus_bus_idx_map::Dict{Int, Int},
    rho::Float64,
)
    s = 0.0
    for part in partitions
        @inbounds for b in part.zone_shared_bus_ids
            j = consensus_bus_idx_map[b]
            s += (consensus[j] - consensus_prev[j])^2
        end
    end
    return rho * sqrt(s)
end

# ------------------------------------------------------------
# Bus dual residual: AC components (separate theta / v_mag)
# ------------------------------------------------------------
# The theta component uses angle wrapping to handle 2π periodicity.
function _dual_residual_bus_ac_components(
    partitions::Partitions{PartitionAC},
    theta_consensus::Vector{Float64},
    v_mag_consensus::Vector{Float64},
    theta_consensus_prev::Vector{Float64},
    v_mag_consensus_prev::Vector{Float64},
    consensus_bus_idx_map::Dict{Int, Int},
    rho_theta::Float64,
    rho_v_mag::Float64,
)
    s_theta = 0.0
    s_v_mag = 0.0
    for part in partitions
        @inbounds for b in part.zone_shared_bus_ids
            j = consensus_bus_idx_map[b]
            s_theta += _wrap_angle_diff(theta_consensus[j] - theta_consensus_prev[j])^2
            s_v_mag += (v_mag_consensus[j] - v_mag_consensus_prev[j])^2
        end
    end
    return rho_theta * sqrt(s_theta), rho_v_mag * sqrt(s_v_mag)
end

# ------------------------------------------------------------
# Branch dual residual: single variable (DC P)
# ------------------------------------------------------------
function _dual_residual_branch_single(
    partitions::Partitions,
    consensus::Dict{Int, Float64},
    consensus_prev::Dict{Int, Float64},
    interzonal::Vector{Int},
    rho::Float64,
)
    s = 0.0
    for part in partitions
        @inbounds for e in interzonal
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            s += (get(consensus, e, 0.0) - get(consensus_prev, e, 0.0))^2
        end
    end
    return rho * sqrt(s)
end

# ------------------------------------------------------------
# Branch dual residual: paired components (SOCP PQ, SOCP RT, AC P/Q)
# ------------------------------------------------------------
# Each block is rho-scaled independently (Boyd per-block form).
function _dual_residual_branch_paired_components(
    partitions::Partitions,
    consensus_a::Dict{Int, Float64},
    consensus_b::Dict{Int, Float64},
    consensus_a_prev::Dict{Int, Float64},
    consensus_b_prev::Dict{Int, Float64},
    interzonal::Vector{Int},
    rho_a::Float64,
    rho_b::Float64,
)
    s_a = 0.0
    s_b = 0.0
    for part in partitions
        @inbounds for e in interzonal
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            diff_a = get(consensus_a, e, 0.0) - get(consensus_a_prev, e, 0.0)
            diff_b = get(consensus_b, e, 0.0) - get(consensus_b_prev, e, 0.0)
            s_a += diff_a^2
            s_b += diff_b^2
        end
    end
    return rho_a * sqrt(s_a), rho_b * sqrt(s_b)
end

# ############################################################
#  DISPATCH: compute_residuals!(ws) → (r, s)
# ############################################################

# ---- DC ----

function compute_residuals!(ws::DCStateWorkspace)
    r = _primal_residual_bus_single(
        ws.partitions,
        ws.theta_consensus,
        ws.consensus_bus_idx_map,
        snapshot_theta!,
    )
    s = _dual_residual_bus_single(
        ws.partitions,
        ws.theta_consensus,
        ws.theta_consensus_prev,
        ws.consensus_bus_idx_map,
        ws.rho_theta[],
    )
    ws.r_primal = r
    ws.r_dual = s
    return nothing
end

function compute_residuals!(ws::DCPowerWorkspace)
    r = _primal_residual_branch_single(
        ws.partitions,
        ws.p_flow_consensus,
        ws.interzonal,
        _dc_flow_value,
    )
    s = _dual_residual_branch_single(
        ws.partitions,
        ws.p_flow_consensus,
        ws.p_flow_consensus_prev,
        ws.interzonal,
        ws.rho_p[],
    )
    ws.r_primal = r
    ws.r_dual = s
    return nothing
end

function compute_residuals!(ws::DCHybridWorkspace)
    # Components (used by rho_update! without recomputation)
    ws.r_primal_bus = _primal_residual_bus_single(
        ws.partitions,
        ws.theta_consensus,
        ws.consensus_bus_idx_map,
        snapshot_theta!,
    )
    ws.r_dual_bus = _dual_residual_bus_single(
        ws.partitions,
        ws.theta_consensus,
        ws.theta_consensus_prev,
        ws.consensus_bus_idx_map,
        ws.rho_theta[],
    )
    ws.r_primal_branch = _primal_residual_branch_single(
        ws.partitions,
        ws.p_flow_consensus,
        ws.interzonal,
        _dc_flow_value,
    )
    ws.r_dual_branch = _dual_residual_branch_single(
        ws.partitions,
        ws.p_flow_consensus,
        ws.p_flow_consensus_prev,
        ws.interzonal,
        ws.rho_p[],
    )

    # Combined (for convergence check, history, print)
    ws.r_primal = hypot(ws.r_primal_bus, ws.r_primal_branch)
    ws.r_dual = hypot(ws.r_dual_bus, ws.r_dual_branch)
    return nothing
end

# ---- AC ----

function compute_residuals!(ws::ACStateWorkspace)
    # Components (used by rho_update! without recomputation)
    ws.r_primal_theta, ws.r_primal_v_mag = _primal_residual_bus_ac_components(
        ws.partitions,
        ws.theta_consensus,
        ws.v_mag_consensus,
        ws.consensus_bus_idx_map,
    )
    ws.r_dual_theta, ws.r_dual_v_mag = _dual_residual_bus_ac_components(
        ws.partitions,
        ws.theta_consensus, ws.v_mag_consensus,
        ws.theta_consensus_prev, ws.v_mag_consensus_prev,
        ws.consensus_bus_idx_map,
        ws.rho_theta[], ws.rho_v[],
    )

    # Combined (for convergence check, history, print)
    ws.r_primal = hypot(ws.r_primal_theta, ws.r_primal_v_mag)
    ws.r_dual = hypot(ws.r_dual_theta, ws.r_dual_v_mag)
    return nothing
end

function compute_residuals!(ws::ACPowerWorkspace)
    # Components (used by rho_update! without recomputation)
    ws.r_primal_p, ws.r_primal_q = _primal_residual_branch_ac_components(
        ws.partitions,
        ws.p_flow_consensus,
        ws.q_flow_consensus,
        ws.interzonal,
    )
    ws.r_dual_p, ws.r_dual_q = _dual_residual_branch_paired_components(
        ws.partitions,
        ws.p_flow_consensus, ws.q_flow_consensus,
        ws.p_flow_consensus_prev, ws.q_flow_consensus_prev,
        ws.interzonal, ws.rho_p[], ws.rho_q[],
    )

    # Combined (for convergence check, history, print)
    ws.r_primal = hypot(ws.r_primal_p, ws.r_primal_q)
    ws.r_dual = hypot(ws.r_dual_p, ws.r_dual_q)
    return nothing
end

function compute_residuals!(ws::ACHybridWorkspace)
    # Bus components (used by rho_update! without recomputation)
    ws.r_primal_theta, ws.r_primal_v_mag = _primal_residual_bus_ac_components(
        ws.partitions,
        ws.theta_consensus,
        ws.v_mag_consensus,
        ws.consensus_bus_idx_map,
    )
    ws.r_dual_theta, ws.r_dual_v_mag = _dual_residual_bus_ac_components(
        ws.partitions,
        ws.theta_consensus, ws.v_mag_consensus,
        ws.theta_consensus_prev, ws.v_mag_consensus_prev,
        ws.consensus_bus_idx_map,
        ws.rho_theta[], ws.rho_v[],
    )

    # Branch components
    ws.r_primal_p, ws.r_primal_q = _primal_residual_branch_ac_components(
        ws.partitions,
        ws.p_flow_consensus,
        ws.q_flow_consensus,
        ws.interzonal,
    )
    ws.r_dual_p, ws.r_dual_q = _dual_residual_branch_paired_components(
        ws.partitions,
        ws.p_flow_consensus, ws.q_flow_consensus,
        ws.p_flow_consensus_prev, ws.q_flow_consensus_prev,
        ws.interzonal, ws.rho_p[], ws.rho_q[],
    )

    # Combined (for convergence check, history, print)
    ws.r_primal = sqrt(
        ws.r_primal_theta^2 + ws.r_primal_v_mag^2 +
        ws.r_primal_p^2 + ws.r_primal_q^2,
    )
    ws.r_dual = sqrt(
        ws.r_dual_theta^2 + ws.r_dual_v_mag^2 +
        ws.r_dual_p^2 + ws.r_dual_q^2,
    )
    return nothing
end

# ---- SOCP ----

function compute_residuals!(ws::SOCPStateWorkspace)
    # Components (used by rho_update! without recomputation)
    ws.r_primal_bus = _primal_residual_bus_single(
        ws.partitions,
        ws.v2_consensus,
        ws.consensus_bus_idx_map,
        snapshot_v2!,
    )
    ws.r_dual_bus = _dual_residual_bus_single(
        ws.partitions,
        ws.v2_consensus,
        ws.v2_consensus_prev,
        ws.consensus_bus_idx_map,
        ws.rho_v2[],
    )
    ws.r_primal_r, ws.r_primal_t = _primal_residual_branch_paired_components(
        ws.partitions,
        ws.r_consensus,
        ws.t_consensus,
        ws.interzonal,
        _socp_rt_values,
    )
    ws.r_dual_r, ws.r_dual_t = _dual_residual_branch_paired_components(
        ws.partitions,
        ws.r_consensus, ws.t_consensus,
        ws.r_consensus_prev, ws.t_consensus_prev,
        ws.interzonal, ws.rho_r[], ws.rho_t[],
    )

    # Combined (for convergence check, history, print)
    ws.r_primal = sqrt(ws.r_primal_bus^2 + ws.r_primal_r^2 + ws.r_primal_t^2)
    ws.r_dual = sqrt(ws.r_dual_bus^2 + ws.r_dual_r^2 + ws.r_dual_t^2)
    return nothing
end

function compute_residuals!(ws::SOCPPowerWorkspace)
    # Components (used by rho_update! without recomputation)
    ws.r_primal_p, ws.r_primal_q = _primal_residual_branch_paired_components(
        ws.partitions,
        ws.p_flow_consensus,
        ws.q_flow_consensus,
        ws.interzonal,
        _socp_pq_values,
    )
    ws.r_dual_p, ws.r_dual_q = _dual_residual_branch_paired_components(
        ws.partitions,
        ws.p_flow_consensus, ws.q_flow_consensus,
        ws.p_flow_consensus_prev, ws.q_flow_consensus_prev,
        ws.interzonal, ws.rho_p[], ws.rho_q[],
    )

    # Combined (for convergence check, history, print)
    ws.r_primal = hypot(ws.r_primal_p, ws.r_primal_q)
    ws.r_dual = hypot(ws.r_dual_p, ws.r_dual_q)
    return nothing
end

function compute_residuals!(ws::SOCPHybridWorkspace)
    # Components (used by rho_update! without recomputation)
    ws.r_primal_bus = _primal_residual_bus_single(
        ws.partitions,
        ws.v2_consensus,
        ws.consensus_bus_idx_map,
        snapshot_v2!,
    )
    ws.r_dual_bus = _dual_residual_bus_single(
        ws.partitions,
        ws.v2_consensus,
        ws.v2_consensus_prev,
        ws.consensus_bus_idx_map,
        ws.rho_v2[],
    )
    ws.r_primal_p, ws.r_primal_q = _primal_residual_branch_paired_components(
        ws.partitions,
        ws.p_flow_consensus,
        ws.q_flow_consensus,
        ws.interzonal,
        _socp_pq_values,
    )
    ws.r_dual_p, ws.r_dual_q = _dual_residual_branch_paired_components(
        ws.partitions,
        ws.p_flow_consensus, ws.q_flow_consensus,
        ws.p_flow_consensus_prev, ws.q_flow_consensus_prev,
        ws.interzonal, ws.rho_p[], ws.rho_q[],
    )
    ws.r_primal_r, ws.r_primal_t = _primal_residual_branch_paired_components(
        ws.partitions,
        ws.r_consensus,
        ws.t_consensus,
        ws.interzonal,
        _socp_rt_values,
    )
    ws.r_dual_r, ws.r_dual_t = _dual_residual_branch_paired_components(
        ws.partitions,
        ws.r_consensus, ws.t_consensus,
        ws.r_consensus_prev, ws.t_consensus_prev,
        ws.interzonal, ws.rho_r[], ws.rho_t[],
    )

    # Combined (for convergence check, history, print)
    ws.r_primal = sqrt(
        ws.r_primal_bus^2 +
        ws.r_primal_p^2 + ws.r_primal_q^2 +
        ws.r_primal_r^2 + ws.r_primal_t^2,
    )
    ws.r_dual = sqrt(
        ws.r_dual_bus^2 +
        ws.r_dual_p^2 + ws.r_dual_q^2 +
        ws.r_dual_r^2 + ws.r_dual_t^2,
    )
    return nothing
end
