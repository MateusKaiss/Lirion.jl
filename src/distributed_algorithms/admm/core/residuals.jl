# ============================================================
# ADMM Residuals — unified for DC · AC · SOCWR
# ============================================================

# ------------------------------------------------------------
# Bus primal residual: single variable (DC theta, SOCWR V²)
# ------------------------------------------------------------
function _primal_residual_bus_single(
    partitions::Partitions,
    consensus::Vector{Float64},
    snapshot_fn,    # (part) -> Vector{Float64}: full-bus buffer
)
    acc = 0.0
    for part in partitions
        buf = snapshot_fn(part)
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            li = part.shared_bus_local_idx[i]
            cj = part.shared_bus_consensus_idx[i]
            acc += (buf[li] - consensus[cj])^2
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
)
    r_theta = 0.0
    r_v_mag = 0.0
    for part in partitions
        theta_buf = snapshot_theta!(part)
        v_mag_buf = snapshot_v_mag!(part)
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            li = part.shared_bus_local_idx[i]
            cj = part.shared_bus_consensus_idx[i]
            r_theta += _wrap_angle_diff(theta_buf[li] - theta_consensus[cj])^2
            r_v_mag += (v_mag_buf[li] - v_mag_consensus[cj])^2
        end
    end
    return sqrt(r_theta), sqrt(r_v_mag)
end

# ------------------------------------------------------------
# Branch primal residual: single variable (DC P)
# ------------------------------------------------------------
function _primal_residual_branch_single(
    partitions::Partitions,
    consensus::Vector{Float64},
    interzonal::Vector{Int},
    extract_fn,    # (part, l) -> Float64
)
    acc = 0.0
    for part in partitions
        @inbounds for (k, e) in enumerate(interzonal)
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            acc += (extract_fn(part, l) - consensus[k])^2
        end
    end
    return sqrt(acc)
end

# ------------------------------------------------------------
# Branch primal residual: paired components (SOCWR PQ, SOCWR RT)
# ------------------------------------------------------------
function _primal_residual_branch_paired_components(
    partitions::Partitions,
    consensus_a::Vector{Float64},
    consensus_b::Vector{Float64},
    interzonal::Vector{Int},
    extract_fn,    # (part, l) -> (Float64, Float64)
)
    acc_a = 0.0
    acc_b = 0.0
    for part in partitions
        @inbounds for (k, e) in enumerate(interzonal)
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            a, b = extract_fn(part, l)
            acc_a += (a - consensus_a[k])^2
            acc_b += (b - consensus_b[k])^2
        end
    end
    return sqrt(acc_a), sqrt(acc_b)
end

# ------------------------------------------------------------
# Branch primal residual: AC P/Q components (with isfinite for NLP)
# ------------------------------------------------------------
function _primal_residual_branch_ac_components(
    partitions::Partitions{PartitionAC},
    p_consensus::Vector{Float64},
    q_consensus::Vector{Float64},
    interzonal::Vector{Int},
)
    acc_p = 0.0
    acc_q = 0.0
    for part in partitions
        @inbounds for (k, e) in enumerate(interzonal)
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            p, q = _pq_flow_from_exprs(part, l; t = 1)
            p_val = JuMP.value(p)
            q_val = JuMP.value(q)
            isfinite(p_val) && (acc_p += (p_val - p_consensus[k])^2)
            isfinite(q_val) && (acc_q += (q_val - q_consensus[k])^2)
        end
    end
    return sqrt(acc_p), sqrt(acc_q)
end

# ------------------------------------------------------------
# Bus dual residual: single vector (DC theta, SOCWR V²)
# ------------------------------------------------------------
function _dual_residual_bus_single(
    partitions::Partitions,
    consensus::Vector{Float64},
    consensus_prev::Vector{Float64},
    rho::Float64,
)
    s = 0.0
    for part in partitions
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            j = part.shared_bus_consensus_idx[i]
            s += (consensus[j] - consensus_prev[j])^2
        end
    end
    return rho * sqrt(s)
end

# ------------------------------------------------------------
# Bus dual residual: AC components (separate theta / v_mag)
# ------------------------------------------------------------
# Theta uses angle wrapping to handle 2π periodicity.
function _dual_residual_bus_ac_components(
    partitions::Partitions{PartitionAC},
    theta_consensus::Vector{Float64},
    v_mag_consensus::Vector{Float64},
    theta_consensus_prev::Vector{Float64},
    v_mag_consensus_prev::Vector{Float64},
    rho_theta::Float64,
    rho_v_mag::Float64,
)
    s_theta = 0.0
    s_v_mag = 0.0
    for part in partitions
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            j = part.shared_bus_consensus_idx[i]
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
    consensus::Vector{Float64},
    consensus_prev::Vector{Float64},
    interzonal::Vector{Int},
    rho::Float64,
)
    s = 0.0
    for part in partitions
        @inbounds for (k, e) in enumerate(interzonal)
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            s += (consensus[k] - consensus_prev[k])^2
        end
    end
    return rho * sqrt(s)
end

# ------------------------------------------------------------
# Branch dual residual: paired components (SOCWR PQ, SOCWR RT, AC P/Q)
# ------------------------------------------------------------
function _dual_residual_branch_paired_components(
    partitions::Partitions,
    consensus_a::Vector{Float64},
    consensus_b::Vector{Float64},
    consensus_a_prev::Vector{Float64},
    consensus_b_prev::Vector{Float64},
    interzonal::Vector{Int},
    rho_a::Float64,
    rho_b::Float64,
)
    s_a = 0.0
    s_b = 0.0
    for part in partitions
        @inbounds for (k, e) in enumerate(interzonal)
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            diff_a = consensus_a[k] - consensus_a_prev[k]
            diff_b = consensus_b[k] - consensus_b_prev[k]
            s_a += diff_a^2
            s_b += diff_b^2
        end
    end
    return rho_a * sqrt(s_a), rho_b * sqrt(s_b)
end

# ---- DC ----

function compute_residuals!(ws::DCStateWorkspace)
    r = _primal_residual_bus_single(
        ws.partitions,
        ws.vars.theta.consensus,
        snapshot_theta!,
    )
    s = _dual_residual_bus_single(
        ws.partitions,
        ws.vars.theta.consensus,
        ws.vars.theta.consensus_prev,
        ws.vars.theta.rho[],
    )
    ws.r_primal = r
    ws.r_dual = s
    return nothing
end

function compute_residuals!(ws::DCPowerWorkspace)
    r = _primal_residual_branch_single(
        ws.partitions,
        ws.vars.p_flow.consensus,
        ws.interzonal,
        _dc_flow_value,
    )
    s = _dual_residual_branch_single(
        ws.partitions,
        ws.vars.p_flow.consensus,
        ws.vars.p_flow.consensus_prev,
        ws.interzonal,
        ws.vars.p_flow.rho[],
    )
    ws.r_primal = r
    ws.r_dual = s
    return nothing
end

function compute_residuals!(ws::DCHybridWorkspace)
    ws.r_primal_bus = _primal_residual_bus_single(
        ws.partitions,
        ws.vars.theta.consensus,
        snapshot_theta!,
    )
    ws.r_dual_bus = _dual_residual_bus_single(
        ws.partitions,
        ws.vars.theta.consensus,
        ws.vars.theta.consensus_prev,
        ws.vars.theta.rho[],
    )
    ws.r_primal_branch = _primal_residual_branch_single(
        ws.partitions,
        ws.vars.p_flow.consensus,
        ws.interzonal,
        _dc_flow_value,
    )
    ws.r_dual_branch = _dual_residual_branch_single(
        ws.partitions,
        ws.vars.p_flow.consensus,
        ws.vars.p_flow.consensus_prev,
        ws.interzonal,
        ws.vars.p_flow.rho[],
    )

    ws.r_primal = hypot(ws.r_primal_bus, ws.r_primal_branch)
    ws.r_dual = hypot(ws.r_dual_bus, ws.r_dual_branch)
    return nothing
end

# ---- AC ----

function compute_residuals!(ws::ACStateWorkspace)
    ws.r_primal_theta, ws.r_primal_v_mag = _primal_residual_bus_ac_components(
        ws.partitions,
        ws.vars.theta.consensus,
        ws.vars.v_mag.consensus,
    )
    ws.r_dual_theta, ws.r_dual_v_mag = _dual_residual_bus_ac_components(
        ws.partitions,
        ws.vars.theta.consensus, ws.vars.v_mag.consensus,
        ws.vars.theta.consensus_prev, ws.vars.v_mag.consensus_prev,
        ws.vars.theta.rho[], ws.vars.v_mag.rho[],
    )

    ws.r_primal = hypot(ws.r_primal_theta, ws.r_primal_v_mag)
    ws.r_dual = hypot(ws.r_dual_theta, ws.r_dual_v_mag)
    return nothing
end

function compute_residuals!(ws::ACPowerWorkspace)
    ws.r_primal_p, ws.r_primal_q = _primal_residual_branch_ac_components(
        ws.partitions,
        ws.vars.p_flow.consensus,
        ws.vars.q_flow.consensus,
        ws.interzonal,
    )
    ws.r_dual_p, ws.r_dual_q = _dual_residual_branch_paired_components(
        ws.partitions,
        ws.vars.p_flow.consensus, ws.vars.q_flow.consensus,
        ws.vars.p_flow.consensus_prev, ws.vars.q_flow.consensus_prev,
        ws.interzonal, ws.vars.p_flow.rho[], ws.vars.q_flow.rho[],
    )

    ws.r_primal = hypot(ws.r_primal_p, ws.r_primal_q)
    ws.r_dual = hypot(ws.r_dual_p, ws.r_dual_q)
    return nothing
end

function compute_residuals!(ws::ACHybridWorkspace)
    ws.r_primal_theta, ws.r_primal_v_mag = _primal_residual_bus_ac_components(
        ws.partitions,
        ws.vars.theta.consensus,
        ws.vars.v_mag.consensus,
    )
    ws.r_dual_theta, ws.r_dual_v_mag = _dual_residual_bus_ac_components(
        ws.partitions,
        ws.vars.theta.consensus, ws.vars.v_mag.consensus,
        ws.vars.theta.consensus_prev, ws.vars.v_mag.consensus_prev,
        ws.vars.theta.rho[], ws.vars.v_mag.rho[],
    )

    ws.r_primal_p, ws.r_primal_q = _primal_residual_branch_ac_components(
        ws.partitions,
        ws.vars.p_flow.consensus,
        ws.vars.q_flow.consensus,
        ws.interzonal,
    )
    ws.r_dual_p, ws.r_dual_q = _dual_residual_branch_paired_components(
        ws.partitions,
        ws.vars.p_flow.consensus, ws.vars.q_flow.consensus,
        ws.vars.p_flow.consensus_prev, ws.vars.q_flow.consensus_prev,
        ws.interzonal, ws.vars.p_flow.rho[], ws.vars.q_flow.rho[],
    )

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

# ---- SOCWR ----

function compute_residuals!(ws::SOCWRStateWorkspace)
    ws.r_primal_bus = _primal_residual_bus_single(
        ws.partitions,
        ws.vars.v2.consensus,
        snapshot_v2!,
    )
    ws.r_dual_bus = _dual_residual_bus_single(
        ws.partitions,
        ws.vars.v2.consensus,
        ws.vars.v2.consensus_prev,
        ws.vars.v2.rho[],
    )
    ws.r_primal_r, ws.r_primal_t = _primal_residual_branch_paired_components(
        ws.partitions,
        ws.vars.r.consensus,
        ws.vars.t.consensus,
        ws.interzonal,
        _socwr_rt_values,
    )
    ws.r_dual_r, ws.r_dual_t = _dual_residual_branch_paired_components(
        ws.partitions,
        ws.vars.r.consensus, ws.vars.t.consensus,
        ws.vars.r.consensus_prev, ws.vars.t.consensus_prev,
        ws.interzonal, ws.vars.r.rho[], ws.vars.t.rho[],
    )

    ws.r_primal = sqrt(ws.r_primal_bus^2 + ws.r_primal_r^2 + ws.r_primal_t^2)
    ws.r_dual = sqrt(ws.r_dual_bus^2 + ws.r_dual_r^2 + ws.r_dual_t^2)
    return nothing
end

function compute_residuals!(ws::SOCWRPowerWorkspace)
    ws.r_primal_p, ws.r_primal_q = _primal_residual_branch_paired_components(
        ws.partitions,
        ws.vars.p_flow.consensus,
        ws.vars.q_flow.consensus,
        ws.interzonal,
        _socwr_pq_values,
    )
    ws.r_dual_p, ws.r_dual_q = _dual_residual_branch_paired_components(
        ws.partitions,
        ws.vars.p_flow.consensus, ws.vars.q_flow.consensus,
        ws.vars.p_flow.consensus_prev, ws.vars.q_flow.consensus_prev,
        ws.interzonal, ws.vars.p_flow.rho[], ws.vars.q_flow.rho[],
    )

    ws.r_primal = hypot(ws.r_primal_p, ws.r_primal_q)
    ws.r_dual = hypot(ws.r_dual_p, ws.r_dual_q)
    return nothing
end

function compute_residuals!(ws::SOCWRHybridWorkspace)
    ws.r_primal_bus = _primal_residual_bus_single(
        ws.partitions,
        ws.vars.v2.consensus,
        snapshot_v2!,
    )
    ws.r_dual_bus = _dual_residual_bus_single(
        ws.partitions,
        ws.vars.v2.consensus,
        ws.vars.v2.consensus_prev,
        ws.vars.v2.rho[],
    )
    ws.r_primal_p, ws.r_primal_q = _primal_residual_branch_paired_components(
        ws.partitions,
        ws.vars.p_flow.consensus,
        ws.vars.q_flow.consensus,
        ws.interzonal,
        _socwr_pq_values,
    )
    ws.r_dual_p, ws.r_dual_q = _dual_residual_branch_paired_components(
        ws.partitions,
        ws.vars.p_flow.consensus, ws.vars.q_flow.consensus,
        ws.vars.p_flow.consensus_prev, ws.vars.q_flow.consensus_prev,
        ws.interzonal, ws.vars.p_flow.rho[], ws.vars.q_flow.rho[],
    )
    ws.r_primal_r, ws.r_primal_t = _primal_residual_branch_paired_components(
        ws.partitions,
        ws.vars.r.consensus,
        ws.vars.t.consensus,
        ws.interzonal,
        _socwr_rt_values,
    )
    ws.r_dual_r, ws.r_dual_t = _dual_residual_branch_paired_components(
        ws.partitions,
        ws.vars.r.consensus, ws.vars.t.consensus,
        ws.vars.r.consensus_prev, ws.vars.t.consensus_prev,
        ws.interzonal, ws.vars.r.rho[], ws.vars.t.rho[],
    )

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
