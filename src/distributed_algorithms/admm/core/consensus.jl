# ============================================================
# ADMM Consensus — unified for DC · AC · SOCWR
# ============================================================

# ------------------------------------------------------------
# Bus consensus: single variable, linear averaging (DC theta, SOCWR V²)
# ------------------------------------------------------------
function _consensus_bus_single!(
    consensus::Vector{Float64},
    counts::Vector{Int},
    partitions::Partitions,
    snapshot_fn,    # (part) -> Vector{Float64}: shared-bus values
)
    fill!(consensus, 0.0)
    fill!(counts, 0)

    for part in partitions
        vals = snapshot_fn(part)
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            j = part.shared_bus_consensus_idx[i]
            consensus[j] += vals[i]
            counts[j] += 1
        end
    end

    @inbounds for j in eachindex(consensus)
        c = counts[j]
        c > 0 && (consensus[j] /= c)
    end
    return nothing
end

# ------------------------------------------------------------
# Bus consensus: AC (circular theta + linear |V|)
# ------------------------------------------------------------
function _consensus_bus_ac!(
    theta_consensus::Vector{Float64},
    v_mag_consensus::Vector{Float64},
    sin_sum::Vector{Float64},
    cos_sum::Vector{Float64},
    counts::Vector{Int},
    partitions::Partitions{PartitionAC},
)
    fill!(sin_sum, 0.0)
    fill!(cos_sum, 0.0)
    fill!(v_mag_consensus, 0.0)
    fill!(counts, 0)

    for part in partitions
        theta_loc = snapshot_theta_shared!(part)
        v_mag_loc = snapshot_v_mag_shared!(part)
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            j = part.shared_bus_consensus_idx[i]
            s, c = sincos(theta_loc[i])
            sin_sum[j] += s
            cos_sum[j] += c
            v_mag_consensus[j] += v_mag_loc[i]
            counts[j] += 1
        end
    end

    @inbounds for j in eachindex(theta_consensus)
        if counts[j] > 0
            theta_consensus[j] = atan(sin_sum[j], cos_sum[j])
            v_mag_consensus[j] /= counts[j]
        end
    end
    return nothing
end

# ------------------------------------------------------------
# Branch consensus: single variable (DC P)
# ------------------------------------------------------------
function _consensus_branch_single!(
    consensus::Vector{Float64},
    counts::Vector{Int},
    partitions::Partitions,
    interzonal::Vector{Int},
    extract_fn,    # (part, l) -> Float64
)
    fill!(consensus, 0.0)
    fill!(counts, 0)

    for part in partitions
        @inbounds for (k, e) in enumerate(interzonal)
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            consensus[k] += extract_fn(part, l)
            counts[k] += 1
        end
    end

    @inbounds for k in eachindex(consensus)
        c = counts[k]
        c > 0 && (consensus[k] /= c)
    end
    return nothing
end

# ------------------------------------------------------------
# Branch consensus: paired variables (SOCWR PQ, SOCWR RT)
# ------------------------------------------------------------
function _consensus_branch_paired!(
    consensus_a::Vector{Float64},
    consensus_b::Vector{Float64},
    counts::Vector{Int},
    partitions::Partitions,
    interzonal::Vector{Int},
    extract_fn,    # (part, l) -> (Float64, Float64)
)
    fill!(consensus_a, 0.0)
    fill!(consensus_b, 0.0)
    fill!(counts, 0)

    for part in partitions
        @inbounds for (k, e) in enumerate(interzonal)
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            a, b = extract_fn(part, l)
            consensus_a[k] += a
            consensus_b[k] += b
            counts[k] += 1
        end
    end

    @inbounds for k in eachindex(consensus_a)
        c = counts[k]
        if c > 0
            consensus_a[k] /= c
            consensus_b[k] /= c
        end
    end
    return nothing
end

# ------------------------------------------------------------
# Branch consensus: AC P/Q (with isfinite checks for NLP)
# ------------------------------------------------------------
# Uses independent P/Q counters: NLP flows can be NaN, and we don't want
# a bad P to poison Q (or vice versa).
function _consensus_branch_ac!(
    p_consensus::Vector{Float64},
    q_consensus::Vector{Float64},
    p_counts::Vector{Int},
    q_counts::Vector{Int},
    partitions::Partitions{PartitionAC},
    interzonal::Vector{Int},
)
    fill!(p_consensus, 0.0)
    fill!(q_consensus, 0.0)
    fill!(p_counts, 0)
    fill!(q_counts, 0)

    for part in partitions
        @inbounds for (k, e) in enumerate(interzonal)
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            p, q = _pq_flow_from_exprs(part, l)
            p_val = JuMP.value(p)
            q_val = JuMP.value(q)
            if isfinite(p_val)
                p_consensus[k] += p_val
                p_counts[k] += 1
            end
            if isfinite(q_val)
                q_consensus[k] += q_val
                q_counts[k] += 1
            end
        end
    end

    @inbounds for k in eachindex(p_consensus)
        np = p_counts[k]
        nq = q_counts[k]
        np > 0 && (p_consensus[k] /= np)
        nq > 0 && (q_consensus[k] /= nq)
    end
    return nothing
end

# ============================================================
# Dispatch: consensus!(ws)
# ============================================================

# ---- DC ----

consensus!(ws::DCStateWorkspace) =
    _consensus_bus_single!(
        ws.theta_consensus,
        ws.bus_counts,
        ws.partitions,
        snapshot_theta_shared!,
    )

consensus!(ws::DCPowerWorkspace) =
    _consensus_branch_single!(
        ws.p_flow_consensus,
        ws.branch_counts,
        ws.partitions,
        ws.interzonal,
        _dc_flow_value,
    )

function consensus!(ws::DCHybridWorkspace)
    _consensus_bus_single!(
        ws.theta_consensus,
        ws.bus_counts,
        ws.partitions,
        snapshot_theta_shared!,
    )
    _consensus_branch_single!(
        ws.p_flow_consensus,
        ws.branch_counts,
        ws.partitions,
        ws.interzonal,
        _dc_flow_value,
    )
    return nothing
end

# ---- AC ----

consensus!(ws::ACStateWorkspace) =
    _consensus_bus_ac!(
        ws.theta_consensus,
        ws.v_mag_consensus,
        ws.sin_sum,
        ws.cos_sum,
        ws.bus_counts,
        ws.partitions,
    )

consensus!(ws::ACPowerWorkspace) =
    _consensus_branch_ac!(
        ws.p_flow_consensus,
        ws.q_flow_consensus,
        ws.p_branch_counts,
        ws.q_branch_counts,
        ws.partitions,
        ws.interzonal,
    )

function consensus!(ws::ACHybridWorkspace)
    _consensus_bus_ac!(
        ws.theta_consensus,
        ws.v_mag_consensus,
        ws.sin_sum,
        ws.cos_sum,
        ws.bus_counts,
        ws.partitions,
    )
    _consensus_branch_ac!(
        ws.p_flow_consensus,
        ws.q_flow_consensus,
        ws.p_branch_counts,
        ws.q_branch_counts,
        ws.partitions,
        ws.interzonal,
    )
    return nothing
end

# ---- SOCWR ----

function consensus!(ws::SOCWRStateWorkspace)
    _consensus_bus_single!(
        ws.v2_consensus,
        ws.bus_counts,
        ws.partitions,
        snapshot_v2_shared!,
    )
    _consensus_branch_paired!(
        ws.r_consensus,
        ws.t_consensus,
        ws.branch_counts,
        ws.partitions,
        ws.interzonal,
        _socwr_rt_values,
    )
    return nothing
end

consensus!(ws::SOCWRPowerWorkspace) =
    _consensus_branch_paired!(
        ws.p_flow_consensus,
        ws.q_flow_consensus,
        ws.branch_counts,
        ws.partitions,
        ws.interzonal,
        _socwr_pq_values,
    )

function consensus!(ws::SOCWRHybridWorkspace)
    _consensus_bus_single!(
        ws.v2_consensus,
        ws.bus_counts,
        ws.partitions,
        snapshot_v2_shared!,
    )
    _consensus_branch_paired!(
        ws.p_flow_consensus,
        ws.q_flow_consensus,
        ws.branch_counts,
        ws.partitions,
        ws.interzonal,
        _socwr_pq_values,
    )
    _consensus_branch_paired!(
        ws.r_consensus,
        ws.t_consensus,
        ws.branch_counts,
        ws.partitions,
        ws.interzonal,
        _socwr_rt_values,
    )
    return nothing
end
