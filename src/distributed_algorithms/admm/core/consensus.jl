# ============================================================
# Helpers
# ============================================================

# ------------------------------------------------------------
# Bus
# ------------------------------------------------------------

function _consensus_bus_single!(
    consensus::Vector{Float64},
    counts::Vector{Int},
    partitions::Partitions,
    snapshot_fn,
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
# Branch
# ------------------------------------------------------------

function _consensus_branch_single!(
    consensus::AbstractVector{Float64},
    counts::AbstractVector{Int},
    partitions::Partitions,
    interzonal::Vector{Int},
    extract_fn,
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

function _consensus_branch_paired!(
    consensus_a::AbstractVector{Float64},
    consensus_b::AbstractVector{Float64},
    counts::AbstractVector{Int},
    partitions::Partitions,
    interzonal::Vector{Int},
    extract_fn,
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

function _consensus_branch_ac!(
    p_consensus::AbstractVector{Float64},
    q_consensus::AbstractVector{Float64},
    p_counts::AbstractVector{Int},
    q_counts::AbstractVector{Int},
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

# ------------------------------------------------------------
# Branch pairs
# ------------------------------------------------------------

function _consensus_branch_pairs!(
    consensus::AbstractVector{Float64},
    counts::AbstractVector{Int},
    partitions::Partitions,
    pairs::Vector{NTuple{4, Int}},
    extract_fn,
)
    isempty(pairs) && return nothing
    fill!(consensus, 0.0)
    fill!(counts, 0)
    @inbounds for (k, (ga, za, gb, zb)) in enumerate(pairs)
        pa = partitions.parts[za]
        la = pa.branch_idx_full_to_zone[ga]
        pb = partitions.parts[zb]
        lb = pb.branch_idx_full_to_zone[gb]
        if la != 0
            consensus[k] += extract_fn(pa, la)
            counts[k] += 1
        end
        if lb != 0
            consensus[k] += extract_fn(pb, lb)
            counts[k] += 1
        end
        counts[k] > 0 && (consensus[k] /= counts[k])
    end
    return nothing
end

function _consensus_branch_pairs_ac!(
    p_consensus::AbstractVector{Float64},
    q_consensus::AbstractVector{Float64},
    counts::AbstractVector{Int},
    partitions::Partitions{PartitionAC},
    pairs::Vector{NTuple{4, Int}},
)
    isempty(pairs) && return nothing
    @inbounds for (k, (ga, za, gb, zb)) in enumerate(pairs)
        pa = partitions.parts[za]
        pb = partitions.parts[zb]
        la = pa.branch_idx_full_to_zone[ga]
        lb = pb.branch_idx_full_to_zone[gb]

        p_sum = 0.0
        q_sum = 0.0
        p_n = 0
        q_n = 0
        if la != 0
            p_expr, q_expr = _pq_flow_from_exprs(pa, la)
            p_val = JuMP.value(p_expr)
            q_val = JuMP.value(q_expr)
            if isfinite(p_val)
                p_sum += p_val
                p_n += 1
            end
            if isfinite(q_val)
                q_sum += q_val
                q_n += 1
            end
        end
        if lb != 0
            p_expr, q_expr = _pq_flow_from_exprs(pb, lb)
            p_val = JuMP.value(p_expr)
            q_val = JuMP.value(q_expr)
            if isfinite(p_val)
                p_sum += p_val
                p_n += 1
            end
            if isfinite(q_val)
                q_sum += q_val
                q_n += 1
            end
        end
        p_consensus[k] = p_n > 0 ? p_sum / p_n : 0.0
        q_consensus[k] = q_n > 0 ? q_sum / q_n : 0.0
        counts[k] = 2
    end
    return nothing
end

# ------------------------------------------------------------
# Bus pairs
# ------------------------------------------------------------

function _consensus_bus_pairs_signed!(
    consensus::AbstractVector{Float64},
    partitions::Partitions,
    bus_pairs::Vector{NTuple{3, Int}},
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}},
)
    isempty(bus_pairs) && return nothing
    fill!(consensus, 0.0)
    @inbounds for (k, (b_id, za, zb)) in enumerate(bus_pairs)
        va = _p_inj_value(partitions.parts[za], p_inj_vars[za], b_id)
        vb = _p_inj_value(partitions.parts[zb], p_inj_vars[zb], b_id)
        consensus[k] = 0.5 * (va - vb)
    end
    return nothing
end

# ============================================================
# DC
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

consensus!(ws::DCStateWorkspace) =
    _consensus_bus_single!(
        ws.vars.theta.consensus,
        ws.bus_counts,
        ws.partitions,
        snapshot_theta_shared!,
    )

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function consensus!(ws::DCPowerWorkspace{TieLine})
    _consensus_branch_single!(
        ws.vars.p_flow.consensus,
        ws.branch_counts,
        ws.partitions,
        ws.interzonal,
        _dc_flow_value,
    )
    return nothing
end

function consensus!(ws::DCPowerWorkspace{DummyBus})
    _consensus_branch_pairs!(
        ws.vars.p_dummy.consensus,
        ws.pair_counts,
        ws.partitions,
        ws.branch_pairs,
        _dc_flow_value,
    )
    return nothing
end

# ------------------------------------------------------------
# Hybrid — TieLine
# ------------------------------------------------------------

function consensus!(ws::DCHybridWorkspace{TieLine})
    _consensus_bus_single!(
        ws.vars.theta.consensus,
        ws.bus_counts,
        ws.partitions,
        snapshot_theta_shared!,
    )
    _consensus_branch_single!(
        ws.vars.p_flow.consensus,
        ws.branch_counts,
        ws.partitions,
        ws.interzonal,
        _dc_flow_value,
    )
    return nothing
end

# ------------------------------------------------------------
# Hybrid — DummyBus
# ------------------------------------------------------------

function consensus!(ws::DCHybridWorkspace{DummyBus})
    _consensus_bus_single!(
        ws.vars.theta.consensus,
        ws.bus_counts,
        ws.partitions,
        snapshot_theta_shared!,
    )
    _consensus_branch_pairs!(
        ws.vars.p_dummy.consensus,
        ws.pair_counts,
        ws.partitions,
        ws.branch_pairs,
        _dc_flow_value,
    )
    return nothing
end

# ------------------------------------------------------------
# Hybrid — SplitBus
# ------------------------------------------------------------

function consensus!(ws::DCHybridWorkspace{SplitBus})
    _consensus_bus_single!(
        ws.vars.theta.consensus,
        ws.bus_counts,
        ws.partitions,
        snapshot_theta_shared!,
    )
    _consensus_inj_groups!(
        ws.vars.p_inj.consensus,
        ws.partitions,
        ws.inj_groups,
        ws.inj_copy_zone,
        ws.p_inj_vars,
    )
    return nothing
end

# ============================================================
# AC
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

consensus!(ws::ACStateWorkspace) =
    _consensus_bus_ac!(
        ws.vars.theta.consensus,
        ws.vars.v_mag.consensus,
        ws.sin_sum,
        ws.cos_sum,
        ws.bus_counts,
        ws.partitions,
    )

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function consensus!(ws::ACPowerWorkspace{TieLine})
    _consensus_branch_ac!(
        ws.vars.p_flow.consensus,
        ws.vars.q_flow.consensus,
        ws.p_branch_counts,
        ws.q_branch_counts,
        ws.partitions,
        ws.interzonal,
    )
    return nothing
end

function consensus!(ws::ACPowerWorkspace{DummyBus})
    _consensus_branch_pairs_ac!(
        ws.vars.p_dummy.consensus,
        ws.vars.q_dummy.consensus,
        ws.pair_counts,
        ws.partitions,
        ws.branch_pairs,
    )
    return nothing
end

# ------------------------------------------------------------
# Hybrid — TieLine
# ------------------------------------------------------------

function consensus!(ws::ACHybridWorkspace{TieLine})
    _consensus_bus_ac!(
        ws.vars.theta.consensus,
        ws.vars.v_mag.consensus,
        ws.sin_sum,
        ws.cos_sum,
        ws.bus_counts,
        ws.partitions,
    )
    _consensus_branch_ac!(
        ws.vars.p_flow.consensus,
        ws.vars.q_flow.consensus,
        ws.p_branch_counts,
        ws.q_branch_counts,
        ws.partitions,
        ws.interzonal,
    )
    return nothing
end

# ------------------------------------------------------------
# Hybrid — DummyBus
# ------------------------------------------------------------

function consensus!(ws::ACHybridWorkspace{DummyBus})
    _consensus_bus_ac!(
        ws.vars.theta.consensus,
        ws.vars.v_mag.consensus,
        ws.sin_sum,
        ws.cos_sum,
        ws.bus_counts,
        ws.partitions,
    )
    _consensus_branch_pairs_ac!(
        ws.vars.p_dummy.consensus,
        ws.vars.q_dummy.consensus,
        ws.pair_counts,
        ws.partitions,
        ws.branch_pairs,
    )
    return nothing
end

# ------------------------------------------------------------
# Hybrid — SplitBus
# ------------------------------------------------------------

function consensus!(ws::ACHybridWorkspace{SplitBus})
    _consensus_bus_ac!(
        ws.vars.theta.consensus,
        ws.vars.v_mag.consensus,
        ws.sin_sum,
        ws.cos_sum,
        ws.bus_counts,
        ws.partitions,
    )
    _consensus_inj_groups!(
        ws.vars.p_inj.consensus, ws.partitions,
        ws.inj_groups, ws.inj_copy_zone, ws.p_inj_vars,
    )
    _consensus_inj_groups!(
        ws.vars.q_inj.consensus, ws.partitions,
        ws.inj_groups, ws.inj_copy_zone, ws.q_inj_vars,
    )
    return nothing
end

# ============================================================
# SOCWR
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function consensus!(ws::SOCWRStateWorkspace)
    _consensus_bus_single!(
        ws.vars.v2.consensus,
        ws.bus_counts,
        ws.partitions,
        snapshot_v2_shared!,
    )
    _consensus_branch_paired!(
        ws.vars.r.consensus,
        ws.vars.t.consensus,
        ws.branch_counts,
        ws.partitions,
        ws.interzonal,
        _socwr_rt_values,
    )
    return nothing
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function consensus!(ws::SOCWRPowerWorkspace{TieLine})
    _consensus_branch_paired!(
        ws.vars.p_flow.consensus,
        ws.vars.q_flow.consensus,
        ws.branch_counts,
        ws.partitions,
        ws.interzonal,
        _socwr_pq_values,
    )
    return nothing
end

function consensus!(ws::SOCWRPowerWorkspace{DummyBus})
    _consensus_branch_pairs!(
        ws.vars.p_dummy.consensus,
        ws.pair_counts,
        ws.partitions,
        ws.branch_pairs,
        _socwr_p_value,
    )
    _consensus_branch_pairs!(
        ws.vars.q_dummy.consensus,
        ws.pair_counts,
        ws.partitions,
        ws.branch_pairs,
        _socwr_q_value,
    )
    return nothing
end

# ------------------------------------------------------------
# Hybrid — TieLine
# ------------------------------------------------------------

function consensus!(ws::SOCWRHybridWorkspace{TieLine})
    _consensus_bus_single!(
        ws.vars.v2.consensus,
        ws.bus_counts,
        ws.partitions,
        snapshot_v2_shared!,
    )
    _consensus_branch_paired!(
        ws.vars.p_flow.consensus,
        ws.vars.q_flow.consensus,
        ws.branch_counts,
        ws.partitions,
        ws.interzonal,
        _socwr_pq_values,
    )
    _consensus_branch_paired!(
        ws.vars.r.consensus,
        ws.vars.t.consensus,
        ws.branch_counts,
        ws.partitions,
        ws.interzonal,
        _socwr_rt_values,
    )
    return nothing
end

# ------------------------------------------------------------
# Hybrid — DummyBus
# ------------------------------------------------------------

function consensus!(ws::SOCWRHybridWorkspace{DummyBus})
    _consensus_bus_single!(
        ws.vars.v2.consensus,
        ws.bus_counts,
        ws.partitions,
        snapshot_v2_shared!,
    )
    _consensus_branch_pairs!(
        ws.vars.p_dummy.consensus,
        ws.pair_counts,
        ws.partitions,
        ws.branch_pairs,
        _socwr_p_value,
    )
    _consensus_branch_pairs!(
        ws.vars.q_dummy.consensus,
        ws.pair_counts,
        ws.partitions,
        ws.branch_pairs,
        _socwr_q_value,
    )
    return nothing
end

# ------------------------------------------------------------
# Hybrid — SplitBus
# ------------------------------------------------------------

function consensus!(ws::SOCWRHybridWorkspace{SplitBus})
    _consensus_bus_single!(
        ws.vars.v2.consensus,
        ws.bus_counts,
        ws.partitions,
        snapshot_v2_shared!,
    )
    _consensus_inj_groups!(
        ws.vars.p_inj.consensus, ws.partitions,
        ws.inj_groups, ws.inj_copy_zone, ws.p_inj_vars,
    )
    _consensus_inj_groups!(
        ws.vars.q_inj.consensus, ws.partitions,
        ws.inj_groups, ws.inj_copy_zone, ws.q_inj_vars,
    )
    return nothing
end
