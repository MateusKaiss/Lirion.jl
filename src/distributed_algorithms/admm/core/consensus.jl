# ============================================================
# ADMM Consensus — unified for DC · AC · SOCP
# ============================================================

# ------------------------------------------------------------
# Bus consensus: single variable, linear averaging (DC theta, SOCP V²)
# ------------------------------------------------------------
function _consensus_bus_single!(
    consensus::Vector{Float64},
    partitions::Partitions,
    consensus_bus_idx_map::Dict{Int, Int},
    snapshot_fn,    # (part) -> Vector{Float64}: shared bus values
)
    fill!(consensus, 0.0)
    counts = zeros(Int, length(consensus))

    for part in partitions
        vals = snapshot_fn(part)
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            j = consensus_bus_idx_map[part.zone_shared_bus_ids[i]]
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
    partitions::Partitions{PartitionAC},
    consensus_bus_idx_map::Dict{Int, Int},
)
    fill!(theta_consensus, 0.0)
    fill!(v_mag_consensus, 0.0)
    n = length(theta_consensus)

    sin_sum = zeros(Float64, n)
    cos_sum = zeros(Float64, n)
    counts = zeros(Int, n)

    for part in partitions
        theta_loc = snapshot_theta_shared!(part)
        v_mag_loc = snapshot_v_mag_shared!(part)
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            j = consensus_bus_idx_map[part.zone_shared_bus_ids[i]]
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
    consensus::Dict{Int, Float64},
    partitions::Partitions,
    interzonal::Vector{Int},
    extract_fn,    # (part, l) -> Float64
)
    sums = Dict{Int, Float64}()
    counts = Dict{Int, Int}()

    for part in partitions
        @inbounds for e in part.zone_edges
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            sums[e] = get(sums, e, 0.0) + extract_fn(part, l)
            counts[e] = get(counts, e, 0) + 1
        end
    end

    @inbounds for e in interzonal
        n = get(counts, e, 0)
        n > 0 && (consensus[e] = sums[e] / n)
    end
    return nothing
end

# ------------------------------------------------------------
# Branch consensus: paired variables (SOCP PQ, SOCP RT)
# ------------------------------------------------------------
function _consensus_branch_paired!(
    consensus_a::Dict{Int, Float64},
    consensus_b::Dict{Int, Float64},
    partitions::Partitions,
    interzonal::Vector{Int},
    extract_fn,    # (part, l) -> (Float64, Float64)
)
    sum_a = Dict{Int, Float64}()
    sum_b = Dict{Int, Float64}()
    counts = Dict{Int, Int}()

    for part in partitions
        @inbounds for e in part.zone_edges
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            a, b = extract_fn(part, l)
            sum_a[e] = get(sum_a, e, 0.0) + a
            sum_b[e] = get(sum_b, e, 0.0) + b
            counts[e] = get(counts, e, 0) + 1
        end
    end

    @inbounds for e in interzonal
        n = get(counts, e, 0)
        if n > 0
            consensus_a[e] = sum_a[e] / n
            consensus_b[e] = sum_b[e] / n
        end
    end
    return nothing
end

# ------------------------------------------------------------
# Branch consensus: AC P/Q (with isfinite checks for NLP)
# ------------------------------------------------------------
# Specialized variant (not using _consensus_branch_paired!) because AC flow
# expressions can return NaN (NLP convergence quirks), requiring independent
# P/Q counters so that a bad P value doesn't discard a good Q value.
function _consensus_branch_ac!(
    p_consensus::Dict{Int, Float64},
    q_consensus::Dict{Int, Float64},
    partitions::Partitions{PartitionAC},
    interzonal::Vector{Int},
)
    sum_p = Dict{Int, Float64}()
    sum_q = Dict{Int, Float64}()
    count_p = Dict{Int, Int}()
    count_q = Dict{Int, Int}()

    for part in partitions
        @inbounds for e in part.zone_edges
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue
            p, q = _pq_flow_from_exprs(part, l)
            p_val = JuMP.value(p)
            q_val = JuMP.value(q)
            if isfinite(p_val)
                sum_p[e] = get(sum_p, e, 0.0) + p_val
                count_p[e] = get(count_p, e, 0) + 1
            end
            if isfinite(q_val)
                sum_q[e] = get(sum_q, e, 0.0) + q_val
                count_q[e] = get(count_q, e, 0) + 1
            end
        end
    end

    @inbounds for e in interzonal
        n_p = get(count_p, e, 0)
        n_q = get(count_q, e, 0)
        n_p > 0 && (p_consensus[e] = sum_p[e] / n_p)
        n_q > 0 && (q_consensus[e] = sum_q[e] / n_q)
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
        ws.partitions,
        ws.consensus_bus_idx_map,
        snapshot_theta_shared!,
    )

consensus!(ws::DCPowerWorkspace) =
    _consensus_branch_single!(
        ws.p_flow_consensus,
        ws.partitions,
        ws.interzonal,
        _dc_flow_value,
    )

function consensus!(ws::DCHybridWorkspace)
    _consensus_bus_single!(
        ws.theta_consensus,
        ws.partitions,
        ws.consensus_bus_idx_map,
        snapshot_theta_shared!,
    )
    _consensus_branch_single!(
        ws.p_flow_consensus,
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
        ws.partitions,
        ws.consensus_bus_idx_map,
    )

consensus!(ws::ACPowerWorkspace) =
    _consensus_branch_ac!(
        ws.p_flow_consensus,
        ws.q_flow_consensus,
        ws.partitions,
        ws.interzonal,
    )

function consensus!(ws::ACHybridWorkspace)
    _consensus_bus_ac!(
        ws.theta_consensus,
        ws.v_mag_consensus,
        ws.partitions,
        ws.consensus_bus_idx_map,
    )
    _consensus_branch_ac!(
        ws.p_flow_consensus,
        ws.q_flow_consensus,
        ws.partitions,
        ws.interzonal,
    )
    return nothing
end

# ---- SOCP ----

function consensus!(ws::SOCPStateWorkspace)
    _consensus_bus_single!(
        ws.v2_consensus,
        ws.partitions,
        ws.consensus_bus_idx_map,
        snapshot_v2_shared!,
    )
    _consensus_branch_paired!(
        ws.r_consensus,
        ws.t_consensus,
        ws.partitions,
        ws.interzonal,
        _socp_rt_values,
    )
    return nothing
end

consensus!(ws::SOCPPowerWorkspace) =
    _consensus_branch_paired!(
        ws.p_flow_consensus,
        ws.q_flow_consensus,
        ws.partitions,
        ws.interzonal,
        _socp_pq_values,
    )

function consensus!(ws::SOCPHybridWorkspace)
    _consensus_bus_single!(
        ws.v2_consensus,
        ws.partitions,
        ws.consensus_bus_idx_map,
        snapshot_v2_shared!,
    )
    _consensus_branch_paired!(
        ws.p_flow_consensus,
        ws.q_flow_consensus,
        ws.partitions,
        ws.interzonal,
        _socp_pq_values,
    )
    _consensus_branch_paired!(
        ws.r_consensus,
        ws.t_consensus,
        ws.partitions,
        ws.interzonal,
        _socp_rt_values,
    )
    return nothing
end
