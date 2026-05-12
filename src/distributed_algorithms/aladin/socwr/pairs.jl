# ============================================================
# ALADIN-SOCWR — Consensus pair types and builders
# ============================================================

# ============================================================
# Pair types (consensus structures)
# ============================================================

# ------------------------------------------------------------
# Shared bus pair
# ------------------------------------------------------------

struct SharedBusPairSOCWR
    consensus_idx::Int
    bus_id::Int
    zone_i::Int
    v2_col_i::Int
    zone_j::Int
    v2_col_j::Int
end

# ------------------------------------------------------------
# Interzonal line pair (TieLine sectioning)
# ------------------------------------------------------------

struct InterzonalLinePairSOCWR
    consensus_idx::Int
    line_global::Int
    zone_i::Int
    line_local_i::Int
    zone_j::Int
    line_local_j::Int
end

# ------------------------------------------------------------
# Split bus pair (SplitBus sectioning)
# ------------------------------------------------------------

struct SplitBusPairSOCWR
    consensus_idx::Int
    bus_id::Int
    zone_i::Int
    p_inj_col_i::Int
    q_inj_col_i::Int
    zone_j::Int
    p_inj_col_j::Int
    q_inj_col_j::Int
end

# ============================================================
# Pair builders
# ============================================================

# ------------------------------------------------------------
# Shared bus pairs
# ------------------------------------------------------------

function _build_shared_bus_pairs(
    partitions::Partitions{PartitionSOCWR},
    zones::Vector{AladinZoneDataSOCWR},
)::Vector{SharedBusPairSOCWR}
    visits = Dict{Int, Vector{Tuple{Int, Int}}}()    # cidx → [(zone_pos, v2_col)]
    bus_id_by_cidx = Dict{Int, Int}()                # cidx → physical bus id

    for (zone_pos, part) in enumerate(partitions.parts)
        zone_data = zones[zone_pos]
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            consensus_idx = part.shared_bus_consensus_idx[i]
            consensus_idx == 0 && continue
            v2_col = zone_data.shared_v2_cols[i]
            cidx_visits = get!(visits, consensus_idx, Tuple{Int, Int}[])
            push!(cidx_visits, (zone_pos, v2_col))
            bus_id_by_cidx[consensus_idx] = part.zone_shared_bus_ids[i]
        end
    end

    pairs = SharedBusPairSOCWR[]
    for consensus_idx in sort!(collect(keys(visits)))
        cidx_visits = visits[consensus_idx]
        length(cidx_visits) >= 2 || continue
        sort!(cidx_visits; by = first)
        root_zone, root_col = cidx_visits[1]
        for k in 2:length(cidx_visits)
            peer_zone, peer_col = cidx_visits[k]
            push!(
                pairs,
                SharedBusPairSOCWR(
                    consensus_idx,
                    bus_id_by_cidx[consensus_idx],
                    root_zone, root_col,
                    peer_zone, peer_col,
                ),
            )
        end
    end
    return pairs
end

# ------------------------------------------------------------
# Interzonal line pairs
# ------------------------------------------------------------

function _build_interzonal_line_pairs(
    partitions::Partitions{PartitionSOCWR},
    interzonal::Vector{Int},
)::Vector{InterzonalLinePairSOCWR}
    pairs = InterzonalLinePairSOCWR[]
    for (consensus_idx, line_global) in enumerate(interzonal)
        line_visits = Tuple{Int, Int}[]
        for (zone_pos, part) in enumerate(partitions.parts)
            part.branch_idx_full_to_zone[line_global] == 0 && continue
            local_idx = findfirst(==(line_global), part.zone_edges)
            local_idx === nothing && continue
            push!(line_visits, (zone_pos, local_idx))
        end
        length(line_visits) == 2 || continue
        zone_i, local_i = line_visits[1]
        zone_j, local_j = line_visits[2]
        push!(
            pairs,
            InterzonalLinePairSOCWR(
                consensus_idx, line_global,
                zone_i, local_i,
                zone_j, local_j,
            ),
        )
    end
    return pairs
end

# ------------------------------------------------------------
# Split bus pairs
# ------------------------------------------------------------

function _build_split_bus_pairs(
    zones::Vector{AladinZoneDataSOCWR},
    sys_bus_pairs::Vector{NTuple{3, Int}},
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}},
    q_inj_vars::Vector{Dict{Int, JuMP.VariableRef}},
)::Vector{SplitBusPairSOCWR}
    pairs = SplitBusPairSOCWR[]
    isempty(sys_bus_pairs) && return pairs
    sizehint!(pairs, length(sys_bus_pairs))

    for (consensus_idx, (bus_id, zone_a, zone_b)) in enumerate(sys_bus_pairs)
        p_inj_var_a = get(p_inj_vars[zone_a], bus_id, nothing)
        q_inj_var_a = get(q_inj_vars[zone_a], bus_id, nothing)
        p_inj_var_b = get(p_inj_vars[zone_b], bus_id, nothing)
        q_inj_var_b = get(q_inj_vars[zone_b], bus_id, nothing)
        (
            p_inj_var_a === nothing || q_inj_var_a === nothing ||
            p_inj_var_b === nothing || q_inj_var_b === nothing
        ) && continue
        p_inj_col_a = zones[zone_a].var_to_col[p_inj_var_a]
        q_inj_col_a = zones[zone_a].var_to_col[q_inj_var_a]
        p_inj_col_b = zones[zone_b].var_to_col[p_inj_var_b]
        q_inj_col_b = zones[zone_b].var_to_col[q_inj_var_b]
        push!(
            pairs,
            SplitBusPairSOCWR(
                consensus_idx, bus_id,
                zone_a, p_inj_col_a, q_inj_col_a,
                zone_b, p_inj_col_b, q_inj_col_b,
            ),
        )
    end
    return pairs
end
