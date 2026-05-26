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

# Net-injection consensus over all copies of a split bus: Sigma_z p_z = 0 and
# Sigma_z q_z = 0. `zones`/`p_inj_cols`/`q_inj_cols` are parallel arrays over
# the copies (>= 2), so a bus shared by N zones is one pair carrying N columns.
struct SplitBusPairSOCWR
    consensus_idx::Int
    bus_id::Int
    zones::Vector{Int}
    p_inj_cols::Vector{Int}
    q_inj_cols::Vector{Int}
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
    bus_groups::Vector{Tuple{Int, Vector{Int}}},
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}},
    q_inj_vars::Vector{Dict{Int, JuMP.VariableRef}},
)::Vector{SplitBusPairSOCWR}
    pairs = SplitBusPairSOCWR[]
    isempty(bus_groups) && return pairs
    sizehint!(pairs, length(bus_groups))

    for (consensus_idx, (bus_id, group_zones)) in enumerate(bus_groups)
        copy_zones = Int[]
        p_cols = Int[]
        q_cols = Int[]
        ok = true
        for zone in group_zones
            p_inj_var = get(p_inj_vars[zone], bus_id, nothing)
            q_inj_var = get(q_inj_vars[zone], bus_id, nothing)
            if p_inj_var === nothing || q_inj_var === nothing
                ok = false
                break
            end
            push!(copy_zones, zone)
            push!(p_cols, zones[zone].var_to_col[p_inj_var])
            push!(q_cols, zones[zone].var_to_col[q_inj_var])
        end
        (ok && length(copy_zones) >= 2) || continue
        push!(pairs, SplitBusPairSOCWR(consensus_idx, bus_id, copy_zones, p_cols, q_cols))
    end
    return pairs
end
