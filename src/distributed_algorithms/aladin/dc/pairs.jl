# ============================================================
# ALADIN-DC — Consensus pair types and builders
# ============================================================

# ============================================================
# Pair types (consensus structures)
# ============================================================

# ------------------------------------------------------------
# Shared bus pair
# ------------------------------------------------------------

struct SharedBusPair
    consensus_idx::Int
    bus_id::Int
    zone_i::Int
    col_i::Int
    zone_j::Int
    col_j::Int
end

# ------------------------------------------------------------
# Interzonal line pair (TieLine sectioning)
# ------------------------------------------------------------

struct InterzonalLinePair
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

struct SplitBusPair
    consensus_idx::Int
    bus_id::Int
    zones::Vector{Int}
    p_inj_cols::Vector{Int}
end

# ------------------------------------------------------------
# Dummy line pair (DummyBus sectioning)
# ------------------------------------------------------------

struct DummyLinePairDC
    consensus_idx::Int
    zone_i::Int
    from_col_i::Int
    to_col_i::Int
    b_i::Float64
    zone_j::Int
    from_col_j::Int
    to_col_j::Int
    b_j::Float64
end

# ============================================================
# Pair builders
# ============================================================

# ------------------------------------------------------------
# Shared bus pairs
# ------------------------------------------------------------

function _build_shared_bus_pairs(
    partitions::Partitions{PartitionDC},
    zones::Vector{AladinZoneData},
)::Vector{SharedBusPair}
    visits = Dict{Int, Vector{Tuple{Int, Int}}}()    # consensus_idx → [(zone_pos, theta_col)]
    bus_id_by_cidx = Dict{Int, Int}()                # consensus_idx → physical bus id

    for (zone_pos, part) in enumerate(partitions.parts)
        zone_data = zones[zone_pos]
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            consensus_idx = part.shared_bus_consensus_idx[i]
            consensus_idx == 0 && continue
            theta_col = zone_data.shared_theta_cols[i]
            cidx_visits = get!(visits, consensus_idx, Tuple{Int, Int}[])
            push!(cidx_visits, (zone_pos, theta_col))
            bus_id_by_cidx[consensus_idx] = part.zone_shared_bus_ids[i]
        end
    end

    pairs = SharedBusPair[]
    for consensus_idx in sort!(collect(keys(visits)))
        cidx_visits = visits[consensus_idx]
        length(cidx_visits) >= 2 || continue
        sort!(cidx_visits; by = first)
        root_zone, root_col = cidx_visits[1]
        for k in 2:length(cidx_visits)
            peer_zone, peer_col = cidx_visits[k]
            push!(
                pairs,
                SharedBusPair(
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
    partitions::Partitions{PartitionDC},
    interzonal::Vector{Int},
)::Vector{InterzonalLinePair}
    pairs = InterzonalLinePair[]

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
            InterzonalLinePair(
                consensus_idx, line_global,
                zone_i, local_i,
                zone_j, local_j,
            ),
        )
    end
    return pairs
end

# ------------------------------------------------------------
# Dummy line pairs
# ------------------------------------------------------------

function _build_dummy_line_pairs(
    partitions::Partitions{PartitionDC},
    zones::Vector{AladinZoneData},
    sys_branch_pairs::Vector{NTuple{4, Int}},
)::Vector{DummyLinePairDC}
    pairs = DummyLinePairDC[]
    isempty(sys_branch_pairs) && return pairs
    sizehint!(pairs, length(sys_branch_pairs))

    for (consensus_idx, (line_global_a, zone_a, line_global_b, zone_b)) in
        enumerate(sys_branch_pairs)
        part_a = partitions.parts[zone_a]
        part_b = partitions.parts[zone_b]
        local_a = part_a.branch_idx_full_to_zone[line_global_a]
        local_b = part_b.branch_idx_full_to_zone[line_global_b]
        (local_a == 0 || local_b == 0) && continue

        zone_data_a = zones[zone_a]
        zone_data_b = zones[zone_b]

        from_a = part_a.branch.from[local_a]
        to_a = part_a.branch.to[local_a]
        from_b = part_b.branch.from[local_b]
        to_b = part_b.branch.to[local_b]

        theta_var_a = part_a.model[:theta]
        theta_var_b = part_b.model[:theta]
        from_col_i = zone_data_a.var_to_col[theta_var_a[from_a, 1]]
        to_col_i = zone_data_a.var_to_col[theta_var_a[to_a, 1]]
        from_col_j = zone_data_b.var_to_col[theta_var_b[from_b, 1]]
        to_col_j = zone_data_b.var_to_col[theta_var_b[to_b, 1]]

        b_i = part_a.branch.b[local_a]
        b_j = part_b.branch.b[local_b]

        push!(
            pairs,
            DummyLinePairDC(
                consensus_idx,
                zone_a, from_col_i, to_col_i, b_i,
                zone_b, from_col_j, to_col_j, b_j,
            ),
        )
    end
    return pairs
end

# ------------------------------------------------------------
# Split bus pairs
# ------------------------------------------------------------

function _build_split_bus_pairs(
    partitions::Partitions{PartitionDC},
    zones::Vector{AladinZoneData},
    bus_groups::Vector{Tuple{Int, Vector{Int}}},
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}},
)::Vector{SplitBusPair}
    pairs = SplitBusPair[]
    isempty(bus_groups) && return pairs
    sizehint!(pairs, length(bus_groups))

    for (consensus_idx, (bus_id, group_zones)) in enumerate(bus_groups)
        copy_zones = Int[]
        copy_cols = Int[]
        ok = true
        for zone in group_zones
            p_inj_var = get(p_inj_vars[zone], bus_id, nothing)
            if p_inj_var === nothing
                ok = false
                break
            end
            push!(copy_zones, zone)
            push!(copy_cols, zones[zone].var_to_col[p_inj_var])
        end
        (ok && length(copy_zones) >= 2) || continue
        push!(pairs, SplitBusPair(consensus_idx, bus_id, copy_zones, copy_cols))
    end
    return pairs
end
