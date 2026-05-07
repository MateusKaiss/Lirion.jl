# ============================================================
# ALADIN-DC — Consensus pair types and builders
# ============================================================
#
# ALADIN imposes the coupling `Σ A_i x_i = b` row-by-row through
# explicit pair structures (one per λ component). This file defines
# the four pair types used across coupling variants and the builders
# that populate the pair lists once at workspace construction time.
#
# Per-partition runtime setup (NLP build, quad_static, p_inj
# variables) lives in `partition.jl`, mirroring `admm/dc/partition.jl`.

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
    zone_i::Int
    p_inj_col_i::Int
    zone_j::Int
    p_inj_col_j::Int
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
#
# Walks every zone, registers (zone_pos, theta_col) for each shared
# bus it touches, then emits a star spanning tree per consensus_idx
# (root = lowest zone_pos). N zones touching a bus → N-1 pairs, which
# preserves LICQ in the coordinator KKT.

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
            push!(pairs, SharedBusPair(
                consensus_idx,
                bus_id_by_cidx[consensus_idx],
                root_zone, root_col,
                peer_zone, peer_col,
            ))
        end
    end
    return pairs
end

# ------------------------------------------------------------
# Interzonal line pairs
# ------------------------------------------------------------
#
# Each interzonal line touches exactly two zones; emit one pair per
# line carrying its global index and the local index in each zone.

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
        push!(pairs, InterzonalLinePair(
            consensus_idx, line_global,
            zone_i, local_i,
            zone_j, local_j,
        ))
    end
    return pairs
end

# ------------------------------------------------------------
# Dummy line pairs
# ------------------------------------------------------------
#
# DummyBus sectioning replaces a shared bus by two phantom lines (one
# per zone, ending at a per-zone dummy bus). For each
# (line_global_a, zone_a, line_global_b, zone_b) tuple from the
# system, emit a pair carrying the θ columns of both endpoints AND
# the susceptances — the consensus row is
# `b_a·(θ_from_a − θ_to_a) = b_b·(θ_from_b − θ_to_b)`.

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
        to_col_i   = zone_data_a.var_to_col[theta_var_a[to_a, 1]]
        from_col_j = zone_data_b.var_to_col[theta_var_b[from_b, 1]]
        to_col_j   = zone_data_b.var_to_col[theta_var_b[to_b, 1]]

        b_i = part_a.branch.b[local_a]
        b_j = part_b.branch.b[local_b]

        push!(pairs, DummyLinePairDC(
            consensus_idx,
            zone_a, from_col_i, to_col_i, b_i,
            zone_b, from_col_j, to_col_j, b_j,
        ))
    end
    return pairs
end

# ------------------------------------------------------------
# Split bus pairs
# ------------------------------------------------------------
#
# SplitBus sectioning duplicates each shared bus with auxiliary
# `p_inj` variables (one per zone). Emit a pair carrying the column
# of each replica so the coordinator can enforce
# `p_inj_a + p_inj_b = 0`.

function _build_split_bus_pairs(
    partitions::Partitions{PartitionDC},
    zones::Vector{AladinZoneData},
    sys_bus_pairs::Vector{NTuple{3, Int}},
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}},
)::Vector{SplitBusPair}
    pairs = SplitBusPair[]
    isempty(sys_bus_pairs) && return pairs
    sizehint!(pairs, length(sys_bus_pairs))

    for (consensus_idx, (bus_id, zone_a, zone_b)) in enumerate(sys_bus_pairs)
        p_inj_var_a = get(p_inj_vars[zone_a], bus_id, nothing)
        p_inj_var_b = get(p_inj_vars[zone_b], bus_id, nothing)
        (p_inj_var_a === nothing || p_inj_var_b === nothing) && continue
        p_inj_col_a = zones[zone_a].var_to_col[p_inj_var_a]
        p_inj_col_b = zones[zone_b].var_to_col[p_inj_var_b]
        push!(pairs, SplitBusPair(
            consensus_idx, bus_id,
            zone_a, p_inj_col_a,
            zone_b, p_inj_col_b,
        ))
    end
    return pairs
end
