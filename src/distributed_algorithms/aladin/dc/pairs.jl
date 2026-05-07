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
# Shared edge pair
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
# Split bus pair
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
# Dummy branch pair
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
    acc = Dict{Int, Vector{Tuple{Int, Int}}}()
    bus_of = Dict{Int, Int}()

    for (zpos, part) in enumerate(partitions.parts)
        zdata = zones[zpos]
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            cidx = part.shared_bus_consensus_idx[i]
            cidx == 0 && continue
            col = zdata.shared_theta_cols[i]
            v = get!(acc, cidx, Tuple{Int, Int}[])
            push!(v, (zpos, col))
            bus_of[cidx] = part.zone_shared_bus_ids[i]
        end
    end

    pairs = SharedBusPair[]
    for cidx in sort!(collect(keys(acc)))
        entries = acc[cidx]
        length(entries) >= 2 || continue
        sort!(entries; by = first)
        zr, cr = entries[1]
        for k in 2:length(entries)
            zk, ck = entries[k]
            push!(pairs, SharedBusPair(cidx, bus_of[cidx], zr, cr, zk, ck))
        end
    end
    return pairs
end

# ------------------------------------------------------------
# Shared edge pairs
# ------------------------------------------------------------

function _build_interzonal_line_pairs(
    partitions::Partitions{PartitionDC},
    interzonal::Vector{Int},
)::Vector{InterzonalLinePair}
    pairs = InterzonalLinePair[]

    for (cidx, e) in enumerate(interzonal)
        entries = Tuple{Int, Int}[]
        for (zpos, part) in enumerate(partitions.parts)
            part.branch_idx_full_to_zone[e] == 0 && continue
            local_k = findfirst(==(e), part.zone_edges)
            local_k === nothing && continue
            push!(entries, (zpos, local_k))
        end
        length(entries) == 2 || continue
        zi, li = entries[1]
        zj, lj = entries[2]
        push!(pairs, InterzonalLinePair(cidx, e, zi, li, zj, lj))
    end
    return pairs
end

# ------------------------------------------------------------
# Dummy branch pairs
# ------------------------------------------------------------

function _build_dummy_line_pairs(
    partitions::Partitions{PartitionDC},
    zones::Vector{AladinZoneData},
    sys_branch_pairs::Vector{NTuple{4, Int}},
)::Vector{DummyLinePairDC}
    pairs = DummyLinePairDC[]
    isempty(sys_branch_pairs) && return pairs
    sizehint!(pairs, length(sys_branch_pairs))

    for (k, (ga, za, gb, zb)) in enumerate(sys_branch_pairs)
        pa = partitions.parts[za]
        pb = partitions.parts[zb]
        la = pa.branch_idx_full_to_zone[ga]
        lb = pb.branch_idx_full_to_zone[gb]
        (la == 0 || lb == 0) && continue

        zda = zones[za]
        zdb = zones[zb]

        fa = pa.branch.from[la]
        ta = pa.branch.to[la]
        fb = pb.branch.from[lb]
        tb = pb.branch.to[lb]

        theta_a = pa.model[:theta]
        theta_b = pb.model[:theta]
        from_col_i = zda.var_to_col[theta_a[fa, 1]]
        to_col_i = zda.var_to_col[theta_a[ta, 1]]
        from_col_j = zdb.var_to_col[theta_b[fb, 1]]
        to_col_j = zdb.var_to_col[theta_b[tb, 1]]

        b_i = pa.branch.b[la]
        b_j = pb.branch.b[lb]

        push!(
            pairs,
            DummyLinePairDC(
                k, za, from_col_i, to_col_i, b_i,
                zb, from_col_j, to_col_j, b_j,
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
    sys_bus_pairs::Vector{NTuple{3, Int}},
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}},
)::Vector{SplitBusPair}
    pairs = SplitBusPair[]
    isempty(sys_bus_pairs) && return pairs
    sizehint!(pairs, length(sys_bus_pairs))

    for (k, (bid, za, zb)) in enumerate(sys_bus_pairs)
        v_a = get(p_inj_vars[za], bid, nothing)
        v_b = get(p_inj_vars[zb], bid, nothing)
        (v_a === nothing || v_b === nothing) && continue
        col_a = zones[za].var_to_col[v_a]
        col_b = zones[zb].var_to_col[v_b]
        push!(pairs, SplitBusPair(k, bid, za, col_a, zb, col_b))
    end
    return pairs
end
