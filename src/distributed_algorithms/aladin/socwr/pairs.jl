# ============================================================
# ALADIN-SOCWR — Consensus pair types and builders
# ============================================================
#
# ALADIN imposes the coupling `Σ A_i x_i = b` row-by-row through
# explicit pair structures (one per λ component). This file defines
# the three pair types used across coupling variants and the builders
# that populate the pair lists once at workspace construction time.
#
# Per-partition runtime setup (NLP build, quad_static, p_inj/q_inj
# variables) lives in `partition.jl`.
#
# Note: SOCWR has no DummyLinePair — DummyBus sectioning is
# explicitly rejected at ALADINConfig time (see `src/types.jl`).

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
# Shared edge pair
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
# Split bus pair
# ------------------------------------------------------------

struct SplitBusPairSOCWR
    consensus_idx::Int
    bus_id::Int
    zone_i::Int
    pinj_col_i::Int
    qinj_col_i::Int
    zone_j::Int
    pinj_col_j::Int
    qinj_col_j::Int
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
    acc = Dict{Int, Vector{Tuple{Int, Int}}}()
    bus_of = Dict{Int, Int}()

    for (zpos, part) in enumerate(partitions.parts)
        zdata = zones[zpos]
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            cidx = part.shared_bus_consensus_idx[i]
            cidx == 0 && continue
            col = zdata.shared_v2_cols[i]
            v = get!(acc, cidx, Tuple{Int, Int}[])
            push!(v, (zpos, col))
            bus_of[cidx] = part.zone_shared_bus_ids[i]
        end
    end

    pairs = SharedBusPairSOCWR[]
    for cidx in sort!(collect(keys(acc)))
        entries = acc[cidx]
        length(entries) >= 2 || continue
        # Spanning tree (star): pair every other zone with the lowest-zpos
        # zone as root. N zones → N-1 pairs, LICQ preserved.
        sort!(entries; by = first)
        zr, cr = entries[1]
        for k in 2:length(entries)
            zk, ck = entries[k]
            push!(pairs, SharedBusPairSOCWR(cidx, bus_of[cidx], zr, cr, zk, ck))
        end
    end
    return pairs
end

# ------------------------------------------------------------
# Shared edge pairs
# ------------------------------------------------------------

function _build_interzonal_line_pairs(
    partitions::Partitions{PartitionSOCWR},
    interzonal::Vector{Int},
)::Vector{InterzonalLinePairSOCWR}
    pairs = InterzonalLinePairSOCWR[]
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
        push!(pairs, InterzonalLinePairSOCWR(cidx, e, zi, li, zj, lj))
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

    for (k, (bid, za, zb)) in enumerate(sys_bus_pairs)
        vpa = get(p_inj_vars[za], bid, nothing)
        vqa = get(q_inj_vars[za], bid, nothing)
        vpb = get(p_inj_vars[zb], bid, nothing)
        vqb = get(q_inj_vars[zb], bid, nothing)
        (vpa === nothing || vqa === nothing ||
         vpb === nothing || vqb === nothing) && continue
        pcol_a = zones[za].var_to_col[vpa]
        qcol_a = zones[za].var_to_col[vqa]
        pcol_b = zones[zb].var_to_col[vpb]
        qcol_b = zones[zb].var_to_col[vqb]
        push!(pairs, SplitBusPairSOCWR(
            k, bid, za, pcol_a, qcol_a, zb, pcol_b, qcol_b,
        ))
    end
    return pairs
end
