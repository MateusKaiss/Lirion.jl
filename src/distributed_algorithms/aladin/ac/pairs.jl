# ============================================================
# ALADIN-AC — Consensus pair types and builders
# ============================================================
#
# ALADIN imposes the coupling `Σ A_i x_i = b` row-by-row through
# explicit pair structures (one per λ component). This file defines
# the four pair types used across coupling variants and the builders
# that populate the pair lists once at workspace construction time.
#
# Per-partition runtime setup (NLP build, quad_static, p_inj/q_inj
# variables) lives in `partition.jl`.

# ============================================================
# Pair types (consensus structures)
# ============================================================

# ------------------------------------------------------------
# Shared bus pair
# ------------------------------------------------------------

struct SharedBusPairAC
    consensus_idx::Int
    bus_id::Int
    zone_i::Int
    theta_col_i::Int
    vm_col_i::Int
    zone_j::Int
    theta_col_j::Int
    vm_col_j::Int
end

# ------------------------------------------------------------
# Shared edge pair
# ------------------------------------------------------------

struct SharedEdgePairAC
    consensus_idx::Int
    edge_global::Int
    zone_i::Int
    edge_local_i::Int
    zone_j::Int
    edge_local_j::Int
end

# ------------------------------------------------------------
# Split bus pair
# ------------------------------------------------------------

struct SplitBusPairAC
    consensus_idx::Int
    bus_id::Int
    zone_i::Int
    pinj_col_i::Int
    qinj_col_i::Int
    zone_j::Int
    pinj_col_j::Int
    qinj_col_j::Int
end

# ------------------------------------------------------------
# Dummy branch pair
# ------------------------------------------------------------

struct DummyBranchPairAC
    consensus_idx::Int
    zone_i::Int
    local_i::Int
    from_theta_col_i::Int
    to_theta_col_i::Int
    from_vm_col_i::Int
    to_vm_col_i::Int
    g_i::Float64
    b_i::Float64
    b_sh_i::Float64
    inv_tap_i::Float64
    tap_angle_i::Float64
    zone_j::Int
    local_j::Int
    from_theta_col_j::Int
    to_theta_col_j::Int
    from_vm_col_j::Int
    to_vm_col_j::Int
    g_j::Float64
    b_j::Float64
    b_sh_j::Float64
    inv_tap_j::Float64
    tap_angle_j::Float64
end

# ============================================================
# Pair builders
# ============================================================

# ------------------------------------------------------------
# Shared bus pairs
# ------------------------------------------------------------

function _build_shared_bus_pairs(
    partitions::Partitions{PartitionAC},
    zones::Vector{AladinZoneDataAC},
)::Vector{SharedBusPairAC}
    acc = Dict{Int, Vector{Tuple{Int, Int, Int}}}()
    bus_of = Dict{Int, Int}()

    for (zpos, part) in enumerate(partitions.parts)
        zdata = zones[zpos]
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            cidx = part.shared_bus_consensus_idx[i]
            cidx == 0 && continue
            tcol = zdata.shared_theta_cols[i]
            vcol = zdata.shared_vm_cols[i]
            v = get!(acc, cidx, Tuple{Int, Int, Int}[])
            push!(v, (zpos, tcol, vcol))
            bus_of[cidx] = part.zone_shared_bus_ids[i]
        end
    end

    pairs = SharedBusPairAC[]
    for cidx in sort!(collect(keys(acc)))
        entries = acc[cidx]
        length(entries) >= 2 || continue
        # Spanning tree (star): pair every other zone with the lowest-zpos
        # zone as root. N zones → N-1 pairs, LICQ preserved.
        sort!(entries; by = first)
        zr, tr, vr = entries[1]
        for k in 2:length(entries)
            zk, tk, vk = entries[k]
            push!(pairs, SharedBusPairAC(cidx, bus_of[cidx], zr, tr, vr, zk, tk, vk))
        end
    end
    return pairs
end

# ------------------------------------------------------------
# Shared edge pairs
# ------------------------------------------------------------

function _build_shared_edge_pairs(
    partitions::Partitions{PartitionAC},
    interzonal::Vector{Int},
)::Vector{SharedEdgePairAC}
    pairs = SharedEdgePairAC[]
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
        push!(pairs, SharedEdgePairAC(cidx, e, zi, li, zj, lj))
    end
    return pairs
end

# ------------------------------------------------------------
# Dummy branch pairs
# ------------------------------------------------------------

function _build_dummy_branch_pairs(
    partitions::Partitions{PartitionAC},
    zones::Vector{AladinZoneDataAC},
    sys_branch_pairs::Vector{NTuple{4, Int}},
)::Vector{DummyBranchPairAC}
    pairs = DummyBranchPairAC[]
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

        theta_a = pa.model[:theta]
        theta_b = pb.model[:theta]
        vm_a = pa.model[:v_magnitude]
        vm_b = pb.model[:v_magnitude]

        fa = pa.branch.from[la]
        ta = pa.branch.to[la]
        fb = pb.branch.from[lb]
        tb = pb.branch.to[lb]

        from_theta_col_i = zda.var_to_col[theta_a[fa, 1]]
        to_theta_col_i = zda.var_to_col[theta_a[ta, 1]]
        from_vm_col_i = zda.var_to_col[vm_a[fa, 1]]
        to_vm_col_i = zda.var_to_col[vm_a[ta, 1]]

        from_theta_col_j = zdb.var_to_col[theta_b[fb, 1]]
        to_theta_col_j = zdb.var_to_col[theta_b[tb, 1]]
        from_vm_col_j = zdb.var_to_col[vm_b[fb, 1]]
        to_vm_col_j = zdb.var_to_col[vm_b[tb, 1]]

        push!(
            pairs,
            DummyBranchPairAC(
                k,
                za, la,
                from_theta_col_i, to_theta_col_i, from_vm_col_i, to_vm_col_i,
                pa.branch.g[la], pa.branch.b[la], pa.branch.b_sh[la],
                pa.branch.inv_tap[la], pa.branch.tap_angle[la],
                zb, lb,
                from_theta_col_j, to_theta_col_j, from_vm_col_j, to_vm_col_j,
                pb.branch.g[lb], pb.branch.b[lb], pb.branch.b_sh[lb],
                pb.branch.inv_tap[lb], pb.branch.tap_angle[lb],
            ),
        )
    end
    return pairs
end

# ------------------------------------------------------------
# Split bus pairs
# ------------------------------------------------------------

function _build_split_bus_pairs(
    zones::Vector{AladinZoneDataAC},
    sys_bus_pairs::Vector{NTuple{3, Int}},
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}},
    q_inj_vars::Vector{Dict{Int, JuMP.VariableRef}},
)::Vector{SplitBusPairAC}
    pairs = SplitBusPairAC[]
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
        push!(pairs, SplitBusPairAC(
            k, bid, za, pcol_a, qcol_a, zb, pcol_b, qcol_b,
        ))
    end
    return pairs
end
