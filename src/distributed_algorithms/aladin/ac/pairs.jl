# ============================================================
# ALADIN-AC — Consensus pair types and builders
# ============================================================

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
# Interzonal line pair (TieLine sectioning)
# ------------------------------------------------------------

struct InterzonalLinePairAC
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

struct SplitBusPairAC
    consensus_idx::Int
    bus_id::Int
    zones::Vector{Int}
    p_inj_cols::Vector{Int}
    q_inj_cols::Vector{Int}
end

# ------------------------------------------------------------
# Dummy line pair (DummyBus sectioning)
# ------------------------------------------------------------

struct DummyLinePairAC
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
    visits = Dict{Int, Vector{Tuple{Int, Int, Int}}}()
    bus_id_by_cidx = Dict{Int, Int}()

    for (zone_pos, part) in enumerate(partitions.parts)
        zone_data = zones[zone_pos]
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            consensus_idx = part.shared_bus_consensus_idx[i]
            consensus_idx == 0 && continue
            theta_col = zone_data.shared_theta_cols[i]
            vm_col = zone_data.shared_vm_cols[i]
            cidx_visits = get!(visits, consensus_idx, Tuple{Int, Int, Int}[])
            push!(cidx_visits, (zone_pos, theta_col, vm_col))
            bus_id_by_cidx[consensus_idx] = part.zone_shared_bus_ids[i]
        end
    end

    pairs = SharedBusPairAC[]
    for consensus_idx in sort!(collect(keys(visits)))
        cidx_visits = visits[consensus_idx]
        if length(cidx_visits) < 2
            @warn "ALADIN-AC pair formation: shared bus $(bus_id_by_cidx[consensus_idx]) " *
                  "(consensus_idx $(consensus_idx)) is present in only " *
                  "$(length(cidx_visits)) zone(s) — coupling DROPPED. Check the partition."
            continue
        end
        sort!(cidx_visits; by = first)
        root_zone, root_theta_col, root_vm_col = cidx_visits[1]
        for k in 2:length(cidx_visits)
            peer_zone, peer_theta_col, peer_vm_col = cidx_visits[k]
            push!(
                pairs,
                SharedBusPairAC(
                    consensus_idx,
                    bus_id_by_cidx[consensus_idx],
                    root_zone, root_theta_col, root_vm_col,
                    peer_zone, peer_theta_col, peer_vm_col,
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
    partitions::Partitions{PartitionAC},
    interzonal::Vector{Int},
)::Vector{InterzonalLinePairAC}
    pairs = InterzonalLinePairAC[]
    for (consensus_idx, line_global) in enumerate(interzonal)
        line_visits = Tuple{Int, Int}[]
        for (zone_pos, part) in enumerate(partitions.parts)
            part.branch_idx_full_to_zone[line_global] == 0 && continue
            local_idx = findfirst(==(line_global), part.zone_edges)
            local_idx === nothing && continue
            push!(line_visits, (zone_pos, local_idx))
        end
        if length(line_visits) != 2
            @warn "ALADIN-AC pair formation: interzonal line $(line_global) " *
                  "(consensus_idx $(consensus_idx)) found in $(length(line_visits)) " *
                  "zone(s), expected 2 — coupling DROPPED. Check the partition."
            continue
        end
        zone_i, local_i = line_visits[1]
        zone_j, local_j = line_visits[2]
        push!(
            pairs,
            InterzonalLinePairAC(
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
    partitions::Partitions{PartitionAC},
    zones::Vector{AladinZoneDataAC},
    sys_branch_pairs::Vector{NTuple{4, Int}},
)::Vector{DummyLinePairAC}
    pairs = DummyLinePairAC[]
    isempty(sys_branch_pairs) && return pairs
    sizehint!(pairs, length(sys_branch_pairs))

    for (consensus_idx, (line_global_a, zone_a, line_global_b, zone_b)) in
        enumerate(sys_branch_pairs)
        part_a = partitions.parts[zone_a]
        part_b = partitions.parts[zone_b]
        local_a = part_a.branch_idx_full_to_zone[line_global_a]
        local_b = part_b.branch_idx_full_to_zone[line_global_b]
        if local_a == 0 || local_b == 0
            @warn "ALADIN-AC pair formation: dummy line pair $(consensus_idx) " *
                  "(lines $(line_global_a) in zone $(zone_a), $(line_global_b) in zone $(zone_b)) " *
                  "has a half missing from its zone (local ids $(local_a), $(local_b)) — " *
                  "coupling DROPPED. Check the partition."
            continue
        end

        zone_data_a = zones[zone_a]
        zone_data_b = zones[zone_b]

        theta_var_a = part_a.local_model[:theta]
        theta_var_b = part_b.local_model[:theta]
        vm_var_a = part_a.local_model[:v_magnitude]
        vm_var_b = part_b.local_model[:v_magnitude]

        from_a = part_a.branch.from[local_a]
        to_a = part_a.branch.to[local_a]
        from_b = part_b.branch.from[local_b]
        to_b = part_b.branch.to[local_b]

        from_theta_col_i = zone_data_a.var_to_col[theta_var_a[from_a, 1]]
        to_theta_col_i = zone_data_a.var_to_col[theta_var_a[to_a, 1]]
        from_vm_col_i = zone_data_a.var_to_col[vm_var_a[from_a, 1]]
        to_vm_col_i = zone_data_a.var_to_col[vm_var_a[to_a, 1]]

        from_theta_col_j = zone_data_b.var_to_col[theta_var_b[from_b, 1]]
        to_theta_col_j = zone_data_b.var_to_col[theta_var_b[to_b, 1]]
        from_vm_col_j = zone_data_b.var_to_col[vm_var_b[from_b, 1]]
        to_vm_col_j = zone_data_b.var_to_col[vm_var_b[to_b, 1]]

        push!(
            pairs,
            DummyLinePairAC(
                consensus_idx,
                zone_a, local_a,
                from_theta_col_i, to_theta_col_i, from_vm_col_i, to_vm_col_i,
                part_a.branch.g[local_a], part_a.branch.b[local_a],
                part_a.branch.b_sh[local_a],
                part_a.branch.inv_tap[local_a], part_a.branch.tap_angle[local_a],
                zone_b, local_b,
                from_theta_col_j, to_theta_col_j, from_vm_col_j, to_vm_col_j,
                part_b.branch.g[local_b], part_b.branch.b[local_b],
                part_b.branch.b_sh[local_b],
                part_b.branch.inv_tap[local_b], part_b.branch.tap_angle[local_b],
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
    bus_groups::Vector{Tuple{Int, Vector{Int}}},
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}},
    q_inj_vars::Vector{Dict{Int, JuMP.VariableRef}},
)::Vector{SplitBusPairAC}
    pairs = SplitBusPairAC[]
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
        if !(ok && length(copy_zones) >= 2)
            @warn "ALADIN-AC pair formation: split bus $(bus_id) (group $(consensus_idx), " *
                  "zones $(group_zones)) has $(length(copy_zones)) usable injection " *
                  "copies$(ok ? "" : " (a zone lacks its p_inj/q_inj variable)") — " *
                  "coupling DROPPED. Check the partition."
            continue
        end
        push!(pairs, SplitBusPairAC(consensus_idx, bus_id, copy_zones, p_cols, q_cols))
    end
    return pairs
end
