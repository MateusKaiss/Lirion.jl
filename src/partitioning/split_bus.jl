# ============================================================
# Sectioning dispatch: Split Bus
# ============================================================

function prepare_distributed_partitions(full::CaseRaw, ::SplitBus)
    _require_interzonal_branches(full, "SplitBus")

    full_bus_ids = full.bus.bus_i
    full_bus_zones = full.bus.zone
    full_f_bus_ids = full.branch.f_bus
    full_t_bus_ids = full.branch.t_bus
    full_branch_in_service_mask = branch_in_service_mask(full.branch)

    n_full_bus = length(full_bus_ids)
    n_full_branches = length(full_f_bus_ids)
    max_bus_id = n_full_bus == 0 ? 0 : maximum(full_bus_ids)
    bus_to_zone = zeros(Int, max_bus_id)
    @inbounds for i in 1:n_full_bus
        bus_to_zone[full_bus_ids[i]] = full_bus_zones[i]
    end

    zone_ids = sort!(unique(full_bus_zones))
    n_zones = length(zone_ids)

    branch_owner, split_set = _compute_branch_owner_and_split_set(
        full_f_bus_ids,
        full_t_bus_ids,
        full_branch_in_service_mask,
        bus_to_zone,
        n_full_branches,
    )

    zone_ghost_buses = Dict{Int, Vector{Int}}()
    @inbounds for zone_id in zone_ids
        zone_ghost_buses[zone_id] = Int[]
    end
    @inbounds for (fbus, zb) in split_set
        push!(zone_ghost_buses[zb], fbus)
    end

    parts = Vector{CaseRaw}(undef, n_zones)
    bus_zone_count = zeros(Int, max_bus_id)
    zone_internal_membership = zeros(Bool, max_bus_id)
    zone_subprob_membership = zeros(Bool, max_bus_id)

    @inbounds for (i, zone_id) in enumerate(zone_ids)
        fill!(zone_internal_membership, false)
        _mark_internal_membership!(zone_internal_membership, full_bus_ids, full_bus_zones, zone_id)

        copyto!(zone_subprob_membership, zone_internal_membership)
        for fbus in zone_ghost_buses[zone_id]
            zone_subprob_membership[fbus] = true
        end

        zone_subprob_bus_idx_full = _collect_member_positions(full_bus_ids, zone_subprob_membership)

        zone_case = _build_zone_meta(full, zone_id)
        zone_case.bus = _slice_bus(full.bus, zone_subprob_bus_idx_full)
        _zero_clone_loads!(zone_case.bus, zone_internal_membership)

        zone_gen_idx = _filter_by_zone_bus(full.gen.gen_bus, zone_internal_membership)
        zone_case.gen = _slice_gen(full.gen, zone_gen_idx)
        _slice_gencosts!(zone_case, full, zone_gen_idx)
        
        zone_branch_idx_in_full = _filter_branches_owned_by_zone(
            full_branch_in_service_mask,
            branch_owner,
            zone_id,
        )
        zone_case.branch = _slice_branch(full.branch, zone_branch_idx_in_full)

        if full.storage !== nothing
            zone_storage_idx =
                _filter_by_zone_bus(full.storage.storage_bus, zone_internal_membership)
            zone_case.storage = _slice_storage(full.storage, zone_storage_idx)
        end

        _filter_timeseries!(zone_case, full, zone_internal_membership, zone_gen_idx)

        zone_case.dist.zone_branch_idx_in_full = zone_branch_idx_in_full
        zone_case.dist.branch_idx_full_to_zone = zeros(Int, n_full_branches)
        for (local_idx, global_idx) in enumerate(zone_branch_idx_in_full)
            zone_case.dist.branch_idx_full_to_zone[global_idx] = local_idx
        end

        add_case_index_maps!(zone_case)

        for b in zone_case.bus.bus_i
            bus_zone_count[b] += 1
        end

        parts[i] = zone_case
    end

    full_shared_bus_ids = Int[]
    @inbounds for b in 1:max_bus_id
        bus_zone_count[b] > 1 && push!(full_shared_bus_ids, b)
    end
    full_shared_membership = zeros(Bool, max_bus_id)
    for b in full_shared_bus_ids
        full_shared_membership[b] = true
    end
    full.dist.full_shared_bus_ids = full_shared_bus_ids

    for zone_case in parts
        zone_case.dist.zone_shared_bus_ids =
            [b for b in zone_case.bus.bus_i if full_shared_membership[b]]
    end
    neighbor_zones_of = Dict{Int, Vector{Int}}()
    for (fbus, zb) in split_set
        push!(get!(neighbor_zones_of, fbus, Int[]), zb)
    end

    zone_id_to_pos = _compute_zone_id_to_pos(zone_ids)

    bus_pairs = NTuple{3, Int}[]
    sizehint!(bus_pairs, length(split_set))
    for (fbus, zb) in split_set
        za = bus_to_zone[fbus]
        push!(
            bus_pairs,
            (fbus, zone_id_to_pos[za], zone_id_to_pos[zb]),
        )
    end
    full.dist.bus_pairs = bus_pairs

    bus_groups = Tuple{Int, Vector{Int}}[]
    sizehint!(bus_groups, length(neighbor_zones_of))
    for (fbus, zbs) in neighbor_zones_of
        za = bus_to_zone[fbus]
        zones_pos = Int[zone_id_to_pos[za]]
        for zb in sort!(unique(zbs))
            push!(zones_pos, zone_id_to_pos[zb])
        end
        push!(bus_groups, (fbus, zones_pos))
    end
    sort!(bus_groups; by = g -> g[1])
    full.dist.bus_groups = bus_groups

    return (zone_ids = zone_ids, parts = parts)
end
