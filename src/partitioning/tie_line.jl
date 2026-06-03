# ============================================================
# Sectioning dispatch: Tie Line
# ============================================================

function prepare_distributed_partitions(full::CaseRaw, ::TieLine)
    _require_interzonal_branches(full, "TieLine")

    full_bus_ids = full.bus.bus_i
    full_bus_zones = full.bus.zone
    full_f_bus_ids = full.branch.f_bus
    full_t_bus_ids = full.branch.t_bus
    full_branch_in_service_mask = branch_in_service_mask(full.branch)

    n_full_bus = length(full_bus_ids)
    n_full_branches = length(full_f_bus_ids)
    max_bus_id = n_full_bus == 0 ? 0 : maximum(full_bus_ids)

    zone_ids = sort!(unique(full_bus_zones))
    n_zones = length(zone_ids)
    parts = Vector{CaseRaw}(undef, n_zones)

    zone_internal_membership = zeros(Bool, max_bus_id)
    zone_subprob_membership = zeros(Bool, max_bus_id)
    bus_zone_count = zeros(Int, max_bus_id)

    @inbounds for (i, zone_id) in enumerate(zone_ids)
        fill!(zone_internal_membership, false)
        _mark_internal_membership!(
            zone_internal_membership,
            full_bus_ids,
            full_bus_zones,
            zone_id,
        )

        copyto!(zone_subprob_membership, zone_internal_membership)
        _collect_external_neighbor_buses!(
            zone_subprob_membership,
            zone_internal_membership,
            full_f_bus_ids,
            full_t_bus_ids,
            full_branch_in_service_mask,
        )

        zone_subprob_bus_idx_full =
            _collect_member_positions(full_bus_ids, zone_subprob_membership)

        zone_case = _build_zone_meta(full, zone_id)

        zone_case.bus = _slice_bus(full.bus, zone_subprob_bus_idx_full)
        _mark_virtual_buses!(zone_case.bus, zone_internal_membership)

        zone_gen_idx = _filter_by_zone_bus(full.gen.gen_bus, zone_internal_membership)
        zone_case.gen = _slice_gen(full.gen, zone_gen_idx)
        _slice_gencosts!(zone_case, full, zone_gen_idx)

        zone_branch_idx_in_full = _filter_branches_incident_to_zone(
            full_f_bus_ids,
            full_t_bus_ids,
            full_branch_in_service_mask,
            zone_internal_membership,
        )
        zone_case.branch = _slice_branch(full.branch, zone_branch_idx_in_full)

        if full.storage !== nothing
            zone_storage_idx =
                _filter_by_zone_bus(full.storage.storage_bus, zone_internal_membership)
            zone_case.storage = _slice_storage(full.storage, zone_storage_idx)
        end

        _filter_timeseries!(
            zone_case,
            full,
            zone_internal_membership,
            zone_gen_idx,
        )

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
    @inbounds for b in full_shared_bus_ids
        full_shared_membership[b] = true
    end

    full.dist.full_shared_bus_ids = full_shared_bus_ids

    @inbounds for zc in parts
        zc.dist.zone_shared_bus_ids = [
            bus_id for bus_id in zc.bus.bus_i if full_shared_membership[bus_id]
        ]
    end

    return (zone_ids = zone_ids, parts = parts)
end
