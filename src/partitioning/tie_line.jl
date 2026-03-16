# ==============================================================================
# Sectioning dispatch: Tie Line
# ==============================================================================

# ============================================================
# Distributed partitioning (tie-line sectioning)
# ============================================================

function prepare_distributed_partitions(full::CaseRaw, ::TieLine)
    full_bus_ids = full.bus.bus_i
    full_bus_zones = full.bus.zone
    full_f_bus_ids = full.branch.f_bus
    full_t_bus_ids = full.branch.t_bus
    full_branch_in_service_mask = branch_in_service_mask(full.branch)

    n_full_branches = length(full_f_bus_ids)

    zone_ids = sort!(unique(full_bus_zones))
    partitions = Dict{Int, CaseRaw}()
    sizehint!(partitions, length(zone_ids))

    for zone_id in zone_ids
        zone_bus_idx = findall(==(zone_id), full_bus_zones)
        zone_bus_ids = @views full_bus_ids[zone_bus_idx]
        zone_internal_bus_set = Set(zone_bus_ids)

        external_neighbor_bus_set = Set{Int}()
        _collect_external_neighbor_buses!(
            external_neighbor_bus_set,
            zone_internal_bus_set,
            full_f_bus_ids,
            full_t_bus_ids,
            full_branch_in_service_mask,
        )

        zone_subprob_bus_set = union(zone_internal_bus_set, external_neighbor_bus_set)
        zone_subprob_bus_idx_full =
            findall(bus_id -> bus_id in zone_subprob_bus_set, full_bus_ids)

        zone_subprob_case = _build_zone_meta(full, zone_id)

        zone_subprob_case.bus = _slice_bus(full.bus, zone_subprob_bus_idx_full)
        _mark_virtual_buses!(zone_subprob_case.bus, zone_internal_bus_set)

        zone_gen_idx = _filter_by_zone_bus(full.gen.gen_bus, zone_internal_bus_set)
        zone_subprob_case.gen = _slice_gen(full.gen, zone_gen_idx)

        _slice_gencosts!(zone_subprob_case, full, zone_gen_idx)

        zone_branch_idx_in_full = _filter_branches_incident_to_zone(
            full_f_bus_ids,
            full_t_bus_ids,
            full_branch_in_service_mask,
            zone_internal_bus_set,
        )
        zone_subprob_case.branch = _slice_branch(full.branch, zone_branch_idx_in_full)

        if full.prosumer !== nothing
            zone_prosumer_idx =
                _filter_by_zone_bus(full.prosumer.bus, zone_internal_bus_set)
            zone_subprob_case.prosumer = _slice_prosumer(full.prosumer, zone_prosumer_idx)
        end

        if full.consumer !== nothing
            zone_consumer_idx =
                _filter_by_zone_bus(full.consumer.bus, zone_internal_bus_set)
            zone_subprob_case.consumer = _slice_consumer(full.consumer, zone_consumer_idx)
        end

        if full.storage !== nothing
            zone_storage_idx =
                _filter_by_zone_bus(full.storage.storage_bus, zone_internal_bus_set)
            zone_subprob_case.storage = _slice_storage(full.storage, zone_storage_idx)
        end

        _filter_timeseries!(
            zone_subprob_case,
            full,
            zone_internal_bus_set,
            zone_gen_idx,
        )

        zone_subprob_case.dist.zone_branch_idx_in_full = zone_branch_idx_in_full

        zone_subprob_case.dist.branch_idx_full_to_zone = fill(0, n_full_branches)
        @inbounds for (local_idx, global_idx) in enumerate(zone_branch_idx_in_full)
            zone_subprob_case.dist.branch_idx_full_to_zone[global_idx] = local_idx
        end

        add_case_index_maps!(zone_subprob_case)

        partitions[zone_id] = zone_subprob_case
    end

    bus_to_zones = Dict{Int, Set{Int}}()

    for (zone_id, zone_case) in partitions
        for bus_id in zone_case.bus.bus_i
            push!(get!(bus_to_zones, bus_id, Set{Int}()), zone_id)
        end
    end

    full_shared_bus_ids = [bus_id for (bus_id, zones) in bus_to_zones if length(zones) > 1]
    full_shared_bus_set = Set(full_shared_bus_ids)

    full.dist.full_shared_bus_ids = full_shared_bus_ids

    for (_, zone_case) in partitions
        zone_case.dist.zone_shared_bus_ids = [
            bus_id for bus_id in zone_case.bus.bus_i if bus_id in full_shared_bus_set
        ]
    end

    return partitions
end
