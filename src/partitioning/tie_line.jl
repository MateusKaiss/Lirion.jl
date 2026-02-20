# ==============================================================================
# Sectioning dispatch: Tie Line
# ==============================================================================

# ============================================================
# Distributed partitioning (tie-line sectioning)
# ============================================================

function prepare_distributed_partitions(full::CaseRaw, ::TieLine)
    full_bus_ids = Int.(full.bus.bus_i)
    full_bus_zones = Int.(full.bus.zone)
    zone_ids = sort!(unique(full_bus_zones))

    full_f_bus = Int.(full.branch.f_bus)
    full_t_bus = Int.(full.branch.t_bus)

    full_branch_in_service_mask = branch_in_service_mask(full.branch)

    partitions = Dict{Int, CaseRaw}()

    for zone_id in zone_ids
        zone_bus_idx = findall(==(zone_id), full_bus_zones)
        zone_bus_ids = @views full_bus_ids[zone_bus_idx]
        zone_bus_set = Set{Int}(zone_bus_ids)

        external_neighbor_bus_set = Set{Int}()

        _collect_external_neighbor_buses!(
            external_neighbor_bus_set,
            zone_bus_set,
            full_f_bus,
            full_t_bus,
            full_branch_in_service_mask,
        )

        all_bus_set = union(zone_bus_set, external_neighbor_bus_set)

        sub_bus_idx = findall(b -> in(b, all_bus_set), full_bus_ids)

        sub = _build_sub_meta(full, zone_id)
        sub.bus = _slice_bus(full.bus, sub_bus_idx)
        _mark_virtual_buses!(sub.bus, zone_bus_set)

        gen_keep = _filter_by_zone_bus(Int.(full.gen.gen_bus), zone_bus_set)
        sub.gen = _slice_gen(full.gen, gen_keep)

        if full.prosumer !== nothing
            keep = _filter_by_zone_bus(Int.(full.prosumer.bus), zone_bus_set)
            sub.prosumer = ProsumerRaw{Float64}(
                full.prosumer.bus[keep],
                full.prosumer.generation[keep],
                full.prosumer.demand[keep],
                full.prosumer.q_generation[keep],
                full.prosumer.q_demand[keep],
                full.prosumer.cost_injection[keep],
                full.prosumer.inverter_capacity[keep],
            )
        end

        if full.consumer !== nothing
            keep = _filter_by_zone_bus(Int.(full.consumer.bus), zone_bus_set)
            sub.consumer = ConsumerRaw{Float64}(
                full.consumer.bus[keep],
                full.consumer.demand[keep],
                full.consumer.q_demand[keep],
            )
        end

        if full.storage !== nothing
            keep = _filter_by_zone_bus(Int.(full.storage.storage_bus), zone_bus_set)
            sub.storage = StorageRaw{Float64}(
                full.storage.storage_bus[keep],
                full.storage.ps[keep], full.storage.qs[keep],
                full.storage.energy[keep], full.storage.energy_rating[keep],
                full.storage.charge_rating[keep], full.storage.discharge_rating[keep],
                full.storage.charge_efficiency[keep],
                full.storage.discharge_efficiency[keep],
                full.storage.thermal_rating[keep],
                full.storage.qmin[keep], full.storage.qmax[keep],
                full.storage.r[keep], full.storage.x[keep],
                full.storage.p_loss[keep], full.storage.q_loss[keep],
                full.storage.status[keep],
            )
        end

        zone_branch_idx = findall(
            i ->
                full_branch_in_service_mask[i] &&
                    (
                        in(Int(full_f_bus[i]), zone_bus_set) ||
                        in(Int(full_t_bus[i]), zone_bus_set)
                    ),
            eachindex(full_f_bus),
        )
        sub.branch = _slice_branch(full.branch, zone_branch_idx)

        sub.dist.branch_full = full.branch
        sub.dist.branch_local_indices = zone_branch_idx
        empty!(sub.dist.branch_global_to_local)
        @inbounds for (loc, glob) in enumerate(zone_branch_idx)
            sub.dist.branch_global_to_local[glob] = loc
        end

        _filter_timeseries!(sub, full, zone_bus_set)
        _filter_gencosts!(sub, full, zone_bus_set)

        add_case_index_maps!(sub)

        sub.dist.bus_set = Set(Int.(sub.bus.bus_i))
        partitions[Int(zone_id)] = sub
    end

    # -------------------------
    # compute shared buses across zones
    # -------------------------
    bus_zone_appearance = Dict{Int, Set{Int}}()
    for (zone_id, sub) in partitions, b in Int.(sub.bus.bus_i)
        push!(get!(bus_zone_appearance, b, Set{Int}()), zone_id)
    end

    shared_buses = [b for (b, zs) in bus_zone_appearance if length(zs) > 1]
    shared_set = Set(shared_buses)

    for (_, sub) in partitions
        sub_bus_ids = Int.(sub.bus.bus_i)
        sub.dist.shared_buses_zone = [b for b in sub_bus_ids if in(b, shared_set)]
    end

    return partitions, shared_buses
end
