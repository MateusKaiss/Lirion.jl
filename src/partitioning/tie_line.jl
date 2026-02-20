
# ==============================================================================
# Sectioning dispatch: Tie Line
# ==============================================================================

# ============================================================
# Distributed partitioning (tie-line sectioning)
# ============================================================

function prepare_distributed_partitions(full::CaseRaw, ::TieLine)
    @assert full.bus !== nothing && full.branch !== nothing

    bus_ids = Int.(full.bus.bus_i)
    zones = Int.(full.bus.zone)
    unique_zones = sort!(unique(zones))

    partitions = Dict{Int, CaseRaw}()

    f_bus_full = Int.(full.branch.f_bus)
    t_bus_full = Int.(full.branch.t_bus)
    branch_status_mask = branch_in_service_mask(full.branch)

    for zone in unique_zones
        zone_mask = findall(==(zone), zones)

        zone_bus_ids = @views bus_ids[zone_mask]
        zone_bus_set = Set{Int}(zone_bus_ids)

        extra = Set{Int}()
        _find_extra_neighbor_buses!(
            extra,
            zone_bus_set,
            f_bus_full,
            t_bus_full,
            branch_status_mask,
        )

        all_bus_set = union(zone_bus_set, extra)
        bus_mask = findall(b -> in(b, all_bus_set), bus_ids)

        sub = _build_sub_meta(full, zone)
        sub.bus = _slice_bus(full.bus, bus_mask)
        _mark_virtual_buses!(sub.bus, zone_bus_set)

        if full.gen !== nothing
            keep = _filter_by_zone_bus(Int.(full.gen.gen_bus), zone_bus_set)
            sub.gen = _slice_gen(full.gen, keep)
        end

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

        keep_all = findall(
            i ->
                branch_status_mask[i] &&
                    (
                        in(Int(f_bus_full[i]), zone_bus_set) ||
                        in(Int(t_bus_full[i]), zone_bus_set)
                    ),
            eachindex(f_bus_full),
        )
        sub.branch = _slice_branch(full.branch, keep_all)

        sub.dist.branch_full = full.branch
        sub.dist.branch_local_indices = keep_all
        empty!(sub.dist.branch_global_to_local)
        @inbounds for (loc, glob) in enumerate(keep_all)
            sub.dist.branch_global_to_local[glob] = loc
        end

        _filter_timeseries!(sub, full, zone_bus_set)
        _filter_gencosts!(sub, full, zone_bus_set)

        add_case_index_maps!(sub)

        sub.dist.bus_set = Set(Int.(sub.bus.bus_i))
        partitions[Int(zone)] = sub
    end

    bus_zone_appearance = Dict{Int, Set{Int}}()
    for (zone_id, sub) in partitions, b in Int.(sub.bus.bus_i)
        push!(get!(bus_zone_appearance, b, Set{Int}()), zone_id)
    end
    shared_buses = [b for (b, zs) in bus_zone_appearance if length(zs) > 1]
    shared_set = Set(shared_buses)

    for (_, sub) in partitions
        local_b = Int.(sub.bus.bus_i)
        sub.dist.shared_buses_zone = [b for b in local_b if in(b, shared_set)]
    end

    return partitions, shared_buses
end
