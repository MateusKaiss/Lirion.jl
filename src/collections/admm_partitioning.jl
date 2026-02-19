# ============================================================
# ADMM preprocessing
# ============================================================

function _find_interzonal_branches(full::CaseRaw)
    @assert full.bus !== nothing && full.branch !== nothing

    zones = full.bus.zone
    bus_ids = full.bus.bus_i
    bus_to_zone = Dict(bus_ids[i] => zones[i] for i in eachindex(bus_ids))

    f_bus = full.branch.f_bus
    t_bus = full.branch.t_bus
    br_status = full.branch.br_status

    interzonal_branches = Int[]
    for i in eachindex(f_bus)
        if br_status[i] == 1 && bus_to_zone[f_bus[i]] != bus_to_zone[t_bus[i]]
            push!(interzonal_branches, i)
        end
    end
    return interzonal_branches
end

@inline function branch_in_service_mask(branch::BranchRaw)
    return branch.br_status .== 1
end

@inline function _build_sub_meta(full::CaseRaw, zone::Int)
    sub = CaseRaw(;
        name = "$(full.name)_zone_$zone",
        source_type = full.source_type,
        source_version = full.source_version,
        baseMVA = full.baseMVA,
        T = full.T,
    )
    return sub
end

@inline function _find_extra_neighbor_buses!(
    extra::Set{Int},
    zone_bus_set::Set{Int},
    f_bus_full::AbstractVector{<:Integer},
    t_bus_full::AbstractVector{<:Integer},
    branch_status_mask::AbstractVector{Bool},
)
    @inbounds for i in eachindex(f_bus_full)
        branch_status_mask[i] || continue
        f_in = in(Int(f_bus_full[i]), zone_bus_set)
        t_in = in(Int(t_bus_full[i]), zone_bus_set)
        if f_in ⊻ t_in
            push!(extra, f_in ? Int(t_bus_full[i]) : Int(f_bus_full[i]))
        end
    end
    return extra
end

@inline function _mark_virtual_buses!(bus::BusRaw{Float64}, zone_bus_set::Set{Int})
    @inbounds for i in eachindex(bus.bus_i)
        b = bus.bus_i[i]
        bus.is_virtual[i] = !in(b, zone_bus_set)
        if bus.is_virtual[i]
            bus.pd[i] = 0.0
            bus.qd[i] = 0.0
            bus.gs[i] = 0.0
            bus.bs[i] = 0.0
            bus.bus_type[i] = 1
        end
    end
    return nothing
end

@inline function _filter_by_zone_bus(bus_vec::Vector{Int}, zone_bus_set::Set{Int})
    return findall(i -> in(Int(bus_vec[i]), zone_bus_set), eachindex(bus_vec))
end

@inline function _slice_bus(fullbus::BusRaw{Float64}, keep::Vector{Int})
    return BusRaw{Float64}(
        fullbus.bus_i[keep],
        fullbus.bus_type[keep],
        fullbus.pd[keep], fullbus.qd[keep],
        fullbus.gs[keep], fullbus.bs[keep],
        fullbus.area[keep],
        fullbus.vm[keep], fullbus.va[keep],
        fullbus.base_kv[keep],
        fullbus.zone[keep],
        fullbus.vmax[keep], fullbus.vmin[keep],
        fullbus.lam_p[keep], fullbus.lam_q[keep],
        fullbus.mu_vmax[keep], fullbus.mu_vmin[keep],
        fullbus.is_virtual[keep],
    )
end

@inline function _slice_gen(fullgen::GenRaw{Float64}, keep::Vector{Int})
    return GenRaw{Float64}(
        fullgen.gen_bus[keep],
        fullgen.pg[keep], fullgen.qg[keep],
        fullgen.qmax[keep], fullgen.qmin[keep],
        fullgen.vg[keep],
        fullgen.mbase[keep],
        fullgen.gen_status[keep],
        fullgen.pmax[keep], fullgen.pmin[keep],
        fullgen.pc1[keep],
        fullgen.pc2[keep],
        fullgen.qc1min[keep], fullgen.qc1max[keep],
        fullgen.qc2min[keep], fullgen.qc2max[keep],
        fullgen.ramp_agc[keep],
        fullgen.ramp_10[keep],
        fullgen.ramp_30[keep],
        fullgen.ramp_q[keep],
        fullgen.apf[keep],
        fullgen.mu_pmax[keep], fullgen.mu_pmin[keep],
        fullgen.mu_qmax[keep], fullgen.mu_qmin[keep],
        fullgen.plant_type[keep],
        fullgen.curtailment_cost[keep],
    )
end

@inline function _slice_branch(fullbr::BranchRaw{Float64}, keep::Vector{Int})
    return BranchRaw{Float64}(
        fullbr.f_bus[keep],
        fullbr.t_bus[keep],
        fullbr.br_r[keep], fullbr.br_x[keep],
        fullbr.br_b[keep],
        fullbr.rate_a[keep],
        fullbr.rate_b[keep],
        fullbr.rate_c[keep],
        fullbr.tap[keep], fullbr.shift[keep],
        fullbr.br_status[keep],
        fullbr.angmin[keep], fullbr.angmax[keep],
        fullbr.pf[keep], fullbr.qf[keep],
        fullbr.pt[keep], fullbr.qt[keep],
        fullbr.mu_sf[keep], fullbr.mu_st[keep],
        fullbr.mu_angmin[keep], fullbr.mu_angmax[keep],
    )
end

@inline function _filter_timeseries!(sub::CaseRaw, full::CaseRaw, zone_bus_set::Set{Int})
    if full.bus_timeseries !== nothing
        ts = full.bus_timeseries
        keep = _filter_by_zone_bus(ts.bus, zone_bus_set)
        sub.bus_timeseries = BusTimeseriesRaw{Float64}(
            ts.time[keep], ts.bus[keep],
            ts.pd[keep], ts.qd[keep], ts.gs[keep], ts.bs[keep],
        )
    end
    if full.prosumer_timeseries !== nothing
        ts = full.prosumer_timeseries
        keep = _filter_by_zone_bus(ts.bus, zone_bus_set)
        sub.prosumer_timeseries = ProsumerTimeseriesRaw{Float64}(
            ts.time[keep], ts.bus[keep],
            ts.p_pv[keep], ts.p_load[keep], ts.p_curtail[keep],
        )
    end
    if full.consumer_timeseries !== nothing
        ts = full.consumer_timeseries
        keep = _filter_by_zone_bus(ts.bus, zone_bus_set)
        sub.consumer_timeseries = ConsumerTimeseriesRaw{Float64}(
            ts.time[keep], ts.bus[keep],
            ts.demand[keep], ts.q_demand[keep],
        )
    end
    return nothing
end

@inline function _filter_gencosts!(sub::CaseRaw, full::CaseRaw, zone_bus_set::Set{Int})
    if full.gen === nothing
        sub.gen = nothing
        sub.gencost = nothing
        sub.gencost_timeseries = nothing
        return nothing
    end

    gen_keep = _filter_by_zone_bus(full.gen.gen_bus, zone_bus_set)
    sub.gen =
        isempty(gen_keep) ? _slice_gen(full.gen, Int[]) : _slice_gen(full.gen, gen_keep)

    if full.gencost !== nothing
        gk = gen_keep
        gc = full.gencost
        sub.gencost = GencostRaw{Float64}(
            gc.model[gk],
            gc.startup[gk],
            gc.shutdown[gk],
            gc.ncost[gk],
            gc.cost_coefficients[gk],
        )
    end

    if full.gencost_timeseries !== nothing
        ts = full.gencost_timeseries
        keep = findall(i -> in(ts.gen_id[i], gen_keep), eachindex(ts.gen_id))
        sub.gencost_timeseries = GencostTimeseriesRaw{Float64}(
            ts.time[keep], ts.gen_id[keep], ts.model[keep], ts.startup[keep],
            ts.shutdown[keep],
            ts.ncost[keep], ts.cost_coefficients[keep], ts.tou[keep],
            ts.feed_in_tariff[keep],
        )
    end

    return nothing
end

# ============================================================
# Prepare ADMM partitions
# ============================================================

function prepare_admm_partitions(full::CaseRaw)
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
            keep = _filter_by_zone_bus(full.gen.gen_bus, zone_bus_set)
            sub.gen = _slice_gen(full.gen, keep)
        end
        if full.prosumer !== nothing
            keep = _filter_by_zone_bus(full.prosumer.bus, zone_bus_set)
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
            keep = _filter_by_zone_bus(full.consumer.bus, zone_bus_set)
            sub.consumer = ConsumerRaw{Float64}(
                full.consumer.bus[keep],
                full.consumer.demand[keep],
                full.consumer.q_demand[keep],
            )
        end
        if full.storage !== nothing
            keep = _filter_by_zone_bus(full.storage.storage_bus, zone_bus_set)
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
                    (in(f_bus_full[i], zone_bus_set) || in(t_bus_full[i], zone_bus_set)),
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
