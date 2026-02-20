function _find_interzonal_branches(full::CaseRaw, ::TieLine)
    full_bus_zones = Int.(full.bus.zone)
    full_bus_ids = Int.(full.bus.bus_i)

    bus_to_zone = Dict{Int, Int}()
    @inbounds for i in eachindex(full_bus_ids)
        bus_to_zone[full_bus_ids[i]] = full_bus_zones[i]
    end

    f_bus = Int.(full.branch.f_bus)
    t_bus = Int.(full.branch.t_bus)
    br_status = full.branch.br_status

    interzonal_branches = Int[]
    @inbounds for i in eachindex(f_bus)
        (br_status[i] == 1) || continue
        if bus_to_zone[f_bus[i]] != bus_to_zone[t_bus[i]]
            push!(interzonal_branches, i)
        end
    end
    return interzonal_branches
end

@inline function branch_in_service_mask(branch::BranchRaw)
    return branch.br_status .== 1
end

@inline function _collect_external_neighbor_buses!(
    external_neighbor_bus_set::Set{Int},
    zone_bus_set::Set{Int},
    full_f_bus::AbstractVector{<:Integer},
    full_t_bus::AbstractVector{<:Integer},
    full_branch_in_service_mask::AbstractVector{Bool},
)
    @inbounds for i in eachindex(full_f_bus)
        full_branch_in_service_mask[i] || continue

        f = Int(full_f_bus[i])
        t = Int(full_t_bus[i])

        f_in = in(f, zone_bus_set)
        t_in = in(t, zone_bus_set)

        if f_in ⊻ t_in
            push!(external_neighbor_bus_set, f_in ? t : f)
        end
    end
    return external_neighbor_bus_set
end

@inline function _build_sub_meta(full::CaseRaw, zone_id::Int)
    return CaseRaw(;
        name = "$(full.name)_zone_$zone_id",
        source_type = full.source_type,
        source_version = full.source_version,
        baseMVA = full.baseMVA,
        T = full.T,
    )
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

@inline function _filter_by_zone_bus(
    bus_vec::AbstractVector{<:Integer},
    zone_bus_set::Set{Int},
)
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
    # generators kept by zone buses
    gen_keep = _filter_by_zone_bus(full.gen.gen_bus, zone_bus_set)
    sub.gen = _slice_gen(full.gen, gen_keep)

    gc = full.gencost
    sub.gencost = GencostRaw{Float64}(
        gc.model[gen_keep],
        gc.startup[gen_keep],
        gc.shutdown[gen_keep],
        gc.ncost[gen_keep],
        gc.cost_coefficients[gen_keep],
    )

    if full.gencost_timeseries !== nothing
        ts = full.gencost_timeseries
        keep = findall(i -> in(Int(ts.gen_id[i]), gen_keep), eachindex(ts.gen_id))
        sub.gencost_timeseries = GencostTimeseriesRaw{Float64}(
            ts.time[keep], ts.gen_id[keep], ts.model[keep],
            ts.startup[keep], ts.shutdown[keep],
            ts.ncost[keep], ts.cost_coefficients[keep],
            ts.tou[keep], ts.feed_in_tariff[keep],
        )
    end

    return nothing
end
