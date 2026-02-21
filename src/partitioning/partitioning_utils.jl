
# ==============================================================================
# Interzonal branches (tie lines) finder
# ==============================================================================

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

# ==============================================================================
# Helpers
# ==============================================================================

@inline function branch_in_service_mask(branch::BranchRaw)
    return branch.br_status .== 1
end

@inline function _collect_external_neighbor_buses!(
    external_neighbor_bus_set::Set{Int},
    zone_bus_set::Set{Int},
    full_f_bus_ids::AbstractVector{Int},
    full_t_bus_ids::AbstractVector{Int},
    full_branch_in_service_mask::AbstractVector{Bool},
)
    @inbounds for i in eachindex(full_f_bus_ids)
        full_branch_in_service_mask[i] || continue

        full_f_bus_id = full_f_bus_ids[i]
        full_t_bus_id = full_t_bus_ids[i]

        f_in_zone = in(full_f_bus_id, zone_bus_set)
        t_in_zone = in(full_t_bus_id, zone_bus_set)

        if f_in_zone ⊻ t_in_zone # XOR Logic - A ⊻ B is true only if exactly one of them is true
            push!(external_neighbor_bus_set, f_in_zone ? full_t_bus_id : full_f_bus_id)
        end
    end
    return external_neighbor_bus_set
end

@inline function _build_zone_meta(full::CaseRaw, zone_id::Int)
    return CaseRaw(;
        name = "$(full.name)_zone_$zone_id",
        source_type = full.source_type,
        source_version = full.source_version,
        baseMVA = full.baseMVA,
        T = full.T,
    )
end

@inline function _slice_bus(full_bus::BusRaw{Float64}, keep::Vector{Int})
    return BusRaw{Float64}(
        full_bus.bus_i[keep],
        full_bus.bus_type[keep],
        full_bus.pd[keep], full_bus.qd[keep],
        full_bus.gs[keep], full_bus.bs[keep],
        full_bus.area[keep],
        full_bus.vm[keep], full_bus.va[keep],
        full_bus.base_kv[keep],
        full_bus.zone[keep],
        full_bus.vmax[keep], full_bus.vmin[keep],
        full_bus.lam_p[keep], full_bus.lam_q[keep],
        full_bus.mu_vmax[keep], full_bus.mu_vmin[keep],
        full_bus.is_virtual[keep],
    )
end

@inline function _mark_virtual_buses!(bus::BusRaw{Float64}, zone_bus_set::Set{Int})
    @inbounds for i in eachindex(bus.bus_i)
        bus_id = bus.bus_i[i]
        bus.is_virtual[i] = !in(bus_id, zone_bus_set)

        if bus.is_virtual[i]
            bus.pd[i] = 0.0
            bus.qd[i] = 0.0
            bus.gs[i] = 0.0
            bus.bs[i] = 0.0
            # bus.bus_type[i] = 1
        end
    end
    return nothing
end

@inline function _filter_by_zone_bus(
    component_bus_ids::AbstractVector{<:Integer},
    zone_internal_bus_set::Set{Int},
)
    return findall(
        idx -> component_bus_ids[idx] in zone_internal_bus_set,
        eachindex(component_bus_ids),
    )
end

@inline function _slice_gen(full_gen::GenRaw{Float64}, keep::Vector{Int})
    return GenRaw{Float64}(
        full_gen.gen_bus[keep],
        full_gen.pg[keep], full_gen.qg[keep],
        full_gen.qmax[keep], full_gen.qmin[keep],
        full_gen.vg[keep],
        full_gen.mbase[keep],
        full_gen.gen_status[keep],
        full_gen.pmax[keep], full_gen.pmin[keep],
        full_gen.pc1[keep],
        full_gen.pc2[keep],
        full_gen.qc1min[keep], full_gen.qc1max[keep],
        full_gen.qc2min[keep], full_gen.qc2max[keep],
        full_gen.ramp_agc[keep],
        full_gen.ramp_10[keep],
        full_gen.ramp_30[keep],
        full_gen.ramp_q[keep],
        full_gen.apf[keep],
        full_gen.mu_pmax[keep], full_gen.mu_pmin[keep],
        full_gen.mu_qmax[keep], full_gen.mu_qmin[keep],
        full_gen.plant_type[keep],
        full_gen.curtailment_cost[keep],
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

@inline function _filter_timeseries!(
    zone_case::CaseRaw,
    full::CaseRaw,
    zone_bus_set::Set{Int},
)
    if full.bus_timeseries !== nothing
        ts = full.bus_timeseries
        zone_ts_idx = _filter_by_zone_bus(ts.bus, zone_bus_set)
        zone_case.bus_timeseries = BusTimeseriesRaw{Float64}(
            ts.time[zone_ts_idx], ts.bus[zone_ts_idx],
            ts.pd[zone_ts_idx], ts.qd[zone_ts_idx], ts.gs[zone_ts_idx],
            ts.bs[zone_ts_idx],
        )
    end

    if full.prosumer_timeseries !== nothing
        ts = full.prosumer_timeseries
        zone_ts_idx = _filter_by_zone_bus(ts.bus, zone_bus_set)
        zone_case.prosumer_timeseries = ProsumerTimeseriesRaw{Float64}(
            ts.time[zone_ts_idx], ts.bus[zone_ts_idx],
            ts.p_pv[zone_ts_idx], ts.p_load[zone_ts_idx], ts.p_curtail[zone_ts_idx],
        )
    end

    if full.consumer_timeseries !== nothing
        ts = full.consumer_timeseries
        zone_ts_idx = _filter_by_zone_bus(ts.bus, zone_bus_set)
        zone_case.consumer_timeseries = ConsumerTimeseriesRaw{Float64}(
            ts.time[zone_ts_idx], ts.bus[zone_ts_idx],
            ts.demand[zone_ts_idx], ts.q_demand[zone_ts_idx],
        )
    end

    return nothing
end

@inline function _filter_gencosts!(
    zone_case::CaseRaw,
    full::CaseRaw,
    zone_bus_set::Set{Int},
)
    zone_gen_idx = _filter_by_zone_bus(full.gen.gen_bus, zone_bus_set)
    zone_case.gen = _slice_gen(full.gen, zone_gen_idx)

    gc = full.gencost
    zone_case.gencost = GencostRaw{Float64}(
        gc.model[zone_gen_idx],
        gc.startup[zone_gen_idx],
        gc.shutdown[zone_gen_idx],
        gc.ncost[zone_gen_idx],
        gc.cost_coefficients[zone_gen_idx],
    )

    if full.gencost_timeseries !== nothing
        ts = full.gencost_timeseries
        zone_gc_ts_idx =
            findall(i -> in(Int(ts.gen_id[i]), zone_gen_idx), eachindex(ts.gen_id))
        zone_case.gencost_timeseries = GencostTimeseriesRaw{Float64}(
            ts.time[zone_gc_ts_idx], ts.gen_id[zone_gc_ts_idx],
            ts.model[zone_gc_ts_idx],
            ts.startup[zone_gc_ts_idx], ts.shutdown[zone_gc_ts_idx],
            ts.ncost[zone_gc_ts_idx], ts.cost_coefficients[zone_gc_ts_idx],
            ts.tou[zone_gc_ts_idx], ts.feed_in_tariff[zone_gc_ts_idx],
        )
    end

    return nothing
end
