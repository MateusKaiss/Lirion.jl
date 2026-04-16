# ==============================================================================
# Interzonal branches (tie lines) finder
# ==============================================================================

# DummyBus sectioning does not have classical tie-lines (each split half is
# internal to one zone, the dummy bus is shared). Branch pairs are carried
# separately via `DistributedRawMeta.branch_pairs`.
_find_interzonal_branches(::CaseRaw, ::DummyBus) = Int[]

function _find_interzonal_branches(full::CaseRaw, ::TieLine)
    full_bus_ids = full.bus.bus_i
    full_bus_zones = full.bus.zone
    f_bus = full.branch.f_bus
    t_bus = full.branch.t_bus
    br_status = full.branch.br_status

    n_bus = length(full_bus_ids)
    max_bus_id = n_bus == 0 ? 0 : maximum(full_bus_ids)

    # Dense map: bus_id -> zone id (0 if absent).
    bus_to_zone = zeros(Int, max_bus_id)
    @inbounds for i in 1:n_bus
        bus_to_zone[full_bus_ids[i]] = full_bus_zones[i]
    end

    interzonal_branches = Int[]
    sizehint!(interzonal_branches, length(f_bus))

    @inbounds for i in eachindex(f_bus)
        br_status[i] == 1 || continue
        fi = f_bus[i]
        ti = t_bus[i]
        if fi < 1 || fi > max_bus_id || ti < 1 || ti > max_bus_id ||
           bus_to_zone[fi] == 0 || bus_to_zone[ti] == 0
            error(
                "Partitioning: branch $i references bus id outside bus table " *
                "(f=$fi, t=$ti, max_bus_id=$max_bus_id). Data inconsistent.",
            )
        end
        if bus_to_zone[fi] != bus_to_zone[ti]
            push!(interzonal_branches, i)
        end
    end

    return interzonal_branches
end

# ------------------------------------------------------------------------------
# Generic helpers
# ------------------------------------------------------------------------------

@inline function branch_in_service_mask(branch::BranchRaw{Float64})
    return branch.br_status .== 1
end

# Mark `zone_internal_membership[bus_id] = true` for all internal buses of
# `zone_id`. Assumes the caller `fill!(membership, false)` before calling.
@inline function _mark_internal_membership!(
    zone_internal_membership::Vector{Bool},
    full_bus_ids::AbstractVector{<:Integer},
    full_bus_zones::AbstractVector{<:Integer},
    zone_id::Int,
)
    @inbounds for i in eachindex(full_bus_ids)
        full_bus_zones[i] == zone_id || continue
        zone_internal_membership[full_bus_ids[i]] = true
    end
    return nothing
end

# Extend `zone_subprob_membership` (which starts as a copy of the internal
# membership) with external neighbors connected via any in-service branch.
@inline function _collect_external_neighbor_buses!(
    zone_subprob_membership::Vector{Bool},
    zone_internal_membership::Vector{Bool},
    full_f_bus_ids::AbstractVector{Int},
    full_t_bus_ids::AbstractVector{Int},
    full_branch_in_service_mask::AbstractVector{Bool},
)
    @inbounds for i in eachindex(full_f_bus_ids)
        full_branch_in_service_mask[i] || continue

        f_bus_id = full_f_bus_ids[i]
        t_bus_id = full_t_bus_ids[i]

        f_in_zone = zone_internal_membership[f_bus_id]
        t_in_zone = zone_internal_membership[t_bus_id]

        if f_in_zone && !t_in_zone
            zone_subprob_membership[t_bus_id] = true
        elseif t_in_zone && !f_in_zone
            zone_subprob_membership[f_bus_id] = true
        end
    end
    return nothing
end

# Collect global bus indices (positions in `full_bus_ids`) whose id is marked
# in `membership`. Output is naturally in ascending position order.
@inline function _collect_member_positions(
    full_bus_ids::AbstractVector{<:Integer},
    membership::Vector{Bool},
)
    out = Int[]
    sizehint!(out, count(membership))
    @inbounds for i in eachindex(full_bus_ids)
        membership[full_bus_ids[i]] && push!(out, i)
    end
    return out
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
        full_bus.pd[keep],
        full_bus.qd[keep],
        full_bus.gs[keep],
        full_bus.bs[keep],
        full_bus.area[keep],
        full_bus.vm[keep],
        full_bus.va[keep],
        full_bus.base_kv[keep],
        full_bus.zone[keep],
        full_bus.vmax[keep],
        full_bus.vmin[keep],
        full_bus.lam_p[keep],
        full_bus.lam_q[keep],
        full_bus.mu_vmax[keep],
        full_bus.mu_vmin[keep],
        full_bus.is_virtual[keep],
    )
end

@inline function _mark_virtual_buses!(
    bus::BusRaw{Float64},
    zone_internal_membership::Vector{Bool},
)
    @inbounds for i in eachindex(bus.bus_i)
        bus_id = bus.bus_i[i]
        is_virtual = !zone_internal_membership[bus_id]

        bus.is_virtual[i] = is_virtual

        if is_virtual
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
    zone_internal_membership::Vector{Bool},
)
    out = Int[]
    @inbounds for i in eachindex(component_bus_ids)
        b = Int(component_bus_ids[i])
        b >= 1 && b <= length(zone_internal_membership) &&
            zone_internal_membership[b] && push!(out, i)
    end
    return out
end

@inline function _slice_gen(full_gen::GenRaw{Float64}, keep::Vector{Int})
    return GenRaw{Float64}(
        full_gen.gen_bus[keep],
        full_gen.pg[keep],
        full_gen.qg[keep],
        full_gen.qmax[keep],
        full_gen.qmin[keep],
        full_gen.vg[keep],
        full_gen.mbase[keep],
        full_gen.gen_status[keep],
        full_gen.pmax[keep],
        full_gen.pmin[keep],
        full_gen.pc1[keep],
        full_gen.pc2[keep],
        full_gen.qc1min[keep],
        full_gen.qc1max[keep],
        full_gen.qc2min[keep],
        full_gen.qc2max[keep],
        full_gen.ramp_agc[keep],
        full_gen.ramp_10[keep],
        full_gen.ramp_30[keep],
        full_gen.ramp_q[keep],
        full_gen.apf[keep],
        full_gen.mu_pmax[keep],
        full_gen.mu_pmin[keep],
        full_gen.mu_qmax[keep],
        full_gen.mu_qmin[keep],
        full_gen.plant_type[keep],
        full_gen.curtailment_cost[keep],
    )
end

@inline function _slice_gencosts!(
    zone_case::CaseRaw,
    full::CaseRaw,
    zone_gen_idx::Vector{Int},
)
    gc = full.gencost

    zone_case.gencost = GencostRaw{Float64}(
        gc.model[zone_gen_idx],
        gc.startup[zone_gen_idx],
        gc.shutdown[zone_gen_idx],
        gc.ncost[zone_gen_idx],
        gc.cost_coefficients[zone_gen_idx],
    )

    return nothing
end

@inline function _filter_branches_incident_to_zone(
    full_f_bus_ids::AbstractVector{<:Integer},
    full_t_bus_ids::AbstractVector{<:Integer},
    full_branch_in_service_mask::AbstractVector{Bool},
    zone_internal_membership::Vector{Bool},
)
    out = Int[]
    sizehint!(out, length(full_f_bus_ids))
    @inbounds for i in eachindex(full_f_bus_ids)
        full_branch_in_service_mask[i] || continue
        f = Int(full_f_bus_ids[i])
        t = Int(full_t_bus_ids[i])
        if zone_internal_membership[f] || zone_internal_membership[t]
            push!(out, i)
        end
    end
    return out
end

@inline function _slice_branch(full_branch::BranchRaw{Float64}, keep::Vector{Int})
    return BranchRaw{Float64}(
        full_branch.f_bus[keep],
        full_branch.t_bus[keep],
        full_branch.br_r[keep],
        full_branch.br_x[keep],
        full_branch.br_b[keep],
        full_branch.rate_a[keep],
        full_branch.rate_b[keep],
        full_branch.rate_c[keep],
        full_branch.tap[keep],
        full_branch.shift[keep],
        full_branch.br_status[keep],
        full_branch.angmin[keep],
        full_branch.angmax[keep],
        full_branch.pf[keep],
        full_branch.qf[keep],
        full_branch.pt[keep],
        full_branch.qt[keep],
        full_branch.mu_sf[keep],
        full_branch.mu_st[keep],
        full_branch.mu_angmin[keep],
        full_branch.mu_angmax[keep],
    )
end

@inline function _slice_storage(full_storage::StorageRaw{Float64}, keep::Vector{Int})
    return StorageRaw{Float64}(
        full_storage.storage_bus[keep],
        full_storage.ps[keep],
        full_storage.qs[keep],
        full_storage.energy[keep],
        full_storage.energy_rating[keep],
        full_storage.charge_rating[keep],
        full_storage.discharge_rating[keep],
        full_storage.charge_efficiency[keep],
        full_storage.discharge_efficiency[keep],
        full_storage.thermal_rating[keep],
        full_storage.qmin[keep],
        full_storage.qmax[keep],
        full_storage.r[keep],
        full_storage.x[keep],
        full_storage.p_loss[keep],
        full_storage.q_loss[keep],
        full_storage.status[keep],
    )
end

@inline function _filter_timeseries!(
    zone_case::CaseRaw,
    full::CaseRaw,
    zone_internal_membership::Vector{Bool},
    zone_gen_idx::Vector{Int},
)
    if full.bus_timeseries !== nothing
        ts = full.bus_timeseries
        zone_ts_idx = _filter_by_zone_bus(ts.bus, zone_internal_membership)

        zone_case.bus_timeseries = BusTimeseriesRaw{Float64}(
            ts.time[zone_ts_idx],
            ts.bus[zone_ts_idx],
            ts.pd[zone_ts_idx],
            ts.qd[zone_ts_idx],
            ts.gs[zone_ts_idx],
            ts.bs[zone_ts_idx],
        )
    end

    if full.gencost_timeseries !== nothing
        ts = full.gencost_timeseries

        # Remap global gen ids to zone-local gen ids (1..length(zone_gen_idx)).
        # `zone_gen_idx[k]` is the GLOBAL gen id of local gen `k`.
        n_global_gen = length(full.gen.gen_bus)
        global_to_local_gen = zeros(Int, n_global_gen)
        @inbounds for (local_k, global_k) in enumerate(zone_gen_idx)
            global_to_local_gen[global_k] = local_k
        end

        zone_gc_ts_idx = Int[]
        sizehint!(zone_gc_ts_idx, length(ts.gen_id))
        @inbounds for i in eachindex(ts.gen_id)
            g = ts.gen_id[i]
            (g >= 1 && g <= n_global_gen && global_to_local_gen[g] != 0) || continue
            push!(zone_gc_ts_idx, i)
        end

        remapped_gen_id = Vector{Int}(undef, length(zone_gc_ts_idx))
        @inbounds for (k, idx) in enumerate(zone_gc_ts_idx)
            remapped_gen_id[k] = global_to_local_gen[ts.gen_id[idx]]
        end

        zone_case.gencost_timeseries = GencostTimeseriesRaw{Float64}(
            ts.time[zone_gc_ts_idx],
            remapped_gen_id,
            ts.model[zone_gc_ts_idx],
            ts.startup[zone_gc_ts_idx],
            ts.shutdown[zone_gc_ts_idx],
            ts.ncost[zone_gc_ts_idx],
            ts.cost_coefficients[zone_gc_ts_idx],
            ts.tou[zone_gc_ts_idx],
            ts.feed_in_tariff[zone_gc_ts_idx],
        )
    end

    return nothing
end
