# ============================================================
# Sectioning dispatch: Dummy Bus
# ============================================================

function prepare_distributed_partitions(full::CaseRaw, ::DummyBus)
    _require_interzonal_branches(full, "DummyBus")

    cas = deepcopy(full)

    bus = cas.bus
    br = cas.branch
    gen = cas.gen
    gen_cost = cas.gencost

    bus_id_to_idx = Dict{Int, Int}()
    for (i, b_id) in enumerate(bus.bus_i)
        bus_id_to_idx[b_id] = i
    end

    tie_idxs = Int[]
    tie_orig_zones = Tuple{Int, Int}[]
    for b in eachindex(br.f_bus)
        br.br_status[b] == 1 || continue
        fi = bus_id_to_idx[br.f_bus[b]]
        ti = bus_id_to_idx[br.t_bus[b]]
        if bus.zone[fi] != bus.zone[ti]
            push!(tie_idxs, b)
            push!(tie_orig_zones, (bus.zone[fi], bus.zone[ti]))
        end
    end

    dummy_bus_ids = Set{Int}()
    next_bus_id = maximum(bus.bus_i) + 1
    half_a_globals = Int[]
    half_b_globals = Int[]

    for b in tie_idxs
        dummy_id = next_bus_id
        next_bus_id += 1
        fi = bus_id_to_idx[br.f_bus[b]]
        _push_dummy_bus!(bus, dummy_id, bus.area[fi], bus.base_kv[fi], bus.zone[fi])
        bus_id_to_idx[dummy_id] = length(bus.bus_i)
        push!(dummy_bus_ids, dummy_id)

        f_orig = br.f_bus[b]
        t_orig = br.t_bus[b]
        r_h = br.br_r[b] / 2.0
        x_h = br.br_x[b] / 2.0
        b_c = br.br_b[b]                       # total π line charging (p.u.)
        has_xfmr = (br.tap[b] != 0.0) || (br.shift[b] != 0.0)

        _push_branch!(br, f_orig, dummy_id, r_h, x_h, 0.0,
            br.rate_a[b], br.rate_b[b], br.rate_c[b],
            has_xfmr ? br.tap[b] : 0.0,
            has_xfmr ? br.shift[b] : 0.0,
            br.angmin[b], br.angmax[b])
        push!(half_a_globals, length(br.f_bus))

        _push_branch!(br, dummy_id, t_orig, r_h, x_h, 0.0,
            br.rate_a[b], br.rate_b[b], br.rate_c[b],
            0.0, 0.0,
            br.angmin[b], br.angmax[b])
        push!(half_b_globals, length(br.f_bus))

        tap_mag = br.tap[b] == 0.0 ? 1.0 : br.tap[b]
        bus.bs[fi] += (b_c / 2.0) / tap_mag^2    
        bus.bs[bus_id_to_idx[t_orig]] += b_c / 2.0 

        br.br_status[b] = 0
        _push_dummy_gen!(gen, gen_cost, dummy_id, cas.baseMVA)
    end

    full_bus_ids = cas.bus.bus_i
    full_bus_zones = cas.bus.zone
    full_f_bus_ids = cas.branch.f_bus
    full_t_bus_ids = cas.branch.t_bus
    full_branch_in_service_mask = branch_in_service_mask(cas.branch)

    n_full_bus = length(full_bus_ids)
    n_full_branches = length(full_f_bus_ids)
    max_bus_id = maximum(full_bus_ids)

    zone_ids = sort!(unique(full_bus_zones))
    n_zones = length(zone_ids)
    parts = Vector{CaseRaw}(undef, n_zones)

    zone_internal_membership = zeros(Bool, max_bus_id)
    zone_subprob_membership = zeros(Bool, max_bus_id)
    bus_zone_count = zeros(Int, max_bus_id)

    @inbounds for (i, zone_id) in enumerate(zone_ids)
        fill!(zone_internal_membership, false)
        _mark_internal_membership!(zone_internal_membership, full_bus_ids, full_bus_zones, zone_id)

        copyto!(zone_subprob_membership, zone_internal_membership)
        _collect_external_dummy_neighbors!(
            zone_subprob_membership,
            zone_internal_membership,
            full_f_bus_ids,
            full_t_bus_ids,
            full_branch_in_service_mask,
            dummy_bus_ids,
        )

        zone_subprob_bus_idx_full = _collect_member_positions(full_bus_ids, zone_subprob_membership)

        zone_case = _build_zone_meta(cas, zone_id)
        zone_case.bus = _slice_bus(cas.bus, zone_subprob_bus_idx_full)
        _mark_virtual_buses!(zone_case.bus, zone_internal_membership)
        
        zone_gen_idx = _filter_by_zone_bus(cas.gen.gen_bus, zone_internal_membership)
        zone_case.gen = _slice_gen(cas.gen, zone_gen_idx)
        _slice_gencosts!(zone_case, cas, zone_gen_idx)

        zone_branch_idx_in_full = _filter_branches_inside_subprob(
            full_f_bus_ids,
            full_t_bus_ids,
            full_branch_in_service_mask,
            zone_subprob_membership,
        )
        zone_case.branch = _slice_branch(cas.branch, zone_branch_idx_in_full)

        if cas.storage !== nothing
            zone_storage_idx = _filter_by_zone_bus(cas.storage.storage_bus, zone_internal_membership)
            zone_case.storage = _slice_storage(cas.storage, zone_storage_idx)
        end

        _filter_timeseries!(zone_case, cas, zone_internal_membership, zone_gen_idx)

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
    cas.dist.full_shared_bus_ids = full_shared_bus_ids
    full.dist.full_shared_bus_ids = full_shared_bus_ids
    for zone_case in parts
        zone_case.dist.zone_shared_bus_ids = [b for b in zone_case.bus.bus_i if full_shared_membership[b]]
    end
    zone_id_to_pos = _compute_zone_id_to_pos(zone_ids)

    branch_pairs = NTuple{4, Int}[]
    for k in eachindex(tie_idxs)
        za_id, zb_id = tie_orig_zones[k]
        half_a_global = half_a_globals[k]
        half_b_global = half_b_globals[k]
        push!(branch_pairs,
            (half_a_global, zone_id_to_pos[za_id], half_b_global, zone_id_to_pos[zb_id]),
        )
    end

    full.dist.branch_pairs = branch_pairs
    cas.dist.branch_pairs = branch_pairs

    return (zone_ids = zone_ids, parts = parts)
end
