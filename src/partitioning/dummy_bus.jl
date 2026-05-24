# ============================================================
# Sectioning dispatch: Dummy Bus
# ============================================================

function prepare_distributed_partitions(full::CaseRaw, ::DummyBus)
    cas = deepcopy(full)

    bus = cas.bus
    br = cas.branch
    gen = cas.gen
    gc = cas.gencost

    bus_id_to_idx = Dict{Int, Int}()
    for (i, bid) in enumerate(bus.bus_i)
        bus_id_to_idx[bid] = i
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

    isempty(tie_idxs) && return prepare_distributed_partitions(cas, TieLine())

    dummy_bus_ids = Set{Int}()
    next_bus_id = maximum(bus.bus_i) + 1
    half_a_globals = Int[]  # global idx of (from → dummy) branch
    half_b_globals = Int[]  # global idx of (dummy → to) branch

    for b in tie_idxs
        dummy_id = next_bus_id
        next_bus_id += 1

        fi = bus_id_to_idx[br.f_bus[b]]

        push!(bus.bus_i, dummy_id)
        push!(bus.bus_type, 1)
        push!(bus.pd, 0.0)
        push!(bus.qd, 0.0)
        push!(bus.gs, 0.0)
        push!(bus.bs, 0.0)
        push!(bus.area, bus.area[fi])
        push!(bus.vm, 1.0)
        push!(bus.va, 0.0)
        push!(bus.base_kv, bus.base_kv[fi])
        push!(bus.zone, bus.zone[fi])
        push!(bus.vmax, 1.1)
        push!(bus.vmin, 0.9)
        push!(bus.lam_p, 0.0)
        push!(bus.lam_q, 0.0)
        push!(bus.mu_vmax, 0.0)
        push!(bus.mu_vmin, 0.0)
        push!(bus.is_virtual, false)
        bus_id_to_idx[dummy_id] = length(bus.bus_i)
        push!(dummy_bus_ids, dummy_id)
        f_orig = br.f_bus[b]
        t_orig = br.t_bus[b]
        r_h = br.br_r[b] / 2.0
        x_h = br.br_x[b] / 2.0
        b_h = br.br_b[b] / 2.0

        has_xfmr = (br.tap[b] != 0.0) || (br.shift[b] != 0.0)

        _push_branch!(br, f_orig, dummy_id, r_h, x_h, b_h, br.rate_a[b],
            br.rate_b[b], br.rate_c[b],
            has_xfmr ? br.tap[b] : 0.0,
            has_xfmr ? br.shift[b] : 0.0,
            br.angmin[b], br.angmax[b])
        push!(half_a_globals, length(br.f_bus))

        _push_branch!(br, dummy_id, t_orig, r_h, x_h, b_h, br.rate_a[b],
            br.rate_b[b], br.rate_c[b], 0.0, 0.0, br.angmin[b], br.angmax[b])
        push!(half_b_globals, length(br.f_bus))

        br.br_status[b] = 0
        _push_dummy_gen!(gen, gc, dummy_id, cas.baseMVA)
    end

    full_bus_ids = cas.bus.bus_i
    full_bus_zones = cas.bus.zone
    full_f = cas.branch.f_bus
    full_t = cas.branch.t_bus
    full_mask = branch_in_service_mask(cas.branch)

    n_full_bus = length(full_bus_ids)
    n_full_branches = length(full_f)
    max_bus_id_val = maximum(full_bus_ids)

    zone_ids = sort!(unique(full_bus_zones))
    n_zones = length(zone_ids)
    parts = Vector{CaseRaw}(undef, n_zones)

    zone_internal = zeros(Bool, max_bus_id_val)
    zone_subprob = zeros(Bool, max_bus_id_val)
    bus_zone_count = zeros(Int, max_bus_id_val)

    @inbounds for (i, zid) in enumerate(zone_ids)
        fill!(zone_internal, false)
        _mark_internal_membership!(zone_internal, full_bus_ids, full_bus_zones, zid)

        copyto!(zone_subprob, zone_internal)

        for bi in eachindex(full_f)
            full_mask[bi] || continue
            fid = full_f[bi]
            tid = full_t[bi]
            f_in = zone_internal[fid]
            t_in = zone_internal[tid]

            if f_in && !t_in && (tid in dummy_bus_ids)
                zone_subprob[tid] = true
            elseif t_in && !f_in && (fid in dummy_bus_ids)
                zone_subprob[fid] = true
            end
        end

        bus_idx_full = _collect_member_positions(full_bus_ids, zone_subprob)

        zc = _build_zone_meta(cas, zid)
        zc.bus = _slice_bus(cas.bus, bus_idx_full)
        _mark_virtual_buses!(zc.bus, zone_internal)
        gen_idx = _filter_by_zone_bus(cas.gen.gen_bus, zone_subprob)
        zc.gen = _slice_gen(cas.gen, gen_idx)
        _slice_gencosts!(zc, cas, gen_idx)

        br_idx = Int[]
        for bi in eachindex(full_f)
            full_mask[bi] || continue
            zone_subprob[full_f[bi]] && zone_subprob[full_t[bi]] && push!(br_idx, bi)
        end
        zc.branch = _slice_branch(cas.branch, br_idx)

        if cas.storage !== nothing
            sidx = _filter_by_zone_bus(cas.storage.storage_bus, zone_internal)
            zc.storage = _slice_storage(cas.storage, sidx)
        end

        _filter_timeseries!(zc, cas, zone_internal, gen_idx)

        zc.dist.zone_branch_idx_in_full = br_idx
        zc.dist.branch_idx_full_to_zone = zeros(Int, n_full_branches)
        for (li, gi) in enumerate(br_idx)
            zc.dist.branch_idx_full_to_zone[gi] = li
        end

        add_case_index_maps!(zc)

        for b in zc.bus.bus_i
            bus_zone_count[b] += 1
        end

        parts[i] = zc
    end

    full_shared = Int[]
    @inbounds for b in 1:max_bus_id_val
        bus_zone_count[b] > 1 && push!(full_shared, b)
    end
    shared_mask = zeros(Bool, max_bus_id_val)
    for b in full_shared
        shared_mask[b] = true
    end
    cas.dist.full_shared_bus_ids = full_shared
    full.dist.full_shared_bus_ids = full_shared
    for zc in parts
        zc.dist.zone_shared_bus_ids = [b for b in zc.bus.bus_i if shared_mask[b]]
    end
    zone_id_to_pos = Dict{Int, Int}()
    for (pos, zid) in enumerate(zone_ids)
        zone_id_to_pos[zid] = pos
    end

    pairs_nt = NTuple{4, Int}[]
    for k in eachindex(tie_idxs)
        za_id, zb_id = tie_orig_zones[k]
        ga = half_a_globals[k]
        gb = half_b_globals[k]
        push!(pairs_nt,
            (ga, zone_id_to_pos[za_id], gb, zone_id_to_pos[zb_id]),
        )
    end

    full.dist.branch_pairs = pairs_nt
    cas.dist.branch_pairs = pairs_nt

    return (zone_ids = zone_ids, parts = parts)
end
