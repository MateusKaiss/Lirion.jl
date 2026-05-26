# ============================================================
# Sectioning dispatch: Split Bus
# ============================================================

function prepare_distributed_partitions(full::CaseRaw, ::SplitBus)
    full_bus_ids = full.bus.bus_i
    full_bus_zones = full.bus.zone
    full_f = full.branch.f_bus
    full_t = full.branch.t_bus
    in_service = branch_in_service_mask(full.branch)

    n_full_bus = length(full_bus_ids)
    n_full_branches = length(full_f)
    max_bus_id = n_full_bus == 0 ? 0 : maximum(full_bus_ids)
    bus_to_zone = zeros(Int, max_bus_id)
    @inbounds for i in 1:n_full_bus
        bus_to_zone[full_bus_ids[i]] = full_bus_zones[i]
    end

    zone_ids = sort!(unique(full_bus_zones))
    n_zones = length(zone_ids)
    zone_id_to_pos = Dict{Int, Int}()
    @inbounds for (pos, zid) in enumerate(zone_ids)
        zone_id_to_pos[zid] = pos
    end

    split_set = Set{Tuple{Int, Int}}()
    branch_owner = zeros(Int, n_full_branches)
    @inbounds for b in eachindex(full_f)
        in_service[b] || continue
        f = full_f[b]
        t = full_t[b]
        zf = bus_to_zone[f]
        zt = bus_to_zone[t]
        if zf == zt
            branch_owner[b] = zf
        else
            branch_owner[b] = zt
            push!(split_set, (f, zt))
        end
    end

    zone_ghost_buses = Dict{Int, Vector{Int}}()
    @inbounds for zid in zone_ids
        zone_ghost_buses[zid] = Int[]
    end
    @inbounds for (fbus, zb) in split_set
        push!(zone_ghost_buses[zb], fbus)
    end

    parts = Vector{CaseRaw}(undef, n_zones)
    bus_zone_count = zeros(Int, max_bus_id)
    zone_internal = zeros(Bool, max_bus_id)
    zone_subprob = zeros(Bool, max_bus_id)

    @inbounds for (i, zid) in enumerate(zone_ids)
        fill!(zone_internal, false)
        _mark_internal_membership!(zone_internal, full_bus_ids, full_bus_zones, zid)

        copyto!(zone_subprob, zone_internal)
        for fbus in zone_ghost_buses[zid]
            zone_subprob[fbus] = true
        end

        zone_bus_idx_full = _collect_member_positions(full_bus_ids, zone_subprob)

        zone_case = _build_zone_meta(full, zid)
        zone_case.bus = _slice_bus(full.bus, zone_bus_idx_full)
        _mark_virtual_buses!(zone_case.bus, zone_internal)
        for k in eachindex(zone_case.bus.bus_i)
            b = zone_case.bus.bus_i[k]
            if !zone_internal[b]
                zone_case.bus.is_virtual[k] = false
                zone_case.bus.pd[k] = 0.0
                zone_case.bus.qd[k] = 0.0
                zone_case.bus.gs[k] = 0.0
                zone_case.bus.bs[k] = 0.0
            end
        end

        zone_gen_idx = _filter_by_zone_bus(full.gen.gen_bus, zone_internal)
        zone_case.gen = _slice_gen(full.gen, zone_gen_idx)
        _slice_gencosts!(zone_case, full, zone_gen_idx)
        zone_branch_idx_in_full = Int[]
        sizehint!(zone_branch_idx_in_full, n_full_branches)
        for b in eachindex(full_f)
            in_service[b] || continue
            branch_owner[b] == zid || continue
            push!(zone_branch_idx_in_full, b)
        end
        zone_case.branch = _slice_branch(full.branch, zone_branch_idx_in_full)

        if full.storage !== nothing
            zone_storage_idx =
                _filter_by_zone_bus(full.storage.storage_bus, zone_internal)
            zone_case.storage = _slice_storage(full.storage, zone_storage_idx)
        end

        _filter_timeseries!(zone_case, full, zone_internal, zone_gen_idx)

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
    shared_mask = zeros(Bool, max_bus_id)
    for b in full_shared_bus_ids
        shared_mask[b] = true
    end
    full.dist.full_shared_bus_ids = full_shared_bus_ids

    for zc in parts
        zc.dist.zone_shared_bus_ids =
            [b for b in zc.bus.bus_i if shared_mask[b]]
    end
    neighbor_zones_of = Dict{Int, Vector{Int}}()
    for (fbus, zb) in split_set
        push!(get!(neighbor_zones_of, fbus, Int[]), zb)
    end

    # N-way SplitBus (a boundary bus shared by >2 zones) is allowed here. Models
    # not yet migrated to the bus_groups consensus guard against it at
    # workspace-build time (see `_assert_bilateral_split`).

    # Pairwise representation (one entry per (split bus, neighbor zone)).
    # Correct for bilateral boundaries; the not-yet-N-way models guard against
    # >2-zone groups at workspace-build time.
    bus_pairs = NTuple{3, Int}[]
    sizehint!(bus_pairs, length(split_set))
    for (fbus, zb) in split_set
        za = bus_to_zone[fbus]
        push!(
            bus_pairs,
            (fbus, zone_id_to_pos[za], zone_id_to_pos[zb]),
        )
    end
    full.dist.bus_pairs = bus_pairs

    # N-way representation: (bus_id, [owner_zone_pos, neighbor_zone_pos...]).
    # Generalizes the pairwise form to a bus shared by more than two zones
    # (net-injection-zero consensus over all copies, Σ p_z = 0).
    bus_groups = Tuple{Int, Vector{Int}}[]
    sizehint!(bus_groups, length(neighbor_zones_of))
    for (fbus, zbs) in neighbor_zones_of
        za = bus_to_zone[fbus]
        zones_pos = Int[zone_id_to_pos[za]]
        for zb in sort!(unique(zbs))
            push!(zones_pos, zone_id_to_pos[zb])
        end
        push!(bus_groups, (fbus, zones_pos))
    end
    sort!(bus_groups; by = g -> g[1])
    full.dist.bus_groups = bus_groups

    return (zone_ids = zone_ids, parts = parts)
end
