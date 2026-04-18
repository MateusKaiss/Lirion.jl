# ==============================================================================
# Sectioning dispatch: Split Bus (Method 3 — Harris et al. 2022, Fig. 1c)
# ------------------------------------------------------------------------------
# For each tie-line (f → t) crossing zones A (= zone(f)) and B (= zone(t)):
#   1. The tie-line is OWNED by zone B as an internal branch.
#   2. Zone A keeps bus `f` with all its physical content (load, gen, shunt).
#   3. Zone B gets a GHOST COPY of bus `f` (is_virtual=false, zeroed physical).
#   4. Shared variables at bus `f`: θ (both zones) and p_inj (virtual injection,
#      with sign-flipped consensus p_A + p_B = 0).
#
# Differs from TieLine (Method 1): no boundary ghost of `t` in zone A — zone A
# doesn't see bus `t` at all. Differs from DummyBus (Method 2): no new bus, no
# artificial generator, tie-line is NOT cut in half.
#
# Implementation note: the p_inj variables are added later at the ADMM layer
# (see `init_hybrid_partition!` for DC). This file only prepares the partition
# geometry and records `bus_pairs` for the consensus bookkeeping.
# ==============================================================================

function prepare_distributed_partitions(full::CaseRaw, ::SplitBus)
    full_bus_ids = full.bus.bus_i
    full_bus_zones = full.bus.zone
    full_f = full.branch.f_bus
    full_t = full.branch.t_bus
    in_service = branch_in_service_mask(full.branch)

    n_full_bus = length(full_bus_ids)
    n_full_branches = length(full_f)
    max_bus_id = n_full_bus == 0 ? 0 : maximum(full_bus_ids)

    # Dense bus_id -> zone_id lookup (0 if bus absent).
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

    # Identify split pairs (f_bus, zone_B) — each unique (from-bus, neighbor-zone)
    # combination across all tie-lines produces exactly one ghost copy in B.
    split_set = Set{Tuple{Int, Int}}()   # (from_bus_id, zone_B_id)
    # Per-branch owner zone: non-tie-lines are owned by their (single) zone;
    # tie-lines (f→t) crossing zones are owned by zone(t) (the `to` side).
    # This is the key piece that prevents a tie-line from being double-counted
    # when both endpoints end up in a zone's subproblem as native/ghost.
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

    # Pre-compute, per zone, which external ghost buses (split points) this zone
    # hosts. Indexed by zone_id.
    zone_ghost_buses = Dict{Int, Vector{Int}}()
    @inbounds for zid in zone_ids
        zone_ghost_buses[zid] = Int[]
    end
    @inbounds for (fbus, zb) in split_set
        push!(zone_ghost_buses[zb], fbus)
    end

    # Tie-lines belong to zone B = zone(t). From zone A's perspective they're
    # "outbound" and are skipped entirely (no ghost of t in A).
    #
    # For each zone, include internal branches (both endpoints internal) PLUS
    # inbound tie-lines (f external ghost, t internal).

    parts = Vector{CaseRaw}(undef, n_zones)
    bus_zone_count = zeros(Int, max_bus_id)

    # Scratch buffers reused across zones.
    zone_internal = zeros(Bool, max_bus_id)
    zone_subprob = zeros(Bool, max_bus_id)

    @inbounds for (i, zid) in enumerate(zone_ids)
        fill!(zone_internal, false)
        _mark_internal_membership!(zone_internal, full_bus_ids, full_bus_zones, zid)

        # Subproblem membership: internal buses + ghost copies of split points
        # whose tie-lines terminate in this zone.
        copyto!(zone_subprob, zone_internal)
        for fbus in zone_ghost_buses[zid]
            zone_subprob[fbus] = true
        end

        zone_bus_idx_full = _collect_member_positions(full_bus_ids, zone_subprob)

        zone_case = _build_zone_meta(full, zid)
        zone_case.bus = _slice_bus(full.bus, zone_bus_idx_full)
        _mark_virtual_buses!(zone_case.bus, zone_internal)
        # Ghost copies of split buses must participate in power balance (they
        # carry the p_inj slack). Override: is_virtual = false, physical zeroed.
        for k in eachindex(zone_case.bus.bus_i)
            b = zone_case.bus.bus_i[k]
            if !zone_internal[b]   # external ⇒ this is a split-bus ghost
                zone_case.bus.is_virtual[k] = false
                zone_case.bus.pd[k] = 0.0
                zone_case.bus.qd[k] = 0.0
                zone_case.bus.gs[k] = 0.0
                zone_case.bus.bs[k] = 0.0
            end
        end

        # Generators: only from internal buses (ghost copies have no generation).
        zone_gen_idx = _filter_by_zone_bus(full.gen.gen_bus, zone_internal)
        zone_case.gen = _slice_gen(full.gen, zone_gen_idx)
        _slice_gencosts!(zone_case, full, zone_gen_idx)

        # Branches: include only those OWNED by this zone (by `branch_owner`).
        # Pure internal branches are owned by their single zone; split tie-lines
        # are owned by the `to` zone (the one hosting the ghost of `from`).
        # Checking ownership explicitly avoids the double-counting trap where a
        # tie-line's endpoints both appear in this zone's subproblem (the `from`
        # as native, the `to` as ghost brought in by a DIFFERENT split bus).
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

    # Shared buses = appearing in ≥ 2 zones = exactly the split buses.
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

    # Build bus_pairs: for each (fbus, zone_B_id) entry in split_set,
    # emit (fbus, zone_A_pos, zone_B_pos). zone_A is the zone natively
    # owning fbus (from `bus_to_zone`).
    # Requirement: each split bus must be shared by exactly 2 zones.
    neighbor_zones_of = Dict{Int, Vector{Int}}()
    for (fbus, zb) in split_set
        push!(get!(neighbor_zones_of, fbus, Int[]), zb)
    end
    for (fbus, zs) in neighbor_zones_of
        if length(unique(zs)) > 1
            error(
                "SplitBus: bus $fbus is a split point for multiple neighbor " *
                "zones $(unique(zs)); multi-zone consensus is not supported. " *
                "Use a partitioning where each boundary bus borders a single " *
                "neighbor zone.",
            )
        end
    end

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

    return (zone_ids = zone_ids, parts = parts)
end
