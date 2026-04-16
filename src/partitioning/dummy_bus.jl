# ==============================================================================
# Sectioning dispatch: Dummy Bus (Method 2 — Harris et al. 2022, Fig. 1b)
# ------------------------------------------------------------------------------
# For each tie-line (i→j) between zones A and B:
#   1. Remove the original tie-line
#   2. Insert a dummy bus d
#   3. Create branch (i→d) in zone A and branch (d→j) in zone B
#   4. Add a zero-cost artificial generator at d
#   5. d appears in BOTH zones as a shared bus
#
# The ONLY connection between zones is through the dummy bus — no ghost copies
# of original buses are created (unlike TieLine / Method 1).
# ==============================================================================

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

    # Identify tie-lines (active branches crossing zone boundaries).
    tie_idxs = Int[]
    for b in eachindex(br.f_bus)
        br.br_status[b] == 1 || continue
        fi = bus_id_to_idx[br.f_bus[b]]
        ti = bus_id_to_idx[br.t_bus[b]]
        bus.zone[fi] != bus.zone[ti] && push!(tie_idxs, b)
    end

    isempty(tie_idxs) && return prepare_distributed_partitions(cas, TieLine())

    dummy_bus_ids = Set{Int}()
    next_bus_id = maximum(bus.bus_i) + 1

    for b in tie_idxs
        dummy_id = next_bus_id
        next_bus_id += 1

        fi = bus_id_to_idx[br.f_bus[b]]

        # ── 1. Create the dummy bus (assigned to from-side zone) ────────
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

        # ── 2. Split branch into two halves ─────────────────────────────
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

        _push_branch!(br, dummy_id, t_orig, r_h, x_h, b_h, br.rate_a[b],
            br.rate_b[b], br.rate_c[b], 0.0, 0.0, br.angmin[b], br.angmax[b])

        # ── 3. Deactivate original tie-line ─────────────────────────────
        br.br_status[b] = 0

        # ── 4. Artificial generator at dummy bus ────────────────────────
        _push_dummy_gen!(gen, gc, dummy_id, cas.baseMVA)
    end

    # ====================================================================
    # Custom partitioning: each zone gets its internal buses + dummy buses
    # that have branches TO this zone. No ghost buses of original endpoints.
    # ====================================================================
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

        # Add ONLY dummy buses from other zones that have branches to this zone.
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

        # Generators: from internal buses + dummy buses in subproblem
        gen_idx = _filter_by_zone_bus(cas.gen.gen_bus, zone_subprob)
        zc.gen = _slice_gen(cas.gen, gen_idx)
        _slice_gencosts!(zc, cas, gen_idx)

        # Branches: both endpoints must be in the subproblem
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

    # Shared buses = those appearing in ≥2 zones (only dummy buses here).
    full_shared = Int[]
    @inbounds for b in 1:max_bus_id_val
        bus_zone_count[b] > 1 && push!(full_shared, b)
    end
    shared_mask = zeros(Bool, max_bus_id_val)
    for b in full_shared
        shared_mask[b] = true
    end
    cas.dist.full_shared_bus_ids = full_shared
    for zc in parts
        zc.dist.zone_shared_bus_ids = [b for b in zc.bus.bus_i if shared_mask[b]]
    end

    return (zone_ids = zone_ids, parts = parts)
end

# ── Helpers ──────────────────────────────────────────────────────────────

function _push_branch!(br, f, t, r, x, b_sh, ra, rb, rc, tap, shift, amin, amax)
    push!(br.f_bus, f)
    push!(br.t_bus, t)
    push!(br.br_r, r)
    push!(br.br_x, x)
    push!(br.br_b, b_sh)
    push!(br.rate_a, ra)
    push!(br.rate_b, rb)
    push!(br.rate_c, rc)
    push!(br.tap, tap)
    push!(br.shift, shift)
    push!(br.br_status, 1)
    push!(br.angmin, amin)
    push!(br.angmax, amax)
    push!(br.pf, 0.0)
    push!(br.qf, 0.0)
    push!(br.pt, 0.0)
    push!(br.qt, 0.0)
    push!(br.mu_sf, 0.0)
    push!(br.mu_st, 0.0)
    push!(br.mu_angmin, 0.0)
    push!(br.mu_angmax, 0.0)
    return nothing
end

function _push_dummy_gen!(gen, gc, bus_id::Int, baseMVA::Float64)
    push!(gen.gen_bus, bus_id)
    push!(gen.pg, 0.0)
    push!(gen.qg, 0.0)
    push!(gen.qmax, 1e4)
    push!(gen.qmin, -1e4)
    push!(gen.vg, 1.0)
    push!(gen.mbase, baseMVA)
    push!(gen.gen_status, 1)
    push!(gen.pmax, 1e4)
    push!(gen.pmin, -1e4)
    push!(gen.pc1, 0.0)
    push!(gen.pc2, 0.0)
    push!(gen.qc1min, 0.0)
    push!(gen.qc1max, 0.0)
    push!(gen.qc2min, 0.0)
    push!(gen.qc2max, 0.0)
    push!(gen.ramp_agc, 0.0)
    push!(gen.ramp_10, 0.0)
    push!(gen.ramp_30, 0.0)
    push!(gen.ramp_q, 0.0)
    push!(gen.apf, 0.0)
    push!(gen.mu_pmax, 0.0)
    push!(gen.mu_pmin, 0.0)
    push!(gen.mu_qmax, 0.0)
    push!(gen.mu_qmin, 0.0)
    push!(gen.plant_type, "artificial")
    push!(gen.curtailment_cost, 0.0)

    push!(gc.model, 2)
    push!(gc.startup, 0.0)
    push!(gc.shutdown, 0.0)
    push!(gc.ncost, 2)
    push!(gc.cost_coefficients, [0.0, 0.0])
    return nothing
end
