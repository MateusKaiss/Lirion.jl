# ============================================================
# MetisExt — graph-based auto-partitioning for distributed OPF
# ============================================================

module MetisExt

using Lirion
import Metis
using Graphs
using SimpleWeightedGraphs

# ============================================================
# Graph construction from CaseRaw
# ============================================================

function _build_impedance_graph(case::Lirion.CaseRaw)
    bus_ids = case.bus.bus_i
    n = length(bus_ids)

    id_to_idx = Dict{Int, Int}()
    sizehint!(id_to_idx, n)
    for i in 1:n
        id_to_idx[bus_ids[i]] = i
    end

    g = SimpleWeightedGraph(n)
    br = case.branch
    for i in eachindex(br.f_bus)
        br.br_status[i] == 1 || continue
        fi = get(id_to_idx, br.f_bus[i], 0)
        ti = get(id_to_idx, br.t_bus[i], 0)
        (fi == 0 || ti == 0 || fi == ti) && continue

        x = abs(br.br_x[i])
        w = x > 1e-8 ? 1.0 / x : 1e8

        tap = isnan(br.tap[i]) ? 0.0 : br.tap[i]
        shift = isnan(br.shift[i]) ? 0.0 : br.shift[i]
        if (tap != 0.0 && tap != 1.0) || shift != 0.0
            w *= 5.0
        end

        rate = br.rate_a[i]
        if rate > 0.0
            w *= 1.0 + 0.002 * rate
        end

        existing = SimpleWeightedGraphs.get_weight(g, fi, ti)
        if existing > 0
            w = max(w, existing)
        end
        add_edge!(g, fi, ti, w)
    end

    return g, id_to_idx
end

function _to_metis_weights(g::SimpleWeightedGraph)
    wts = Int[]
    for e in edges(g)
        w = weight(e)
        push!(wts, clamp(round(Int, w * 100), 1, 50_000))
    end
    return wts
end

# ============================================================
# METIS k-way partition
# ============================================================

function _metis_partition(g::SimpleWeightedGraph, k::Int)
    # Metis.jl requires integer edge weights, so rebuild the graph with
    # rounded / clamped integer weights (preserving topology and relative
    # magnitudes via the same scaling used by `_to_metis_weights`).
    n = Graphs.nv(g)
    g_int = SimpleWeightedGraph(n)
    for e in edges(g)
        w_float = weight(e)
        w_int = clamp(round(Int, w_float * 100), 1, 50_000)
        add_edge!(g_int, src(e), dst(e), w_int)
    end
    mg = Metis.graph(g_int; weights = true)
    parts = Metis.partition(mg, k; alg = :KWAY)
    return Int.(parts)   # Metis returns Int32; downstream helpers assume Int
end

# ============================================================
# Post-processing: connectivity guarantee (BFS merge)
# ============================================================

function _fix_connectivity!(zone::Vector{Int}, adj::Vector{Vector{Int}})
    n = length(zone)
    k = maximum(zone)
    visited = falses(n)

    for z in 1:k
        fill!(visited, false)
        members = findall(==(z), zone)
        isempty(members) && continue

        queue = [members[1]]
        visited[members[1]] = true
        head = 1
        while head <= length(queue)
            u = queue[head]
            head += 1
            for v in adj[u]
                if !visited[v] && zone[v] == z
                    visited[v] = true
                    push!(queue, v)
                end
            end
        end

        for m in members
            visited[m] && continue
            best_z = z
            best_count = 0
            for v in adj[m]
                zv = zone[v]
                zv == z && continue
                c = count(==(zv), @view(zone[adj[m]]))
                if c > best_count
                    best_count = c
                    best_z = zv
                end
            end
            zone[m] = best_z
        end
    end
    return nothing
end

# ============================================================
# Post-processing: ensure each zone has ≥1 generator (Pmax>0)
# ============================================================

function _ensure_generation!(case::Lirion.CaseRaw, zone::Vector{Int},
    id_to_idx::Dict{Int, Int})
    k = maximum(zone)
    gen = case.gen
    gen_per_zone = zeros(Int, k)
    pmax_per_zone = zeros(Float64, k)

    for i in eachindex(gen.gen_bus)
        gen.gen_status[i] == 1 || continue
        gen.pmax[i] > 0 || continue
        idx = get(id_to_idx, gen.gen_bus[i], 0)
        idx == 0 && continue
        z = zone[idx]
        gen_per_zone[z] += 1
        pmax_per_zone[z] += gen.pmax[i]
    end

    n = length(zone)
    br = case.branch
    for z in 1:k
        gen_per_zone[z] > 0 && continue

        best_bus = 0
        best_dist = typemax(Int)

        for i in eachindex(gen.gen_bus)
            gen.gen_status[i] == 1 || continue
            gen.pmax[i] > 0 || continue
            gi = get(id_to_idx, gen.gen_bus[i], 0)
            gi == 0 && continue

            for j in 1:n
                zone[j] == z || continue
                for b in eachindex(br.f_bus)
                    br.br_status[b] == 1 || continue
                    fi = get(id_to_idx, br.f_bus[b], 0)
                    ti = get(id_to_idx, br.t_bus[b], 0)
                    if (fi == gi && zone[ti] == z) || (ti == gi && zone[fi] == z)
                        if best_bus == 0
                            best_bus = gi
                        end
                    end
                end
            end
        end

        if best_bus > 0
            zone[best_bus] = z
        end
    end
    return nothing
end

# ============================================================
# Refinement: greedy boundary-bus relocation to reduce coupling
# ============================================================

function _count_shared_buses(zone::Vector{Int}, adj::Vector{Vector{Int}})
    count = 0
    for i in eachindex(zone)
        z = zone[i]
        for v in adj[i]
            if zone[v] != z
                count += 1
                break
            end
        end
    end
    return count
end

function _count_tie_lines(zone::Vector{Int}, adj::Vector{Vector{Int}})
    count = 0
    for i in eachindex(zone)
        for v in adj[i]
            if v > i && zone[v] != zone[i]
                count += 1
            end
        end
    end
    return count
end

function _zone_connected_after_remove(zone::Vector{Int}, adj::Vector{Vector{Int}},
    node::Int, z::Int)
    members = [m for m in eachindex(zone) if zone[m] == z && m != node]
    isempty(members) && return true
    length(members) == 1 && return true

    visited = Set{Int}()
    push!(visited, members[1])
    queue = [members[1]]
    head = 1
    while head <= length(queue)
        u = queue[head]
        head += 1
        for v in adj[u]
            if zone[v] == z && v != node && !(v in visited)
                push!(visited, v)
                push!(queue, v)
            end
        end
    end
    return length(visited) == length(members)
end

function _has_gen_after_remove(case::Lirion.CaseRaw, zone::Vector{Int},
    id_to_idx::Dict{Int, Int}, node::Int, z::Int)
    gen = case.gen
    for i in eachindex(gen.gen_bus)
        gen.gen_status[i] == 1 || continue
        gen.pmax[i] > 0 || continue
        gi = get(id_to_idx, gen.gen_bus[i], 0)
        (gi == 0 || gi == node) && continue
        zone[gi] == z && return true
    end
    return false
end

function _build_adjacency(case::Lirion.CaseRaw, id_to_idx::Dict{Int, Int})
    n = length(case.bus.bus_i)
    adj = [Int[] for _ in 1:n]
    br = case.branch
    for i in eachindex(br.f_bus)
        br.br_status[i] == 1 || continue
        fi = get(id_to_idx, br.f_bus[i], 0)
        ti = get(id_to_idx, br.t_bus[i], 0)
        (fi == 0 || ti == 0 || fi == ti) && continue
        push!(adj[fi], ti)
        push!(adj[ti], fi)
    end
    for a in adj
        sort!(unique!(a))
    end
    return adj
end

function _refine_partition!(case::Lirion.CaseRaw, zone::Vector{Int},
    adj::Vector{Vector{Int}}, id_to_idx::Dict{Int, Int};
    max_passes::Int = 50)
    improved = true
    pass = 0
    while improved && pass < max_passes
        improved = false
        pass += 1
        baseline_shared = _count_shared_buses(zone, adj)
        baseline_ties = _count_tie_lines(zone, adj)

        for i in eachindex(zone)
            z_orig = zone[i]
            any(v -> zone[v] != z_orig, adj[i]) || continue

            neighbor_zones = Set{Int}()
            for v in adj[i]
                zone[v] != z_orig && push!(neighbor_zones, zone[v])
            end

            best_z = z_orig
            best_score = baseline_shared + 0.5 * baseline_ties

            for zn in neighbor_zones
                !_zone_connected_after_remove(zone, adj, i, z_orig) && continue
                !_has_gen_after_remove(case, zone, id_to_idx, i, z_orig) && continue

                zone[i] = zn
                new_shared = _count_shared_buses(zone, adj)
                new_ties = _count_tie_lines(zone, adj)
                new_score = new_shared + 0.5 * new_ties

                if new_score < best_score
                    best_score = new_score
                    best_z = zn
                end
                zone[i] = z_orig
            end

            if best_z != z_orig
                zone[i] = best_z
                improved = true
            end
        end
    end
    return nothing
end

# ============================================================
# Auto-choose k
# ============================================================

function _auto_k(n_bus::Int; kmin::Int = 2, kmax::Int = 8)
    k = max(kmin, round(Int, n_bus / 40))
    return clamp(k, kmin, kmax)
end

# ============================================================
# Remap zone labels to 1:k contiguous
# ============================================================

function _remap_zones!(zone::Vector{Int})
    labels = sort!(unique(zone))
    mapping = Dict{Int, Int}()
    for (i, z) in enumerate(labels)
        mapping[z] = i
    end
    for i in eachindex(zone)
        zone[i] = mapping[zone[i]]
    end
    return nothing
end

# ============================================================
# Main entry point (overrides the stub in Lirion)
# ============================================================

function _do_partition!(
    case::Lirion.CaseRaw;
    k::Union{Int, Nothing},
    refine::Bool,
    max_refine_passes::Int,
    seed::Int,
)
    n_bus = length(case.bus.bus_i)
    kval = k === nothing ? _auto_k(n_bus) : k
    kval = clamp(kval, 2, n_bus)

    Random.seed!(seed)

    g, id_to_idx = _build_impedance_graph(case)
    zone = _metis_partition(g, kval)
    adj = _build_adjacency(case, id_to_idx)

    _fix_connectivity!(zone, adj)
    _ensure_generation!(case, zone, id_to_idx)
    _remap_zones!(zone)

    if refine
        _refine_partition!(case, zone, adj, id_to_idx;
            max_passes = max_refine_passes)
        _fix_connectivity!(zone, adj)
        _remap_zones!(zone)
    end

    for i in eachindex(case.bus.bus_i)
        idx = id_to_idx[case.bus.bus_i[i]]
        case.bus.zone[i] = zone[idx]
    end

    k_actual = maximum(zone)
    n_shared = _count_shared_buses(zone, adj)
    n_ties = _count_tie_lines(zone, adj)
    @info "Partition created" zones=k_actual buses=n_bus shared_buses=n_shared tie_lines=n_ties

    return case
end

function Lirion.create_partition(
    case::Lirion.CaseRaw;
    k::Union{Int, Nothing} = nothing,
    refine::Bool = true,
    max_refine_passes::Int = 50,
    seed::Int = 42,
)
    return _do_partition!(case; k, refine, max_refine_passes, seed)
end

function Lirion.create_partition(
    path::AbstractString;
    k::Union{Int, Nothing} = nothing,
    save_path::Union{AbstractString, Nothing} = nothing,
    refine::Bool = true,
    max_refine_passes::Int = 50,
    seed::Int = 42,
)
    case = Lirion.load_input_data(path)
    _do_partition!(case; k, refine, max_refine_passes, seed)

    k_actual = maximum(case.bus.zone)
    if save_path === nothing
        base = replace(basename(path), r"\.m$" => "")
        dir = dirname(abspath(path))
        out = joinpath(dir, "$(base)_$(k_actual)zones.m")
    else
        out = String(save_path)
    end

    Lirion.write_matpower(case, out)
    return out
end

# ============================================================
# SplitBus-compatible partition (Method 3 bilateral constraint)
# ------------------------------------------------------------
# Guarantees: every bus has neighbors in AT MOST ONE foreign zone. Equivalently,
# no bus sits at the boundary of 3+ zones. Required by `Lirion.SplitBus()`.
# ============================================================

# Count how many buses violate the bilateral constraint (≥ 2 foreign zones).
function _count_split_bus_violations(zone::Vector{Int}, adj::Vector{Vector{Int}})
    violations = 0
    for i in eachindex(zone)
        z = zone[i]
        foreign = Set{Int}()
        for v in adj[i]
            zv = zone[v]
            zv != z && push!(foreign, zv)
            length(foreign) >= 2 && break
        end
        length(foreign) >= 2 && (violations += 1)
    end
    return violations
end

# Iteratively move violating buses to their dominant foreign zone until no
# violations remain (or max_passes exhausted).
function _enforce_split_bus_compatible!(
    case::Lirion.CaseRaw,
    zone::Vector{Int},
    adj::Vector{Vector{Int}},
    id_to_idx::Dict{Int, Int};
    max_passes::Int = 100,
)
    n = length(zone)
    foreign_set = Set{Int}()
    zone_count = Dict{Int, Int}()

    for pass in 1:max_passes
        changed = false
        for i in 1:n
            empty!(foreign_set)
            z = zone[i]
            for v in adj[i]
                zv = zone[v]
                zv != z && push!(foreign_set, zv)
            end
            length(foreign_set) >= 2 || continue

            # Pick dominant foreign zone (most neighbors of i located there).
            empty!(zone_count)
            for v in adj[i]
                zv = zone[v]
                zv != z && (zone_count[zv] = get(zone_count, zv, 0) + 1)
            end
            best_zone = z
            best_count = 0
            for (fz, c) in zone_count
                if c > best_count
                    best_count = c
                    best_zone = fz
                end
            end
            best_zone == z && continue

            # Safety: do not orphan a zone's generator and do not disconnect it.
            _zone_connected_after_remove(zone, adj, i, z) || continue
            _has_gen_after_remove(case, zone, id_to_idx, i, z) || continue

            zone[i] = best_zone
            changed = true
        end
        changed || break
    end

    # Final validation: raise an explicit error if any violation remains.
    for i in 1:n
        empty!(foreign_set)
        z = zone[i]
        for v in adj[i]
            zv = zone[v]
            zv != z && push!(foreign_set, zv)
        end
        if length(foreign_set) >= 2
            bus_id = case.bus.bus_i[i]
            fz_sorted = sort(collect(foreign_set))
            error(
                "create_partition_split_bus: no SplitBus-compatible partition " *
                "found at this k. Bus $bus_id still borders multiple foreign " *
                "zones $fz_sorted after $max_passes repair passes. Try a " *
                "smaller k (k=2 is always feasible).",
            )
        end
    end
    return nothing
end

function _do_partition_split_bus!(
    case::Lirion.CaseRaw;
    k::Union{Int, Nothing},
    refine::Bool,
    max_refine_passes::Int,
    max_repair_passes::Int,
    seed::Int,
)
    n_bus = length(case.bus.bus_i)
    kval = k === nothing ? _auto_k(n_bus) : k
    kval = clamp(kval, 2, n_bus)

    Random.seed!(seed)

    g, id_to_idx = _build_impedance_graph(case)
    zone = _metis_partition(g, kval)
    adj = _build_adjacency(case, id_to_idx)

    _fix_connectivity!(zone, adj)
    _ensure_generation!(case, zone, id_to_idx)
    _remap_zones!(zone)

    if refine
        _refine_partition!(case, zone, adj, id_to_idx;
            max_passes = max_refine_passes)
        _fix_connectivity!(zone, adj)
        _remap_zones!(zone)
    end

    # SplitBus-specific repair: moves buses until bilateral constraint holds.
    _enforce_split_bus_compatible!(
        case, zone, adj, id_to_idx;
        max_passes = max_repair_passes,
    )
    _fix_connectivity!(zone, adj)
    _remap_zones!(zone)

    # Write zones back to the bus table.
    for i in eachindex(case.bus.bus_i)
        idx = id_to_idx[case.bus.bus_i[i]]
        case.bus.zone[i] = zone[idx]
    end

    k_actual = maximum(zone)
    n_shared = _count_shared_buses(zone, adj)
    n_ties = _count_tie_lines(zone, adj)
    @info "Partition created (SplitBus-compatible)" zones = k_actual buses = n_bus shared_buses = n_shared tie_lines = n_ties

    return case
end

function Lirion.create_partition_split_bus(
    case::Lirion.CaseRaw;
    k::Union{Int, Nothing} = 2,
    refine::Bool = true,
    max_refine_passes::Int = 50,
    max_repair_passes::Int = 100,
    seed::Int = 42,
)
    return _do_partition_split_bus!(
        case;
        k = k,
        refine = refine,
        max_refine_passes = max_refine_passes,
        max_repair_passes = max_repair_passes,
        seed = seed,
    )
end

function Lirion.create_partition_split_bus(
    path::AbstractString;
    k::Union{Int, Nothing} = 2,
    save_path::Union{AbstractString, Nothing} = nothing,
    refine::Bool = true,
    max_refine_passes::Int = 50,
    max_repair_passes::Int = 100,
    seed::Int = 42,
)
    case = Lirion.load_input_data(path)
    _do_partition_split_bus!(
        case;
        k = k,
        refine = refine,
        max_refine_passes = max_refine_passes,
        max_repair_passes = max_repair_passes,
        seed = seed,
    )

    k_actual = maximum(case.bus.zone)
    if save_path === nothing
        base = replace(basename(path), r"\.m$" => "")
        dir = dirname(abspath(path))
        out = joinpath(dir, "$(base)_$(k_actual)zones_splitbus.m")
    else
        out = String(save_path)
    end

    Lirion.write_matpower(case, out)
    return out
end

using Random

end
