# ============================================================
# paper/scripts/lib/common.jl
# ============================================================
# Shared helpers for the prepare scripts (01, 02). Pure data/topology
# utilities + PGLib path resolution. No solver.
# ============================================================

using Lirion
using TOML
using Dates
using Printf

# ---- Paths ------------------------------------------------------------

const PAPER_DIR        = normpath(joinpath(@__DIR__, "..", ".."))
const CONFIG_DIR       = joinpath(PAPER_DIR, "configs")
const CASES_DIR        = joinpath(PAPER_DIR, "cases")
const RAW_DIR          = joinpath(CASES_DIR, "raw", "pglib")
const PARTITIONED_DIR  = joinpath(CASES_DIR, "partitioned")
const RESULTS_DIR      = joinpath(PAPER_DIR, "results")
const LOG_DIR          = joinpath(RESULTS_DIR, "logs")

"""Resolve the PGLib-OPF root: env var `PGLIB_OPF_ROOT`, else the default."""
function pglib_root()
    default = "/home/mateus-kaiss/Downloads/papers/pglib-opf-master"
    root = get(ENV, "PGLIB_OPF_ROOT", default)
    isdir(root) || error("PGLib root not found: $root (set PGLIB_OPF_ROOT)")
    return root
end

"""Source `.m` path for a given case name + variant in the PGLib tree."""
function pglib_variant_path(root::AbstractString, name::AbstractString, variant::AbstractString)
    if variant == "std"
        return joinpath(root, "pglib_opf_$(name).m")
    elseif variant == "api"
        return joinpath(root, "api", "pglib_opf_$(name)__api.m")
    elseif variant == "sad"
        return joinpath(root, "sad", "pglib_opf_$(name)__sad.m")
    else
        error("unknown variant: $variant")
    end
end

const VARIANTS = ("std", "api", "sad")

# ---- Config loading ---------------------------------------------------

load_selected_cases() = TOML.parsefile(joinpath(CASES_DIR, "selected_cases.toml"))
load_partitioning_config() = TOML.parsefile(joinpath(CONFIG_DIR, "partitioning.toml"))

# ---- Topology signature ----------------------------------------------

"""
Canonical topology fingerprint: sorted bus ids + sorted undirected
in-service branch endpoint pairs. Same fingerprint => safe to share a
bus->zone map across variants.
"""
function topology_signature(case::Lirion.CaseRaw)
    bus_ids = sort(collect(case.bus.bus_i))
    edges = Tuple{Int,Int}[]
    br = case.branch
    for i in eachindex(br.f_bus)
        br.br_status[i] == 0 && continue
        f, t = br.f_bus[i], br.t_bus[i]
        push!(edges, (min(f, t), max(f, t)))
    end
    sort!(edges)
    return (bus_ids = bus_ids, edges = edges, hash = hash((bus_ids, edges)))
end

topology_equal(a, b) = (a.bus_ids == b.bus_ids && a.edges == b.edges)

# ---- Counts -----------------------------------------------------------

n_bus(c::Lirion.CaseRaw)    = length(c.bus.bus_i)
n_branch(c::Lirion.CaseRaw) = length(c.branch.f_bus)
n_gen(c::Lirion.CaseRaw)    = c.gen === nothing ? 0 : length(c.gen.gen_bus)

# ---- Partition metrics (from bus.zone + branch topology) -------------

function bus_zone_map(c::Lirion.CaseRaw)
    d = Dict{Int,Int}()
    for i in eachindex(c.bus.bus_i)
        d[c.bus.bus_i[i]] = c.bus.zone[i]
    end
    return d
end

"""Tie-lines (cut edges): in-service branches whose ends are in different zones."""
function tie_lines(c::Lirion.CaseRaw)
    z = bus_zone_map(c)
    ties = Vector{Vector{Int}}()   # [f, t, zone_f, zone_t]
    br = c.branch
    for i in eachindex(br.f_bus)
        br.br_status[i] == 0 && continue
        f, t = br.f_bus[i], br.t_bus[i]
        z[f] == z[t] && continue
        push!(ties, [f, t, z[f], z[t]])
    end
    return ties
end

function boundary_buses(ties::Vector{Vector{Int}})
    s = Set{Int}()
    for e in ties
        push!(s, e[1]); push!(s, e[2])
    end
    return sort!(collect(s))
end

function zone_members(c::Lirion.CaseRaw)
    d = Dict{Int,Vector{Int}}()
    for i in eachindex(c.bus.bus_i)
        push!(get!(d, c.bus.zone[i], Int[]), c.bus.bus_i[i])
    end
    for v in values(d); sort!(v); end
    return d
end

function balance_metrics(members::Dict{Int,Vector{Int}})
    sizes = [length(v) for v in values(members)]
    n = length(sizes)
    mean = sum(sizes) / n
    var = n > 1 ? sum((s - mean)^2 for s in sizes) / (n - 1) : 0.0
    sd = sqrt(var)
    return Dict(
        "n_zones" => n,
        "min" => minimum(sizes),
        "max" => maximum(sizes),
        "mean" => mean,
        "std" => sd,
        "cv" => mean > 0 ? sd / mean : 0.0,
        "imbalance_max_over_mean" => mean > 0 ? maximum(sizes) / mean : 0.0,
        "sizes_by_zone" => Dict(string(z) => length(v) for (z, v) in members),
    )
end

"""Per-zone connectivity using only intra-zone in-service branches (BFS)."""
function zone_connectivity(c::Lirion.CaseRaw)
    z = bus_zone_map(c)
    members = zone_members(c)
    br = c.branch
    adj = Dict{Int,Vector{Int}}()
    for b in c.bus.bus_i; adj[b] = Int[]; end
    for i in eachindex(br.f_bus)
        br.br_status[i] == 0 && continue
        f, t = br.f_bus[i], br.t_bus[i]
        z[f] == z[t] || continue
        push!(adj[f], t); push!(adj[t], f)
    end
    result = Dict{String,Bool}()
    all_conn = true
    for (zone, buses) in members
        if length(buses) <= 1
            result[string(zone)] = true
            continue
        end
        seen = Set{Int}(); stack = [buses[1]]; push!(seen, buses[1])
        while !isempty(stack)
            u = pop!(stack)
            for v in adj[u]
                if v ∉ seen
                    push!(seen, v); push!(stack, v)
                end
            end
        end
        conn = length(seen) == length(buses)
        result[string(zone)] = conn
        all_conn &= conn
    end
    return (per_zone = result, all_connected = all_conn)
end

"""Per-zone generation presence (dispatchable: status==1 && pmax>0)."""
function zone_generation(c::Lirion.CaseRaw)
    z = bus_zone_map(c)
    members = zone_members(c)
    has = Dict(string(zone) => false for zone in keys(members))
    if c.gen !== nothing
        for i in eachindex(c.gen.gen_bus)
            c.gen.gen_status[i] == 1 || continue
            c.gen.pmax[i] > 0 || continue
            gb = c.gen.gen_bus[i]
            haskey(z, gb) && (has[string(z[gb])] = true)
        end
    end
    return (per_zone = has, all_zones_have_gen = all(values(has)))
end

# ---- TOML writer ------------------------------------------------------

function write_toml(path::AbstractString, data::AbstractDict)
    mkpath(dirname(path))
    open(path, "w") do io
        TOML.print(io, data; sorted = true)
    end
    return path
end

function log_line(io, msg)
    stamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    println(io, "[$stamp] $msg"); flush(io)
    println(msg)
end
