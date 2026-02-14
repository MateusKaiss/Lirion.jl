############################################################
# Particionador METIS + Refinamento (BUS-ADMM / edge-cut / hybrid)
#  - Refinamento SEM copy(zone) (delta local)
#  - objective = :bus_admm | :edge_cut | :hybrid  (default :hybrid)
#  - Métricas ADMM-friendly: #shared, boundary-degree, cut-strength, edge_cut_ratio
#  - Modo "shared partition": se houver 2+ arquivos no glob, o 1º vira referência
#    e os demais recebem os MESMOS zones por bus_id (API/SAD/STD idênticos).
############################################################

using Printf, LinearAlgebra, SparseArrays, Random, Glob
using Lirion, Statistics
using Metis
using SimpleWeightedGraphs
using Graphs


const DEFAULT_DIR = normpath(@__DIR__, "..", "test", "matpower_cases")
function branch_counts_by_zone(case::Dict{String,Any}, zone::Vector{Int})
    buses = case["bus"]["bus_i"]
    pos = Dict(buses[i] => i for i in eachindex(buses))

    f = case["branch"]["f_bus"]
    t = case["branch"]["t_bus"]
    on = case["branch"]["br_status"]

    k = maximum(zone)
    internal = zeros(Int, k)      # linhas 100% dentro da zona
    incident = zeros(Int, k)      # linhas que tocam a zona (tie-line conta pros dois lados)
    cut = 0

    for i in eachindex(f)
        on[i] == 1 || continue
        zf = zone[pos[f[i]]]
        zt = zone[pos[t[i]]]

        incident[zf] += 1
        incident[zt] += 1

        if zf == zt
            internal[zf] += 1
        else
            cut += 1
        end
    end

    return internal, incident, cut
end

# ---------------------------
# utilidades
# ---------------------------
function require_lirion()
    try
        @eval using Lirion
    catch e
        error("""
        Não consegui carregar `Lirion`.
        Certifique-se de:
          1) Estar no ambiente certo (julia --project=.)
          2) Ter `Lirion.jl` no ambiente (ex: Pkg.develop(path=".../Lirion.jl"))
        Erro: $e
        """)
    end
end

load_case(path::String) = Lirion.load_input_data(path)

function has_multiple_zones(case::Dict{String,Any})
    if haskey(case, "bus") && haskey(case["bus"], "zone")
        z = case["bus"]["zone"]
        return length(unique(z)) > 1
    end
    return false
end

# ---------------------------
# writer MATPOWER com zones
# ---------------------------
function write_matpower_with_zones(path_out::String, case::Dict{String,Any})
    open(path_out, "w") do io
        name = get(case, "name", "case_generated")
        baseMVA = case["baseMVA"]
        println(io, "function mpc = $name")
        println(io, "mpc.version = '2';")
        println(io, "mpc.baseMVA = $(baseMVA);")

        # BUS
        println(io, "mpc.bus = [")
        b = case["bus"]
        @assert haskey(b, "zone") "Coluna 'zone' não encontrada em bus. Aplique a partição antes de salvar."
        for i in eachindex(b["bus_i"])
            @printf(io, " %d %d %.10g %.10g %.10g %.10g %d %.10g %.10g %.10g %d %.10g %.10g 0 0 0 0;\n",
                b["bus_i"][i], b["bus_type"][i], b["pd"][i], b["qd"][i], b["gs"][i], b["bs"][i],
                b["area"][i], b["vm"][i], b["va"][i], b["base_kv"][i], b["zone"][i],
                b["vmax"][i], b["vmin"][i])
        end
        println(io, "];")

        # GEN
        println(io, "mpc.gen = [")
        if haskey(case, "gen")
            g = case["gen"]
            for i in eachindex(g["gen_bus"])
                @printf(io, " %d %.10g %.10g %.10g %.10g %.10g %.10g %d %.10g %.10g 0 0 0 0 0 0 0 0 0 0;\n",
                    g["gen_bus"][i], g["pg"][i], g["qg"][i], g["qmax"][i], g["qmin"][i], g["vg"][i],
                    g["mbase"][i], g["gen_status"][i], g["pmax"][i], g["pmin"][i])
            end
        end
        println(io, "];")

        # BRANCH (angmin/angmax → graus)  ⚠️ garanta que o loader guarda em rad
        println(io, "mpc.branch = [")
        br = case["branch"]
        for i in eachindex(br["f_bus"])
            angmin_deg = br["angmin"][i] * 180/π
            angmax_deg = br["angmax"][i] * 180/π
            @printf(io, " %d %d %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %d %.10g %.10g 0 0 0 0 0 0 0;\n",
                br["f_bus"][i], br["t_bus"][i], br["br_r"][i], br["br_x"][i], br["br_b"][i],
                br["rate_a"][i], br["rate_b"][i], br["rate_c"][i], br["tap"][i], br["shift"][i],
                br["br_status"][i], angmin_deg, angmax_deg)
        end
        println(io, "];")

        # GENCOST
        println(io, "mpc.gencost = [")
        if haskey(case, "gencost")
            gc = case["gencost"]
            for i in eachindex(gc["model"])
                coefs = gc["cost_coefficients"][i]
                @printf(io, " %d %.10g %.10g %d", gc["model"][i], gc["startup"][i], gc["shutdown"][i], gc["ncost"][i])
                for c in coefs
                    @printf(io, " %.10g", c)
                end
                println(io, ";")
            end
        end
        println(io, "];")
    end
    return path_out
end

# ------------------------------------------------------------------
# Métricas ADMM-friendly
# ------------------------------------------------------------------
function count_shared_buses(case::Dict{String,Any}, zone::Vector{Int})
    buses = case["bus"]["bus_i"]
    pos = Dict(buses[i] => i for i in eachindex(buses))
    f = case["branch"]["f_bus"]; t = case["branch"]["t_bus"]; on = case["branch"]["br_status"]
    shared = falses(length(buses))
    for i in eachindex(f)
        on[i] == 1 || continue
        u = pos[f[i]]; v = pos[t[i]]
        if zone[u] != zone[v]
            shared[u] = true
            shared[v] = true
        end
    end
    return count(shared), shared
end

function edge_cut_strength(case::Dict{String,Any}, zone::Vector{Int})
    f = case["branch"]["f_bus"]; t = case["branch"]["t_bus"]; on = case["branch"]["br_status"]
    x = case["branch"]["br_x"]; tap = get(case["branch"], "tap", fill(1.0, length(f)))
    shift = get(case["branch"], "shift", fill(0.0, length(f)))
    ratea = get(case["branch"], "rate_a", fill(0.0, length(f)))
    buses = case["bus"]["bus_i"]; pos = Dict(buses[i]=>i for i in eachindex(buses))
    S = 0.0
    for i in eachindex(f)
        on[i]==1 || continue
        if zone[pos[f[i]]] != zone[pos[t[i]]]
            w = (abs(x[i]) > 1e-10 ? 1/abs(x[i]) : 1e6)
            if abs(tap[i]-1.0) > 1e-6 || abs(shift[i]) > 1e-6; w *= 50.0; end
            if ratea[i] > 0; w *= (1 + ratea[i]); end
            S += w
        end
    end
    return S
end

function edge_cut_ratio(case::Dict{String,Any}, zone::Vector{Int})
    f = case["branch"]["f_bus"]; t = case["branch"]["t_bus"]
    x = case["branch"]["br_x"]; on = case["branch"]["br_status"]
    tap = get(case["branch"], "tap", fill(1.0, length(f)))
    shift = get(case["branch"], "shift", fill(0.0, length(f)))
    ratea = get(case["branch"], "rate_a", fill(0.0, length(f)))
    buses = case["bus"]["bus_i"]; pos = Dict(buses[i]=>i for i in eachindex(buses))

    w_tot = 0.0; w_cut = 0.0
    for i in eachindex(f)
        on[i]==1 || continue
        u, v = f[i], t[i]
        w = (abs(x[i]) > 1e-10 ? 1/abs(x[i]) : 1e6)
        if abs(tap[i]-1.0) > 1e-6 || abs(shift[i]) > 1e-6; w *= 50.0; end
        if ratea[i] > 0; w *= (1 + ratea[i]); end
        w_tot += w
        if zone[pos[u]] != zone[pos[v]]
            w_cut += w
        end
    end
    return w_cut / max(w_tot, eps())
end

# ------------------------------------------------------------------
# Pós-processos: conectividade + geração/slack
# ------------------------------------------------------------------
function fix_connectivity_and_slack!(case::Dict{String,Any}, zone::Vector{Int}, k::Int)
    buses = case["bus"]["bus_i"]
    f = case["branch"]["f_bus"]; t = case["branch"]["t_bus"]; on = case["branch"]["br_status"]

    neigh = Dict{Int,Vector{Int}}()
    for i in eachindex(f)
        on[i]==1 || continue
        push!(get!(neigh, f[i], Int[]), t[i])
        push!(get!(neigh, t[i], Int[]), f[i])
    end
    pos = Dict(buses[i]=>i for i in eachindex(buses))

    # conectividade por zona
    for z in 1:k
        nodes = [buses[i] for i in eachindex(buses) if zone[i]==z]
        isempty(nodes) && continue
        seen = Set{Int}(); comps = Vector{Vector{Int}}()
        for u in nodes
            u in seen && continue
            comp = Int[]; Q = [u]; push!(seen,u)
            while !isempty(Q)
                w = popfirst!(Q); push!(comp, w)
                for v in get(neigh, w, Int[])
                    if zone[pos[v]]==z && !(v in seen)
                        push!(seen, v); push!(Q, v)
                    end
                end
            end
            push!(comps, comp)
        end
        length(comps) <= 1 && continue
        sort!(comps, by=length, rev=true)
        for comp in comps[2:end]
            border = [(u,v) for u in comp for v in get(neigh,u,Int[]) if zone[pos[v]] != z]
            isempty(border) && continue
            counts = Dict{Int,Int}()
            for (_,v) in border
                zv = zone[pos[v]]
                counts[zv] = get(counts, zv, 0) + 1
            end
            zbest = let ks = collect(keys(counts)); ks[argmax(collect(values(counts)))] end
            for u in comp
                zone[pos[u]] = zbest
            end
        end
    end

    # ao menos um gerador por zona (pmax>0)
    gen = get(case, "gen", Dict{String,Any}())
    isempty(gen) && return zone

    genbus = gen["gen_bus"]; status = gen["gen_status"]; pmax = gen["pmax"]
    pos2 = Dict(case["bus"]["bus_i"][i]=>i for i in eachindex(case["bus"]["bus_i"]))

    pmax_by_zone = Dict(z=>0.0 for z in 1:k)
    for i in eachindex(genbus)
        status[i]==1 || continue
        b = genbus[i]; z = zone[pos2[b]]
        pmax_by_zone[z] = get(pmax_by_zone, z, 0.0) + max(pmax[i],0.0)
    end

    neighZ = Dict{Int,Set{Int}}()
    for i in eachindex(f)
        on[i]==1 || continue
        zf = zone[pos[f[i]]]; zt = zone[pos[t[i]]]
        if zf != zt
            push!(get!(neighZ, zf, Set{Int}()), zt)
            push!(get!(neighZ, zt, Set{Int}()), zf)
        end
    end

    gens_in_zone = Dict(z=>Int[] for z in 1:k)
    for i in eachindex(genbus)
        status[i]==1 || continue
        if pmax[i] > 1e-6
            bz = zone[pos2[genbus[i]]]
            push!(gens_in_zone[bz], i)
        end
    end

    for z in 1:k
        get(pmax_by_zone, z, 0.0) > 1e-9 && continue
        moved = false
        for zv in sort(collect(get(neighZ, z, Set{Int}())))
            if !isempty(get(gens_in_zone, zv, Int[]))
                gi = pop!(gens_in_zone[zv]); b = genbus[gi]
                zone[pos2[b]] = z
                pmax_by_zone[z]  = get(pmax_by_zone, z, 0.0) + max(pmax[gi],0.0)
                pmax_by_zone[zv] = max(get(pmax_by_zone, zv, 0.0) - max(pmax[gi],0.0), 0.0)
                push!(gens_in_zone[z], gi)
                moved = true
                break
            end
        end
        moved && continue
        for zz in 1:k
            zz == z && continue
            if !isempty(get(gens_in_zone, zz, Int[]))
                gi = pop!(gens_in_zone[zz]); b = genbus[gi]
                zone[pos2[b]] = z
                pmax_by_zone[z]  += max(pmax[gi],0.0)
                pmax_by_zone[zz]  = max(pmax_by_zone[zz] - max(pmax[gi],0.0), 0.0)
                push!(gens_in_zone[z], gi)
                break
            end
        end
    end
    return zone
end

function apply_zones!(case::Dict{String,Any}, zone::Vector{Int})
    buses = case["bus"]["bus_i"]
    @assert length(zone) == length(buses)
    case["bus"]["zone"] = copy(zone)
    return case
end

# ------------------------------------------------------------------
# METIS particionamento
# ------------------------------------------------------------------
function metis_partition(case::Dict{String,Any}; k::Int)
    @assert k >= 1
    @assert haskey(case, "bus") && haskey(case, "branch")

    buses  = case["bus"]["bus_i"]
    present = Set(buses)
    f   = case["branch"]["f_bus"]; t = case["branch"]["t_bus"]
    x   = case["branch"]["br_x"];  on  = case["branch"]["br_status"]
    tap   = get(case["branch"], "tap",   fill(1.0, length(f)))
    shift = get(case["branch"], "shift", fill(0.0, length(f)))
    ratea = get(case["branch"], "rate_a", fill(0.0, length(f)))

    glob2loc = Dict{Int,Int}(); n = length(buses)
    for (i,b) in enumerate(buses)
        glob2loc[b] = i
    end

    # pesos normalizados por aresta (acumulando paralelas)
    raw_w = Float64[]
    for i in eachindex(f)
        on[i]==1 || continue
        (f[i] in present && t[i] in present) || continue
        w = (abs(x[i]) > 1e-10 ? 1/abs(x[i]) : 1e6)
        if abs(tap[i]-1.0) > 1e-6 || abs(shift[i]) > 1e-6
            w *= 50.0
        end
        if ratea[i] > 0
            w *= (1 + ratea[i])
        end
        push!(raw_w, w)
    end
    isempty(raw_w) && (raw_w = [1.0])
    maxw = maximum(raw_w)

    wdict = Dict{Tuple{Int,Int}, Int}()
    for i in eachindex(f)
        on[i]==1 || continue
        ug, vg = f[i], t[i]
        (ug in present && vg in present) || continue
        wfloat = (abs(x[i]) > 1e-10 ? 1/abs(x[i]) : 1e6)
        if abs(tap[i]-1.0) > 1e-6 || abs(shift[i]) > 1e-6
            wfloat *= 50.0
        end
        if ratea[i] > 0
            wfloat *= (1 + ratea[i])
        end
        w = max(1, Int(round(1 + 10_000 * (wfloat / maxw))))
        u = glob2loc[ug]; v = glob2loc[vg]
        key = u < v ? (u,v) : (v,u)
        wdict[key] = get(wdict, key, 0) + w
    end

    G = SimpleWeightedGraph(n)
    for ((u,v), w) in wdict
        add_edge!(G, u, v, w)
    end

    parts = Metis.partition(G, k; alg = :KWAY)  # rótulos 1..k
    Z = fill(1, length(buses))
    @inbounds for (i, bid) in enumerate(buses)
        Z[i] = parts[glob2loc[bid]]
    end
    fix_connectivity_and_slack!(case, Z, k)
    return Z
end

# ------------------------------------------------------------------
# Refinamento sem copy(zone)
# ------------------------------------------------------------------
mutable struct BusADJ
    neigh::Vector{Vector{Int}}
    edge_w::Vector{Dict{Int,Float64}}
end

function build_adjacency_with_weights(case::Dict{String,Any})
    buses = case["bus"]["bus_i"]
    pos = Dict(buses[i]=>i for i in eachindex(buses))

    f = case["branch"]["f_bus"]; t = case["branch"]["t_bus"]; on = case["branch"]["br_status"]
    x = case["branch"]["br_x"]
    tap = get(case["branch"], "tap", fill(1.0, length(f)))
    shift = get(case["branch"], "shift", fill(0.0, length(f)))
    ratea = get(case["branch"], "rate_a", fill(0.0, length(f)))

    neigh = [Int[] for _ in buses]
    ew = [Dict{Int,Float64}() for _ in buses]

    for i in eachindex(f)
        on[i]==1 || continue
        u = pos[f[i]]; v = pos[t[i]]
        w = (abs(x[i]) > 1e-10 ? 1/abs(x[i]) : 1e6)
        if abs(tap[i]-1.0) > 1e-6 || abs(shift[i]) > 1e-6
            w *= 50.0
        end
        if ratea[i] > 0
            w *= (1 + ratea[i])
        end
        push!(neigh[u], v); push!(neigh[v], u)
        ew[u][v] = get(ew[u], v, 0.0) + w
        ew[v][u] = get(ew[v], u, 0.0) + w
    end

    return BusADJ(neigh, ew)
end

# Conectividade (origem) após remover 1 nó
function zone_connected_after_remove(zone::Vector{Int}, z::Int, node::Int, adj::BusADJ)
    nodes = Int[]
    for (i, zi) in enumerate(zone)
        if zi == z && i != node
            push!(nodes, i)
        end
    end
    isempty(nodes) && return true

    seen = Set{Int}([nodes[1]])
    Q = [nodes[1]]
    while !isempty(Q)
        u = popfirst!(Q)
        for v in adj.neigh[u]
            if v != node && zone[v] == z && !(v in seen)
                push!(seen, v); push!(Q, v)
            end
        end
    end
    return length(seen) == length(nodes)
end

# pmax por bus
function pmax_per_bus(case::Dict{String,Any})
    g = get(case, "gen", Dict{String,Any}())
    isempty(g) && return zeros(length(case["bus"]["bus_i"]))

    buses = case["bus"]["bus_i"]
    pos = Dict(buses[i]=>i for i in eachindex(buses))

    genbus = g["gen_bus"]; status = g["gen_status"]; pmax = g["pmax"]
    v = zeros(Float64, length(buses))
    for i in eachindex(genbus)
        status[i]==1 || continue
        v[pos[genbus[i]]] += max(pmax[i], 0.0)
    end
    return v
end

# shared(v) "antes" e "depois" quando só u muda de zona
@inline function shared_before_after_for_neighbor(v::Int, u::Int,
    z_u_before::Int, z_u_after::Int, zone::Vector{Int}, adj::BusADJ)

    zv = zone[v]
    has_other_cross = false
    for w in adj.neigh[v]
        w == u && continue
        if zone[w] != zv
            has_other_cross = true
            break
        end
    end
    before = has_other_cross || (z_u_before != zv)
    after  = has_other_cross || (z_u_after  != zv)
    return before, after
end

"""
delta ao mover u para znew:
  Δshared  = mudança no #shared global, computada localmente (u + vizinhos)
  Δbdeg    = mudança no boundary-degree (soma do grau dos shared)
  Δcut     = mudança no cut-strength local (arestas incidentes a u)
Sem copy(zone).
"""
function local_deltas_move_nocopy(u::Int, znew::Int, zone::Vector{Int}, adj::BusADJ)
    zold = zone[u]
    zold == znew && return (0, 0, 0.0)

    # u: shared antes/depois
    shared_before_u = false
    shared_after_u  = false
    for v in adj.neigh[u]
        (zone[v] != zold) && (shared_before_u = true)
        (zone[v] != znew) && (shared_after_u  = true)
        (shared_before_u && shared_after_u) && break
    end
    Δshared_u = (shared_after_u ? 1 : 0) - (shared_before_u ? 1 : 0)

    deg_u = length(adj.neigh[u])
    Δbdeg_u = (shared_after_u ? deg_u : 0) - (shared_before_u ? deg_u : 0)

    # vizinhos: só eles podem mudar shared por causa de u
    Δshared_neigh = 0
    Δbdeg_neigh = 0
    for v in adj.neigh[u]
        before, after = shared_before_after_for_neighbor(v, u, zold, znew, zone, adj)
        if before && !after
            Δshared_neigh -= 1
            Δbdeg_neigh   -= length(adj.neigh[v])
        elseif !before && after
            Δshared_neigh += 1
            Δbdeg_neigh   += length(adj.neigh[v])
        end
    end

    # cut strength local (arestas incidentes em u)
    cut_before = 0.0
    cut_after  = 0.0
    for v in adj.neigh[u]
        w = adj.edge_w[u][v]
        (zone[v] != zold) && (cut_before += w)
        (zone[v] != znew) && (cut_after  += w)
    end
    Δcut = cut_after - cut_before

    return (Δshared_u + Δshared_neigh, Δbdeg_u + Δbdeg_neigh, Δcut)
end

# score por objetivo (default = :hybrid)
@inline function move_score(objective::Symbol, Δs::Int, Δbd::Int, Δc::Float64;
    w_shared::Float64, w_bdeg::Float64, w_cut::Float64)

    if objective == :bus_admm
        # lexicográfico: shared -> boundary-degree -> cut
        return (Δs, Δbd, Δc)
    elseif objective == :edge_cut
        # lexicográfico: cut -> shared -> boundary-degree
        return (Δc, Δs, Δbd)
    else
        # escalar: shared e boundary-degree pesam muito; cut refina
        return (w_shared*Δs + w_bdeg*Δbd + w_cut*Δc)
    end
end

"""
refine_partition!(case, zone; objective=:hybrid, passes=4, max_zone_size=..., verbose=true, weights...)
"""
function refine_partition!(case::Dict{String,Any}, zone::Vector{Int};
    objective::Symbol = :hybrid,
    passes::Int = 4,
    max_zone_size::Int = typemax(Int),
    verbose::Bool = true,
    w_shared::Float64 = 1.0,
    w_bdeg::Float64 = 0.15,
    w_cut::Float64 = 0.002,
)
    k = maximum(zone)
    adj = build_adjacency_with_weights(case)

    # estruturas por zona
    sizeZ = [count(==(z), zone) for z in 1:k]
    pmax_bus = pmax_per_bus(case)
    pmaxZ = [0.0 for _ in 1:k]
    for (i, z) in enumerate(zone)
        pmaxZ[z] += pmax_bus[i]
    end

    zone_neighbors(u) = unique([zone[v] for v in adj.neigh[u] if zone[v] != zone[u]])

    for pass in 1:passes
        n_shared, shared_mask = count_shared_buses(case, zone)
        cutS = edge_cut_strength(case, zone)

        degs = [length(adj.neigh[i]) for i in eachindex(zone)]
        bdeg_sum = sum(degs[i] for i in eachindex(zone) if shared_mask[i])
        bdeg_max = maximum([shared_mask[i] ? degs[i] : 0 for i in eachindex(zone)])
        bdeg_avg = (n_shared > 0) ? (bdeg_sum / n_shared) : 0.0

        verbose && @info @sprintf(
            "[refine pass %d] shared=%d | bdeg_sum=%d (avg=%.2f, max=%d) | cut=%.2f",
            pass, n_shared, bdeg_sum, bdeg_avg, bdeg_max, cutS
        )

        candidates = findall(shared_mask)
        isempty(candidates) && break
        shuffle!(candidates)

        moved_any = false

        for u in candidates
            z0 = zone[u]

            # não mover se for o único “pmax útil” da zona
            if pmax_bus[u] > 0 && pmaxZ[z0] - pmax_bus[u] <= 1e-9
                continue
            end

            nbr_zones = zone_neighbors(u)
            isempty(nbr_zones) && continue

            # conectividade: checar uma vez (origem) para este u
            if !zone_connected_after_remove(zone, z0, u, adj)
                continue
            end

            best_z = z0
            best_score = nothing

            for zt in nbr_zones
                if sizeZ[zt] + 1 > max_zone_size || sizeZ[z0] - 1 <= 0
                    continue
                end

                Δs, Δbd, Δc = local_deltas_move_nocopy(u, zt, zone, adj)
                sc = move_score(objective, Δs, Δbd, Δc; w_shared=w_shared, w_bdeg=w_bdeg, w_cut=w_cut)

                if best_score === nothing || sc < best_score
                    best_score = sc
                    best_z = zt
                end
            end

            if best_z != z0
                zone[u] = best_z
                sizeZ[z0] -= 1; sizeZ[best_z] += 1
                pmaxZ[z0] = max(pmaxZ[z0] - pmax_bus[u], 0.0)
                pmaxZ[best_z] += pmax_bus[u]
                moved_any = true
            end
        end

        fix_connectivity_and_slack!(case, zone, k)

        if !moved_any
            verbose && @info "[refine pass $pass] sem movimentos benéficos."
            break
        end
    end

    return zone
end

# ------------------------------------------------------------------
# Escolha automática de k
# ------------------------------------------------------------------
function choose_k(case; max_zone::Int=5000, kmin::Int=2, kmax::Int=2000)
    n = length(case["bus"]["bus_i"])
    k = max(kmin, ceil(Int, n / max_zone))
    return min(max(k, kmin), kmax)
end

# ============================================================
# Replicar a MESMA partição (zones por bus_id) em vários casos
# ============================================================

function build_bus_zone_map(case_ref::Dict{String,Any}, zone_ref::Vector{Int})
    buses = case_ref["bus"]["bus_i"]
    @assert length(buses) == length(zone_ref)
    m = Dict{Int,Int}()
    @inbounds for i in eachindex(buses)
        m[buses[i]] = zone_ref[i]
    end
    return m
end

function apply_bus_zone_map!(
    case_target::Dict{String,Any},
    bus2zone::Dict{Int,Int};
    strict::Bool = true,
    default_zone::Int = 1,
)
    buses = case_target["bus"]["bus_i"]
    zone = Vector{Int}(undef, length(buses))
    missing = Int[]

    @inbounds for i in eachindex(buses)
        b = buses[i]
        if haskey(bus2zone, b)
            zone[i] = bus2zone[b]
        else
            push!(missing, b)
            zone[i] = default_zone
        end
    end

    if strict && !isempty(missing)
        error("Buses sem zona no mapeamento de referência: $(missing[1:min(end,20)])" *
              (length(missing) > 20 ? " ..." : ""))
    end

    uniq = sort(unique(zone))
    relabel = Dict(uniq[i] => i for i in eachindex(uniq))
    @inbounds for i in eachindex(zone)
        zone[i] = relabel[zone[i]]
    end

    case_target["bus"]["zone"] = copy(zone)
    return zone
end

"""
partition_and_reuse_zones:
- Se 1 arquivo: particiona+refina e salva.
- Se 2+: particiona+refina o 1º e aplica nos demais via bus_id=>zone.
"""
function partition_and_reuse_zones(indir::AbstractString=DEFAULT_DIR;
    glob::AbstractString = "*.m",
    overwrite::Bool = false,
    ks::Union{Vector{Int},Symbol} = :auto,
    max_zone_size::Int = 5000,
    refine_passes::Int = 4,
    objective::Symbol = :hybrid,
    w_shared::Float64 = 1.0,
    w_bdeg::Float64 = 0.15,
    w_cut::Float64 = 0.002,
    strict::Bool = true,
    seed::Integer = 42,
)
    require_lirion()
    Random.seed!(seed)

    files = sort(Glob.glob(glob, indir))
    isempty(files) && error("Nenhum arquivo encontrado: glob=$glob em $indir")

    ks_list_for(case) = (ks === :auto) ? Int[choose_k(case; max_zone=max_zone_size)] : collect(Int, ks)

    # --- 1 arquivo: pipeline normal
    if length(files) == 1
        f = files[1]
        base = basename(f)
        case = load_case(f)

        for k in ks_list_for(case)
            outdir = joinpath(indir, "case $(k) zones (metis-refined)")
            isdir(outdir) || mkpath(outdir)
            outpath = joinpath(outdir, base)
            if isfile(outpath) && !overwrite
                @info "PULANDO (destino existe): $outpath"
                continue
            end

            z = metis_partition(case; k=k)
            z = refine_partition!(case, z;
                objective=objective,
                passes=refine_passes,
                max_zone_size=max_zone_size,
                verbose=true,
                w_shared=w_shared, w_bdeg=w_bdeg, w_cut=w_cut,
            )

            apply_zones!(case, z)
            write_matpower_with_zones(outpath, case)
            @info "OK (1 arquivo): $base → $outpath"
        end
        return nothing
    end

    # --- 2+ arquivos: 1º é referência
    f_ref = files[1]
    base_ref = basename(f_ref)
    case_ref = load_case(f_ref)

    for k in ks_list_for(case_ref)
        @info "Referência: $base_ref (k=$k) — gerando partição única para reaproveitar"

        z_ref = metis_partition(case_ref; k=k)
        z_ref = refine_partition!(case_ref, z_ref;
            objective=objective,
            passes=refine_passes,
            max_zone_size=max_zone_size,
            verbose=true,
            w_shared=w_shared, w_bdeg=w_bdeg, w_cut=w_cut,
        )
        internal, incident, cut = branch_counts_by_zone(case_ref, z_ref)
        @info "Branches: cut=$cut"
        for z in 1:maximum(z_ref)
            @info "Zone $z: internal=$(internal[z]) | incident=$(incident[z])"
        end


        bus2zone = build_bus_zone_map(case_ref, z_ref)

        outdir = joinpath(indir, "case $(k) zones (metis-refined-shared)")
        isdir(outdir) || mkpath(outdir)

        # salva referência
        out_ref = joinpath(outdir, base_ref)
        if !(isfile(out_ref) && !overwrite)
            apply_zones!(case_ref, z_ref)
            write_matpower_with_zones(out_ref, case_ref)
            @info "OK (ref): $base_ref → $out_ref"
        else
            @info "PULANDO (ref destino existe): $out_ref"
        end

        # aplica nos demais
        for f in files[2:end]
            base = basename(f)
            outpath = joinpath(outdir, base)
            if isfile(outpath) && !overwrite
                @info "PULANDO (destino existe): $outpath"
                continue
            end

            case_i = load_case(f)
            z_i = apply_bus_zone_map!(case_i, bus2zone; strict=strict)

            # garante conectividade/geração por zona (pode ajustar levemente zones em componentes soltos)
            fix_connectivity_and_slack!(case_i, z_i, maximum(z_i))

            write_matpower_with_zones(outpath, case_i)
            @info "OK (reuse zones): $base → $outpath"
        end
    end

    @info "Concluído: partição compartilhada para $(length(files)) arquivos."
    return nothing
end

# --------------------------------------------
# Execução (shared partition se houver 2+ arquivos no glob)
# --------------------------------------------
partition_and_reuse_zones(
    raw"/home/mateus-kaiss/Desktop/experiments";
    glob = "*.m",
    ks = [30],
    max_zone_size = 200,
    refine_passes = 14,
    objective = :hybrid,
    w_shared = 2.0,
    w_bdeg = 0.30,
    w_cut = 0.001,
    strict = true,
    overwrite = false,
    seed = 42,
)

# partition_and_reuse_zones(
#     raw"/home/mateus-kaiss/Desktop/experiments";
#     glob = "*.m",
#     ks = [5],
#     max_zone_size = 200,
#     refine_passes = 3,      # <- melhor que 6 na maioria dos casos
#     objective = :hybrid,
#     w_shared = 2.0,
#     w_bdeg = 0.30,
#     w_cut = 0.01,
#     strict = true,
#     overwrite = false,
#     seed = 42,
# )


# using Lirion, Statistics

# file = raw"/home/mateus-kaiss/Desktop/experiments/case 3 zones (metis-refined-shared)/pglib_opf_case73_ieee_rts__sad.m"
# case = Lirion.load_input_data(file)

# z = case["bus"]["zone"]          # vetor de zonas por bus (mesma ordem de bus_i)
# bus_i = case["bus"]["bus_i"]

# for k in sort(unique(z))
#     nb = count(==(k), z)
#     println("Zona $k: $nb buses")
# end
# println("Total buses = ", length(bus_i))
# # 