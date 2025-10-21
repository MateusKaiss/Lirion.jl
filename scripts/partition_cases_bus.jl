############################################################
# Particionador por METIS (foco BUS-ADMM)
############################################################
using Printf, LinearAlgebra, SparseArrays, Graphs, Random, Glob
using Lirion, Statistics
using Metis
using SimpleWeightedGraphs

const DEFAULT_DIR = normpath(@__DIR__, "..", "test", "matpower_cases")

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
        z = case["bus"]["zone"]; return length(unique(z)) > 1
    end
    return false
end

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

        # BRANCH  (angmin/angmax → graus)
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
# Métricas
# ------------------------------------------------------------------
# nós de fronteira (barramentos com vizinho em outra zona)
function count_shared_buses(case::Dict{String,Any}, zone::Vector{Int})
    buses = case["bus"]["bus_i"]
    pos = Dict(buses[i]=>i for i in eachindex(buses))
    f = case["branch"]["f_bus"]; t = case["branch"]["t_bus"]; on = case["branch"]["br_status"]
    shared = falses(length(buses))
    for i in eachindex(f)
        on[i]==1 || continue
        u=pos[f[i]]; v=pos[t[i]]
        if zone[u] != zone[v]
            shared[u] = true; shared[v] = true
        end
    end
    return count(shared), shared
end

# força do corte global
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
# Base METIS (idêntica ao seu particionador “branch”)
# ------------------------------------------------------------------
function metis_partition(case::Dict{String,Any}; k::Int, balance_tol::Float64=0.10)
    @assert k >= 1
    @assert haskey(case, "bus") && haskey(case["bus"], "bus_i")
    @assert haskey(case, "branch")

    buses  = case["bus"]["bus_i"]
    present = Set(buses)
    f   = case["branch"]["f_bus"]; t = case["branch"]["t_bus"]
    x   = case["branch"]["br_x"];  on  = case["branch"]["br_status"]
    tap   = get(case["branch"], "tap",   fill(1.0, length(f)))
    shift = get(case["branch"], "shift", fill(0.0, length(f)))
    ratea = get(case["branch"], "rate_a", fill(0.0, length(f)))

    glob2loc = Dict{Int,Int}(); n = length(buses)
    for (i,b) in enumerate(buses) glob2loc[b] = i end

    # pesos normalizados por aresta (acumulando paralelas)
    raw_w = Float64[]
    for i in eachindex(f)
        on[i]==1 || continue
        (f[i] in present && t[i] in present) || continue
        w = (abs(x[i]) > 1e-10 ? 1/abs(x[i]) : 1e6)
        if abs(tap[i]-1.0) > 1e-6 || abs(shift[i]) > 1e-6; w *= 50.0; end
        if ratea[i] > 0; w *= (1 + ratea[i]); end
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
        if abs(tap[i]-1.0) > 1e-6 || abs(shift[i]) > 1e-6; wfloat *= 50.0; end
        if ratea[i] > 0; wfloat *= (1 + ratea[i]); end
        w = max(1, Int(round(1 + 10_000 * (wfloat / maxw))))
        u = glob2loc[ug]; v = glob2loc[vg]
        key = u < v ? (u,v) : (v,u)
        wdict[key] = get(wdict, key, 0) + w
    end

    G = SimpleWeightedGraph(n)
    for ((u,v), w) in wdict
        add_edge!(G, u, v, w)
    end

    parts = Metis.partition(G, k; alg = :KWAY)   # rótulos em 1..k

    Z = fill(1, length(buses))
    @inbounds for (i, bid) in enumerate(buses)
        Z[i] = parts[glob2loc[bid]]
    end

    fix_connectivity_and_slack!(case, Z, k)
    return Z
end

# ------------------------------------------------------------------
# Refinamento BUS-ADMM (greedy local)
# ------------------------------------------------------------------
# Estruturas auxiliares para deltas locais
mutable struct BusADJ
    neigh::Vector{Vector{Int}}   # vizinhos por índice local do bus
    edge_w::Vector{Dict{Int,Float64}} # pesos (por vizinho local)
end

function build_adjacency_with_weights(case::Dict{String,Any})
    buses = case["bus"]["bus_i"]; pos = Dict(buses[i]=>i for i in eachindex(buses))
    f = case["branch"]["f_bus"]; t = case["branch"]["t_bus"]; on = case["branch"]["br_status"]
    x = case["branch"]["br_x"]; tap = get(case["branch"], "tap", fill(1.0, length(f)))
    shift = get(case["branch"], "shift", fill(0.0, length(f)))
    ratea = get(case["branch"], "rate_a", fill(0.0, length(f)))

    neigh = [Int[] for _ in buses]
    ew = [Dict{Int,Float64}() for _ in buses]
    for i in eachindex(f)
        on[i]==1 || continue
        u = pos[f[i]]; v = pos[t[i]]
        w = (abs(x[i]) > 1e-10 ? 1/abs(x[i]) : 1e6)
        if abs(tap[i]-1.0) > 1e-6 || abs(shift[i]) > 1e-6; w *= 50.0; end
        if ratea[i] > 0; w *= (1 + ratea[i]); end
        push!(neigh[u], v); push!(neigh[v], u)
        ew[u][v] = get(ew[u], v, 0.0) + w
        ew[v][u] = get(ew[v], u, 0.0) + w
    end
    return BusADJ(neigh, ew)
end

# checagem rápida de conectividade (BFS) para uma zona após remover/adicionar 1 nó
function zone_connected_after_move(zone::Vector{Int}, z::Int, node::Int, add_to::Union{Nothing,Int},
                                   adj::BusADJ)
    # constrói lista de nós da zona z após remover 'node' (se zone[node]==z)
    nodes = Int[]
    for (i, zi) in enumerate(zone)
        if zi==z && i!=node
            push!(nodes, i)
        end
    end
    isempty(nodes) && return true
    # BFS dentro da zona z
    seen = Set{Int}([nodes[1]])
    Q = [nodes[1]]
    while !isempty(Q)
        u = popfirst!(Q)
        for v in adj.neigh[u]
            if v != node && zone[v]==z && !(v in seen)
                push!(seen, v); push!(Q, v)
            end
        end
    end
    return length(seen) == length(nodes)
end

# pmax por bus (para garantir 1 gerador útil por zona)
function pmax_per_bus(case::Dict{String,Any})
    g = get(case, "gen", Dict{String,Any}())
    isempty(g) && return zeros(length(case["bus"]["bus_i"]))
    buses = case["bus"]["bus_i"]; pos = Dict(buses[i]=>i for i in eachindex(buses))
    genbus = g["gen_bus"]; status = g["gen_status"]; pmax = g["pmax"]
    v = zeros(Float64, length(buses))
    for i in eachindex(genbus)
        status[i]==1 || continue
        v[pos[genbus[i]]] += max(pmax[i],0.0)
    end
    return v
end

# delta (#shared) e delta (força de corte local) ao mover 'u' para zona 'znew'
function local_deltas_move(u::Int, znew::Int, zone::Vector{Int}, adj::BusADJ)
    zold = zone[u]
    zold == znew && return (0, 0.0)

    function is_shared_after(b::Int, zonelocal::Vector{Int})
        zb = zonelocal[b]
        for v in adj.neigh[b]
            if zonelocal[v] != zb; return true; end
        end
        return false
    end

    # estados atuais
    shared_before_u = any(zone[v]!=zold for v in adj.neigh[u])
    cut_before = 0.0
    for v in adj.neigh[u]
        if zone[v] != zold
            cut_before += adj.edge_w[u][v]
        end
    end
    # contrib dos vizinhos para #shared antes
    shared_before_neigh = 0
    for v in adj.neigh[u]
        if any(zone[w]!=zone[v] for w in adj.neigh[v])
            shared_before_neigh += 0  # contaremos pelo delta abaixo (mudança local)
        end
    end

    # aplicar movimento virtual
    ztmp = copy(zone); ztmp[u] = znew

    shared_after_u = any(ztmp[v]!=znew for v in adj.neigh[u])
    cut_after = 0.0
    for v in adj.neigh[u]
        if ztmp[v] != znew
            cut_after += adj.edge_w[u][v]
        end
    end

    # variação do status de shared nos vizinhos causada por trocar z de u
    Δshared_neigh = 0
    for v in adj.neigh[u]
        before = any(zone[w] != zone[v] for w in adj.neigh[v])
        after  = any(ztmp[w] != ztmp[v] for w in adj.neigh[v])
        if before && !after; Δshared_neigh -= 1
        elseif !before && after; Δshared_neigh += 1
        end
    end

    Δshared_u = (shared_after_u ? 1 : 0) - (shared_before_u ? 1 : 0)
    Δshared = Δshared_u + Δshared_neigh
    Δcut = (cut_after - cut_before)
    return (Δshared, Δcut)
end

"""
    refine_bus_boundary!(case, zone; passes=4, max_zone_size=typemax(Int), verbose=true)

Refina rótulos `zone` reduzindo o nº de barramentos de fronteira.
Respeita conectividade, teto de tamanho por zona e pelo menos um gerador útil por zona.
"""
function refine_bus_boundary!(case::Dict{String,Any}, zone::Vector{Int};
    passes::Int=4, max_zone_size::Int=typemax(Int), verbose::Bool=true)

    buses = case["bus"]["bus_i"]
    k = maximum(zone)
    adj = build_adjacency_with_weights(case)

    # estruturas por zona
    sizeZ = [count(==(z), zone) for z in 1:k]
    pmax_bus = pmax_per_bus(case)
    pmaxZ = [0.0 for _ in 1:k]
    for (i, z) in enumerate(zone); pmaxZ[z] += pmax_bus[i]; end

    # zonas vizinhas de um nó (candidatas para movimento)
    zone_neighbors = function(u)
        Z = Set{Int}()
        for v in adj.neigh[u]
            if zone[v] != zone[u]; push!(Z, zone[v]); end
        end
        return collect(Z)
    end

    for pass in 1:passes
        moved_any = false
        total_shared, shared_mask = count_shared_buses(case, zone)
        verbose && @info "[refine pass $pass] shared=$total_shared"

        candidates = findall(shared_mask)
        isempty(candidates) && break
        shuffle!(candidates)

        for u in candidates
            z0 = zone[u]

            # se este nó é o único que garante pmax>0 na zona atual, não mova
            if pmax_bus[u] > 0 && pmaxZ[z0] - pmax_bus[u] <= 1e-9
                continue
            end

            nbr_zones = zone_neighbors(u)
            isempty(nbr_zones) && continue  # sem destino possível

            # mantém SEMPRE um NamedTuple
            best = (Δs = 0, Δc = 0.0, zbest = z0)

            for zt in nbr_zones
                # teto de tamanho e não esvaziar origem
                if sizeZ[zt] + 1 > max_zone_size || sizeZ[z0] - 1 <= 0
                    continue
                end
                # conectividade da origem após remover u
                if !zone_connected_after_move(zone, z0, u, nothing, adj)
                    continue
                end

                Δs, Δc = local_deltas_move(u, zt, zone, adj)

                # regra: reduzir #shared; em empate, reduzir força do corte
                if (Δs < best.Δs) || (Δs == best.Δs && Δc < best.Δc)
                    best = (Δs = Δs, Δc = Δc, zbest = zt)  # <-- mantém NamedTuple
                end
            end

            if best.zbest != z0
                newZ = best.zbest

                # checagem p/ não zerar pmax da origem
                if pmax_bus[u] > 0 && pmaxZ[z0] - pmax_bus[u] <= 1e-9
                    continue
                end

                # commit do movimento
                zone[u] = newZ
                sizeZ[z0] -= 1; sizeZ[newZ] += 1
                pmaxZ[z0] = max(pmaxZ[z0] - pmax_bus[u], 0.0)
                pmaxZ[newZ] += pmax_bus[u]
                moved_any = true
            end
        end

        # reparos (caso algum movimento tenha deixado algo frágil)
        fix_connectivity_and_slack!(case, zone, k)

        if !moved_any
            verbose && @info "[refine pass $pass] sem movimentos benéficos."
            break
        end
    end

    return zone
end

# ------------------------------------------------------------------
# Escolha automática de k (inclui alvo para shared buses)
# ------------------------------------------------------------------
function choose_k_bus(case; max_zone::Int=5000, kmin::Int=2, kmax::Int=2000)
    n = length(case["bus"]["bus_i"])
    k = max(kmin, ceil(Int, n / max_zone))
    return min(max(k, kmin), kmax)
end

# ------------------------------------------------------------------
# Pipeline BUS-ADMM
# ------------------------------------------------------------------
function partition_cases_metis_bus(indir::AbstractString=DEFAULT_DIR;
    ks::Union{Vector{Int},Symbol} = :auto,
    seed::Integer = 42,
    glob::AbstractString = "*.m",
    overwrite::Bool = false,
    balance_tol::Float64 = 0.10,
    max_zone_size::Int = 5000,
    target_shared_ratio::Float64 = 0.02,    # até ~2% dos buses na fronteira
    refine_passes::Int = 4,
)
    require_lirion()
    Random.seed!(seed)

    files = sort(Glob.glob(glob, indir))
    isempty(files) && @info "Nenhum arquivo $glob encontrado em $indir."

    for f in files
        base = basename(f)
        case = try
            load_case(f)
        catch e
            @warn "Falha ao carregar $base: $e"; continue
        end
        has_multiple_zones(case) && (@info "PULANDO (já tem >1 zona): $base"; continue)

        # decidir k(s)
        local ks_list::Vector{Int}
        if ks === :auto
            ks_list = Int[ choose_k_bus(case; max_zone=max_zone_size) ]
        else
            ks_list = collect(Int, ks)
        end

        for k in ks_list
            outdir = joinpath(indir, "case $(k) zones (metis-bus)")
            isdir(outdir) || mkpath(outdir)
            outpath = joinpath(outdir, base)
            if isfile(outpath) && !overwrite
                @info "PULANDO (destino existe): $outpath"; continue
            end

            # 1) METIS bruto
            z = metis_partition(case; k=k, balance_tol=balance_tol)

            # 2) Refinamento focado BUS-ADMM
            z = refine_bus_boundary!(case, z; passes=refine_passes, max_zone_size=max_zone_size)

            # 3) Métricas
            counts = Dict{Int,Int}(); for zz in z; counts[zz] = get(counts, zz, 0) + 1; end
            total_buses = length(z)
            max_block = maximum(values(counts))
            n_shared, _ = count_shared_buses(case, z)
            shared_ratio = n_shared / max(1,total_buses)
            cutS = edge_cut_strength(case, z)

            @info "Resumo zonas (BUS-ADMM) $base — k=$k, total=$total_buses"
            for zz in sort(collect(keys(counts)))
                @info "  zone $zz: $(counts[zz]) buses"
            end
            @info @sprintf("  shared buses: %d (%.2f%%) | cut strength≈ %.2f",
                           n_shared, 100*shared_ratio, cutS)

            # 4) Reforço de teto / shared_ratio (aumenta k se necessário)
            tries = 0
            while ks === :auto && (max_block > max_zone_size || shared_ratio > target_shared_ratio) && tries < 6
                k = max(k + 1, ceil(Int, total_buses / max_zone_size))
                @info "→ Repetindo com k=$k (teto=$max_zone_size, alvo shared≤$(round(100*target_shared_ratio,digits=1))%)"
                z = metis_partition(case; k=k, balance_tol=balance_tol)
                z = refine_bus_boundary!(case, z; passes=refine_passes, max_zone_size=max_zone_size)
                counts = Dict{Int,Int}(); for zz in z; counts[zz] = get(counts, zz, 0) + 1; end
                max_block = maximum(values(counts))
                n_shared, _ = count_shared_buses(case, z)
                shared_ratio = n_shared / max(1,total_buses)
                cutS = edge_cut_strength(case, z)
                @info @sprintf("  → shared: %d (%.2f%%) | maior bloco=%d | cut≈%.2f",
                               n_shared, 100*shared_ratio, max_block, cutS)
                tries += 1
            end

            apply_zones!(case, z)
            write_matpower_with_zones(outpath, case)
            @info "OK: $(base) → $(outpath)"
        end
    end
    @info "Concluído (METIS BUS-ADMM)."
    return nothing
end

# --------------------------------------------
# Execução direta (edite os caminhos/flags se quiser)
# --------------------------------------------
# partition_cases_metis_bus(
#     raw"C:\\Users\\franc\\Lirion.jl\\test\\zones";
#     ks = :auto,                     # ou ks=[25,50,100]
#     max_zone_size = 1500,           # teto desejado por zona
#     target_shared_ratio = 0.02,     # ≤2% de nós de fronteira (ajuste p/ sistemas grandes)
#     balance_tol = 0.10,             # ±10%
#     refine_passes = 4,
#     glob = "*.m",
#     overwrite = false,
# )

partition_cases_metis_bus(
    raw"/home/mateus-kaiss/Projects/Lirion.jl/test/zones";
    ks = :auto,
    max_zone_size = 200,        # força k≈ceil(2000/450)=5
    balance_tol = 0.08,
    refine_passes = 6,
    target_shared_ratio = 0.012,
    glob = "*.m",
    overwrite = false,
)
