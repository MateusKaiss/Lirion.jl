############################################################
# Particionamento uniforme (uma partição → aplica em todos)
# - Estratégias: :bus (METIS + refinamento fronteira) | :edgecut (METIS puro)
############################################################
using Printf, LinearAlgebra, SparseArrays, Graphs, Random, Glob
using Statistics
using Metis
using SimpleWeightedGraphs

# --- Lirion (carregar do ambiente atual) ---
try
    @eval using Lirion
catch e
    error("""
    Não consegui carregar `Lirion`.
    Certifique-se de:
      1) Rodar com o --project correto
      2) Ter Lirion.jl no ambiente (Pkg.develop/instantiate)
    Erro: $e
    """)
end

# ---------------------------
# utilidades básicas
# ---------------------------
load_case(path::String) = Lirion.load_input_data(path)

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
        @assert haskey(b, "zone") "Coluna 'zone' não encontrada em bus."
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

        # BRANCH  (angmin/angmax em graus)
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

# ---------------------------
# Métricas auxiliares
# ---------------------------
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
    return count(shared)
end

# ---------------------------
# Helpers p/ rótulos (garantir ≥2 zonas e compactar)
# ---------------------------
# ---------------------------
# Helpers p/ rótulos (garantir ≥2 zonas e compactar)
# ---------------------------
function _relabel_compact!(Z::Vector{Int})
    labs = sort(unique(Z))
    m = length(labs)
    if m == 0
        return 0
    end
    map = Dict(labs[i] => i for i in eachindex(labs))
    @inbounds for i in eachindex(Z)
        Z[i] = map[Z[i]]
    end
    return m
end

function _ensure_k_labels!(case::Dict{String,Any}, Z::Vector{Int}, K::Int)
    n = length(Z)
    target = min(K, n)
    m = _relabel_compact!(Z)
    m >= target && return Z

    # maior bloco atual
    counts = Dict(z => count(==(z), Z) for z in unique(Z))
    zbig   = argmax(counts)  # rótulo com maior cardinalidade

    # preferir mover um bus com gerador a partir de zbig
    buses = case["bus"]["bus_i"]
    pos   = Dict(buses[i]=>i for i in eachindex(buses))
    gen   = get(case, "gen", Dict{String,Any}())
    genbus = get(gen, "gen_bus", Int[])

    # CORREÇÃO: usar haskey(pos, b) (e não b in pos)
    candidates = [pos[b] for b in genbus if haskey(pos, b) && Z[pos[b]] == zbig]

    # se não houver gerador em zbig, qualquer nó do bloco
    if isempty(candidates)
        candidates = [i for i in eachindex(Z) if Z[i] == zbig]
    end
    isempty(candidates) && return Z  # nada a fazer

    need = target - m
    # garanta que haja 'need' candidatos (repete do bloco maior se preciso)
    if length(candidates) < need
        extra = [i for i in eachindex(Z) if Z[i] == zbig && !(i in candidates)]
        append!(candidates, extra)
    end
    need = min(need, length(candidates))

    for j in 1:need
        idx = candidates[j]
        Z[idx] = m + j  # cria novo rótulo (m+1, m+2, ...)
    end
    _relabel_compact!(Z)
    return Z
end

# ---------------------------
# Grafo ponderado e METIS comum
# ---------------------------
function _build_weighted_graph(case::Dict{String,Any})
    buses  = case["bus"]["bus_i"]
    present = Set(buses)
    f   = case["branch"]["f_bus"]; t = case["branch"]["t_bus"]
    x   = case["branch"]["br_x"];  on  = case["branch"]["br_status"]
    tap   = get(case["branch"], "tap",   fill(1.0, length(f)))
    shift = get(case["branch"], "shift", fill(0.0, length(f)))
    ratea = get(case["branch"], "rate_a", fill(0.0, length(f)))

    glob2loc = Dict{Int,Int}(); n = length(buses)
    for (i,b) in enumerate(buses) glob2loc[b] = i end

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

    G = SimpleWeightedGraph(length(buses))
    for ((u,v), w) in wdict
        add_edge!(G, u, v, w)
    end
    return G, buses, glob2loc
end

function metis_partition_edgecut(case::Dict{String,Any}; k::Int)
    G, buses, glob2loc = _build_weighted_graph(case)
    parts = Metis.partition(G, k; alg = :KWAY)   # rótulos em 1..k (pode "degenerar")
    Z = fill(1, length(buses))
    @inbounds for (i, bid) in enumerate(buses)
        Z[i] = parts[glob2loc[bid]]
    end
    _ensure_k_labels!(case, Z, k)   # garante pelo menos min(k,n) zonas não-vazias
    return Z
end

# ---------------------------
# Refinamento p/ BUS-ADMM (reduz fronteira)
# ---------------------------
function fix_connectivity!(case::Dict{String,Any}, zone::Vector{Int})
    buses = case["bus"]["bus_i"]
    f = case["branch"]["f_bus"]; t = case["branch"]["t_bus"]; on = case["branch"]["br_status"]
    neigh = Dict{Int,Vector{Int}}()
    for i in eachindex(f)
        on[i]==1 || continue
        push!(get!(neigh, f[i], Int[]), t[i])
        push!(get!(neigh, t[i], Int[]), f[i])
    end
    pos = Dict(buses[i]=>i for i in eachindex(buses))
    k = maximum(zone)
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
            ks = collect(keys(counts))
            vs = [counts[k] for k in ks]
            zbest = ks[argmax(vs)]
            for u in comp
                zone[pos[u]] = zbest
            end
        end
    end
    return zone
end

function metis_partition_bus(case::Dict{String,Any}; k::Int, refine_passes::Int=3)
    Z = metis_partition_edgecut(case; k=k)
    Z = _ensure_k_labels!(case, Z, k)  # reforço extra
    # refinamento leve: algumas passagens de conectividade e redução de shared
    for _ in 1:refine_passes
        fix_connectivity!(case, Z)
    end
    return Z
end

# ---------------------------
# Escolha automática de k
# ---------------------------
choose_k(n_buses::Int; max_zone::Int=5000, kmin::Int=2, kmax::Int=2000) =
    min(max(max(kmin, ceil(Int, n_buses / max_zone)), kmin), kmax)

# ---------------------------
# Aplicação uniforme na subpasta
# ---------------------------
"""
uniform_partition_subfolder(indir;
    strategy = :bus,            # :bus | :edgecut
    k = :auto,                  # :auto | Int
    seed = 42,
    refine_passes = 4,          # só afeta :bus
    max_zone_size = 5000,       # para k=:auto → k≈ceil(n/max_zone_size)
    overwrite = false,
) -> outdir

- Lê todos os .m de `indir`
- Decide um arquivo referência (prioriza *_lmbd.m)
- Cria partição **uma vez** no ref e reaplica (mesmo bus_i→zone) nos demais
- Salva em `indir/partitions_uniform_(strategy)_kK/arquivo.m`
"""
function uniform_partition_subfolder(indir::AbstractString;
    strategy::Symbol = :bus,
    k::Union{Int,Symbol} = :auto,
    seed::Integer = 42,
    refine_passes::Int = 4,
    max_zone_size::Int = 5000,
    overwrite::Bool = false,
)
    @assert isdir(indir) "Diretório não existe: $indir"
    Random.seed!(seed)

    files = sort(Glob.glob("*.m", indir))
    isempty(files) && error("Nenhum .m encontrado em $indir")

    # referência: prioriza *_lmbd.m (STD). Se não houver, pega o primeiro.
    ref = begin
        cand = filter(f -> occursin(r"_lmbd\.m$", basename(f)), files)
        isempty(cand) ? first(files) : first(cand)
    end
    ref_case = load_case(ref)
    n = length(ref_case["bus"]["bus_i"])

    # decide k
    K = k === :auto ? choose_k(n; max_zone=max_zone_size) : Int(k)
    @info "Referência: $(basename(ref)) | n_buses=$n | k=$K | estratégia=$(String(strategy))"

    # constrói partição no ref
    Z = if strategy == :bus
        metis_partition_bus(ref_case; k=K, refine_passes=refine_passes)
    elseif strategy == :edgecut
        metis_partition_edgecut(ref_case; k=K)
    else
        error("strategy deve ser :bus ou :edgecut")
    end
    # métricas da ref
    shared_ref = count_shared_buses(ref_case, Z)
    @info @sprintf("Ref: shared buses = %d (%.2f%%)", shared_ref, 100*shared_ref/max(n,1))

    # mapa bus->zone da referência (por bus_i)
    buses_ref = ref_case["bus"]["bus_i"]
    zmap = Dict(buses_ref[i] => Z[i] for i in eachindex(buses_ref))

    # pasta de saída
    outdir = joinpath(indir, @sprintf("partitions_uniform_%s_k%d", String(strategy), K))
    isdir(outdir) || mkpath(outdir)

    # aplica a TODOS os .m da subpasta
    for f in files
        base = basename(f)
        outpath = joinpath(outdir, base)
        if isfile(outpath) && !overwrite
            @info "PULANDO (existe): $outpath"
            continue
        end

        c = load_case(f)
        buses = c["bus"]["bus_i"]

        # valida conjunto de buses
        if Set(buses) != Set(buses_ref)
            @warn "Topo/buses divergente de $(basename(ref)) — mantendo comparação fiel, vou PULAR." file=base
            continue
        end

        # aplica por bus_i
        zone = [ zmap[b] for b in buses ]
        c["bus"]["zone"] = zone

        # salvar
        write_matpower_with_zones(outpath, c)

        # métricas por arquivo
        shared = count_shared_buses(c, zone)
        @info @sprintf("OK: %-30s → zones=%d | shared=%d (%.2f%%)",
                       base, length(unique(zone)), shared, 100*shared/length(zone))
    end

    @info "Concluído. Saída em: $outdir"
    return outdir
end

# ---------------------------
# EXEMPLOS DE USO
# ---------------------------
# 1) Particionar 03-bus com estratégia BUS (refinamento), k automático (~ n/max_zone_size)
# uniform_partition_subfolder("/home/mateus-kaiss/Projects/experiments/03-bus";
#     strategy=:bus, k=:auto, max_zone_size=50, refine_passes=4, overwrite=false)

# 2) Mesma coisa porém usando METIS “puro” (edge-cut balanceado)
# uniform_partition_subfolder("/home/mateus-kaiss/Projects/experiments/03-bus";
#     strategy=:edgecut, k=3, overwrite=true)

# 3 buses cases (forçando k=2)
uniform_partition_subfolder("/home/mateus-kaiss/Projects/experiments/03-bus";
    strategy=:bus, k=2, refine_passes=3, overwrite=true)
