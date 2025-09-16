using Printf
using LinearAlgebra
using SparseArrays
using Graphs
using Random
using Arpack
using Glob 
using PowerOptimization
using Statistics

const DEFAULT_DIR = normpath(@__DIR__, "..", "test", "matpower_cases")

function require_poweroptimization()
    try
        @eval using PowerOptimization
    catch e
        error("""
        Não consegui carregar `PowerOptimization`.
        Certifique-se de:
          1) Estar no ambiente certo (julia --project=.)
          2) Ter `PowerOptimization.jl` no ambiente (ex: Pkg.develop(path=".../PowerOptimization.jl"))
        Erro: $e
        """)
    end
end

load_case(path::String) = PowerOptimization.load_input_data(path)

function has_multiple_zones(case::Dict{String,Any})
    if haskey(case, "bus") && haskey(case["bus"], "zone")
        z = case["bus"]["zone"]
        return length(unique(z)) > 1
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

        # BRANCH  (angmin/angmax → GRAUS!)
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

# --------------------------------------------
# Construção do grafo
# --------------------------------------------
"Retorna (g, edge_weights, vertex_weights, present, posmap)."
function build_bus_graph(case::Dict{String,Any}; edge_weight::Symbol=:susceptance)
    f = case["branch"]["f_bus"]; t = case["branch"]["t_bus"]
    x = case["branch"]["br_x"];  r = case["branch"]["br_r"]
    on = case["branch"]["br_status"]
    tap = get(case["branch"], "tap", fill(1.0, length(f)))
    shift = get(case["branch"], "shift", fill(0.0, length(f)))
    ratea = get(case["branch"], "rate_a", fill(0.0, length(f)))

    buses = case["bus"]["bus_i"]
    maxbus = maximum(buses)
    present = falses(maxbus)
    for b in buses; present[b] = true; end

    g = SimpleGraph(maxbus)
    raw_w = Float64[]
    for i in eachindex(f)
        on[i] == 1 || continue
        u, v = f[i], t[i]
        present[u] && present[v] || continue
        add_edge!(g, u, v)

        w = edge_weight === :susceptance ? (abs(x[i]) > 1e-10 ? 1/abs(x[i]) : 1e6) :
            edge_weight === :impedance   ? 1/(x[i]^2 + r[i]^2 + 1e-12) :
            1.0

        if abs(tap[i] - 1.0) > 1e-6 || abs(shift[i]) > 1e-6
            w *= 50.0
        end
        if ratea[i] > 0
            w *= (1 + ratea[i])
        end

        push!(raw_w, w)
    end
    isempty(raw_w) && (raw_w = [1.0])
    maxw = maximum(raw_w)
    ew = round.(Int, 1 .+ 10_000 .* (raw_w ./ maxw))

    vw = ones(Float64, maxbus)
    pd_by_id = Dict{Int,Float64}()
    if haskey(case["bus"], "pd")
        for i in eachindex(buses)
            pd_by_id[buses[i]] = max(case["bus"]["pd"][i], 0.0)
        end
    end
    gen = get(case, "gen", Dict{String,Any}())
    if !isempty(gen)
        for i in eachindex(gen["gen_bus"])
            b = gen["gen_bus"][i]
            vw[b] += max(gen["pmax"][i], 0.0)
        end
    end
    for (b,pd) in pd_by_id
        vw[b] += pd
    end
    for b in 1:maxbus
        present[b] || (vw[b] = 0)
    end
    sum(vw) == 0 && (vw[present] .= 1)

    posmap = Dict(buses[i]=>i for i in eachindex(buses))
    return g, ew, vw, present, posmap
end

# --- pesos de corte: punir cortar trafo/phase shifter/linhas fortes ---
function _build_cut_weights(case::Dict{String,Any}, g::SimpleGraph)
    f = case["branch"]["f_bus"]; t = case["branch"]["t_bus"]
    x = case["branch"]["br_x"];  r = case["branch"]["br_r"]
    on = case["branch"]["br_status"]
    tap   = get(case["branch"], "tap",   fill(1.0, length(f)))
    shift = get(case["branch"], "shift", fill(0.0, length(f)))
    ratea = get(case["branch"], "rate_a", fill(0.0, length(f)))

    cutw = Dict{Tuple{Int,Int},Float64}()
    for i in eachindex(f)
        on[i]==1 || continue
        u, v = f[i], t[i]
        has_edge(g, u, v) || continue
        w = (abs(x[i]) > 1e-10) ? 1/abs(x[i]) : 1e6      # forte = caro cortar
        if abs(tap[i]-1.0) > 1e-6 || abs(shift[i]) > 1e-6
            w *= 100.0                                   # trafo/phase shifter MUITO caro cortar
        end
        if ratea[i] > 0
            w *= (1 + 0.01*ratea[i])                     # capacidade alta => mais caro
        end
        key = u < v ? (u,v) : (v,u)
        cutw[key] = get(cutw, key, 0.0) + w
    end
    return cutw
end

_edge_key(u::Int,v::Int) = u < v ? (u,v) : (v,u)

# --------------------------------------------
# Particionamento orientado a ADMM (min-cut ponderado + mínimo de 15%)
# --------------------------------------------

# Fiedler seguro
function _fiedler_vector(L::SparseMatrixCSC{Float64,Int})
    n = size(L,1)
    if n <= 2
        return nothing
    end
    if n <= 200
        vals, vecs = eigen(Symmetric(Matrix(L)))
        return length(vals) >= 2 ? vecs[:,2] : nothing
    else
        try
            vals, vecs = Arpack.eigs(L; nev=2, which=:SR)
            return size(vecs,2) >= 2 ? vecs[:,2] : nothing
        catch
            return nothing
        end
    end
end

# Laplaciano do subgrafo induzido por "nodes" (IDs absolutos)
function _laplacian_of_induced(g::SimpleGraph, nodes::Vector{Int})
    mappos = Dict(nodes[i]=>i for i in eachindex(nodes))
    n = length(nodes)
    I = Int[]; J = Int[]; V = Float64[]
    deg = zeros(Float64, n)
    for (i, u) in enumerate(nodes)
        for v in neighbors(g, u)
            haskey(mappos, v) || continue
            j = mappos[v]
            if i != j
                push!(I, i); push!(J, j); push!(V, -1.0)
                deg[i] += 1.0
            end
        end
    end
    for i in 1:n
        push!(I, i); push!(J, i); push!(V, deg[i])
    end
    return sparse(I, J, V, n, n)
end

"""
    _balanced_bisection_cut!(g, nodes, cutw; min_zone_ratio=0.15)

Divide `nodes` em (left, right) minimizando custo de corte ponderado (`cutw`)
e garantindo que ambas as partes tenham ao menos `min_zone_ratio * length(nodes)`.
Retorna (left, right). Se não houver split viável, retorna (nodes, Int[]).
"""
function _balanced_bisection_cut!(g::SimpleGraph, nodes::Vector{Int}, cutw::Dict{Tuple{Int,Int},Float64};
                                  min_zone_ratio::Float64=0.15)
    n = length(nodes)
    min_size = max(1, ceil(Int, n * min_zone_ratio))
    if n < 2*min_size
        return nodes, Int[]
    end

    # ordenação por Fiedler quando possível (fallback: grau)
    L = _laplacian_of_induced(g, nodes)
    f = _fiedler_vector(L)
    order = if f === nothing || all(abs.(f .- f[1]) .< 1e-12)
        degs = [degree(g, u) for u in nodes]
        sortperm(degs; rev=true)
    else
        sortperm(f)
    end

    s_lo = min_size
    s_hi = n - min_size
    step = max(1, fld(n, 100))  # ~100 amostras
    best_s = 0
    best_cost = Inf

    function cut_cost(s::Int)
        leftset = Set(nodes[order[1:s]])
        cost = 0.0
        @inbounds for i in (s+1):n
            u = nodes[order[i]]
            for v in neighbors(g, u)
                if v in leftset
                    cost += get(cutw, _edge_key(u,v), 1.0)
                end
            end
        end
        return cost
    end

    for s in s_lo:step:s_hi
        c = cut_cost(s)
        if c < best_cost
            best_cost = c; best_s = s
        end
    end
    if best_s == 0
        return nodes, Int[]
    end
    # refino local fino
    for s in max(s_lo, best_s-20):min(s_hi, best_s+20)
        c = cut_cost(s)
        if c < best_cost
            best_cost = c; best_s = s
        end
    end

    left  = [nodes[order[i]] for i in 1:best_s]
    right = [nodes[order[i]] for i in (best_s+1):n]
    return left, right
end

# Rebalance leve por tamanho (move nós de fronteira)
function _rebalance_by_size!(g::SimpleGraph, labels::Dict{Int,Int}, k::Int; balance_tol::Float64=0.10)
    nodes = collect(keys(labels))
    n = length(nodes)
    target = round(Int, n / k)
    minZ = max(1, round(Int, target*(1 - balance_tol)))
    maxZ = round(Int, target*(1 + balance_tol))

    zones = Dict(z=>Int[] for z in 1:k)
    for u in nodes
        push!(zones[labels[u]], u)
    end

    over = [z for z in 1:k if length(zones[z]) > maxZ]
    under = [z for z in 1:k if length(zones[z]) < minZ]
    if isempty(over) && isempty(under)
        return
    end

    for zbig in over
        while length(zones[zbig]) > maxZ
            moved = false
            for u in zones[zbig]
                nbr = neighbors(g, u)
                nbrzones = unique([labels[v] for v in nbr])
                for zsmall in under
                    if zsmall in nbrzones && length(zones[zsmall]) < minZ
                        labels[u] = zsmall
                        push!(zones[zsmall], u)
                        filter!(x -> x != u, zones[zbig])
                        moved = true
                        break
                    end
                end
                if moved; break; end
            end
            if !moved
                degs = [(degree(g,u),u) for u in zones[zbig]]
                sort!(degs)            # menor grau primeiro
                u = last(first(degs))  # nó de menor grau
                idx_small = findfirst(z -> length(zones[z]) < minZ, under)
                idx_small === nothing && break
                zsmall = under[idx_small]
                labels[u] = zsmall
                push!(zones[zsmall], u)
                filter!(x -> x != u, zones[zbig])
            end
            over  = [z for z in 1:k if length(zones[z]) > maxZ]
            under = [z for z in 1:k if length(zones[z]) < minZ]
            isempty(over) && break
        end
    end
end

"""
    partition_spectral(case, g, k, present, posmap;
                       min_zone_ratio=0.15, balance_tol=0.10)

Particionamento k-way orientado a **minimizar corte ponderado** e garantindo
tamanho mínimo por zona de `min_zone_ratio` do tamanho do bloco no momento da divisão.
"""
function partition_spectral(case, g::SimpleGraph, k::Int, present::BitVector, posmap::Dict{Int,Int};
                            min_zone_ratio::Float64=0.15, balance_tol::Float64=0.10)

    all_buses = sort(collect(keys(posmap)))
    cutw = _build_cut_weights(case, g)

    # fila começa com um único bloco (todas as barras)
    queue = [all_buses]
    buckets = Vector{Vector{Int}}()

    # --- LOOP ROBUSTO: não tente argmax em fila vazia
    while (length(buckets) + length(queue) < k) && !isempty(queue)
        # escolha o maior bloco remanescente para tentar dividir
        idx = argmax(length.(queue))
        nodes = queue[idx]
        deleteat!(queue, idx)

        left, right = _balanced_bisection_cut!(g, nodes, cutw; min_zone_ratio=min_zone_ratio)
        if isempty(right)
            # não dá pra dividir respeitando o mínimo -> fixa em buckets
            push!(buckets, nodes)
        else
            # continua quebrando o MAIOR pedaço primeiro
            if length(left) >= length(right)
                push!(queue, left); push!(buckets, right)
            else
                push!(queue, right); push!(buckets, left)
            end
        end
    end

    # Se ainda não atingiu k e a fila acabou, aplica "force-split" seguro
    # (relaxa o min_zone_ratio de forma adaptativa e aceita singletons se precisar)
    while length(buckets) + length(queue) < k
        isempty(buckets) && break
        # pega o maior bucket para fatiar
        big_idx = argmax(length.(buckets))
        nodes = splice!(buckets, big_idx)

        n = length(nodes)
        if n <= 1
            # já é singleton, reempilha e sai do laço se só restam singletons
            push!(buckets, nodes)
            # se todos os buckets são tamanho 1, não há mais o que fazer
            all_singletons = all(length(b) == 1 for b in buckets)
            all_singletons && break
            continue
        end

        # relaxa a razão mínima conforme o bloco encolhe
        # (quanto menor o bloco, menor a exigência mínima)
        local_min_ratio = max(0.0, min(min_zone_ratio, 0.5*sqrt(1/n)))

        left, right = _balanced_bisection_cut!(g, nodes, cutw; min_zone_ratio=local_min_ratio)
        if isempty(right)
            # split falhou mesmo relaxando -> faz um split trivial 1 | (n-1)
            # preservando conectividade quando possível (move um nó de menor grau)
            if n >= 2
                degs = [(degree(g,u), u) for u in nodes]
                sort!(degs)                 # menor grau primeiro
                u = last(first(degs))       # nó de menor grau
                left = [u]
                right = [v for v in nodes if v != u]
            else
                push!(buckets, nodes)
                break
            end
        end

        # empurrar de volta: continuar quebrando o maior
        if length(left) >= length(right)
            push!(queue, left); push!(buckets, right)
        else
            push!(queue, right); push!(buckets, left)
        end
    end

    # Junta fila e buckets e, se passou de k por arredondamentos, combina menores
    append!(buckets, queue)
    while length(buckets) > k
        sort!(buckets, by=length)           # dois menores
        a = popfirst!(buckets); b = popfirst!(buckets)
        push!(buckets, vcat(a,b))
    end

    # Rotulagem
    label_of = Dict{Int,Int}()
    for (z, nodes) in enumerate(buckets)
        for u in nodes
            label_of[u] = z
        end
    end

    # Rebalance leve
    _rebalance_by_size!(g, label_of, k; balance_tol=balance_tol)

    # Monta vetor zone alinhado ao case
    buses = case["bus"]["bus_i"]
    zone = fill(1, length(buses))
    for (i,b) in enumerate(buses)
        zone[i] = get(label_of, b, 1)
    end

    # Conectividade + garantir geração/slack por zona
    fix_connectivity_and_slack!(case, zone, k)
    return zone
end

# --------------------------------------------
# Pós-processos (conectividade + geração/slack)
# --------------------------------------------
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

    for z in 1:k
        nodes = [buses[i] for i in eachindex(buses) if zone[i]==z]
        isempty(nodes) && continue
        seen = Set{Int}()
        comps = Vector{Vector{Int}}()
        for u in nodes
            u in seen && continue
            comp = Int[]
            Q = [u]; push!(seen,u)
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
            if isempty(border)
                continue
            else
                counts = Dict{Int,Int}()
                for (_,v) in border
                    counts[zone[pos[v]]] = get(counts, zone[pos[v]], 0) + 1
                end
                zbest = let ks = collect(keys(counts)); ks[argmax(collect(values(counts)))] end
                for u in comp
                    zone[pos[u]] = zbest
                end
            end
        end
    end

    gen = get(case, "gen", Dict{String,Any}())
    if !isempty(gen)
        genbus = gen["gen_bus"]
        status = gen["gen_status"]
        pmax   = gen["pmax"]
        pos = Dict(case["bus"]["bus_i"][i]=>i for i in eachindex(case["bus"]["bus_i"]))

        pmax_by_zone = Dict(z=>0.0 for z in 1:k)
        for i in eachindex(genbus)
            status[i]==1 || continue
            b  = genbus[i]
            z  = zone[pos[b]]
            pM = max(pmax[i], 0.0)
            pmax_by_zone[z] = get(pmax_by_zone, z, 0.0) + pM
        end

        f = case["branch"]["f_bus"]; t = case["branch"]["t_bus"]; on = case["branch"]["br_status"]
        neigh = Dict{Int,Set{Int}}()
        for i in eachindex(f)
            on[i]==1 || continue
            zf = zone[pos[f[i]]]; zt = zone[pos[t[i]]]
            if zf != zt
                push!(get!(neigh, zf, Set{Int}()), zt)
                push!(get!(neigh, zt, Set{Int}()), zf)
            end
        end

        gens_in_zone = Dict(z=>Int[] for z in 1:k)
        for i in eachindex(genbus)
            status[i]==1 || continue
            if pmax[i] > 1e-6
                bz = zone[pos[genbus[i]]]
                push!(gens_in_zone[bz], i)
            end
        end

        for z in 1:k
            if get(pmax_by_zone, z, 0.0) > 1e-9
                continue  # já tem pmax
            end

            moved = false
            for zv in sort(collect(get(neigh, z, Set{Int}())))
                if !isempty(get(gens_in_zone, zv, Int[]))
                    gi = pop!(gens_in_zone[zv]) 
                    b  = genbus[gi]
                    zone[pos[b]] = z                
                    pmax_by_zone[z]  = get(pmax_by_zone, z, 0.0) + max(pmax[gi],0.0)
                    pmax_by_zone[zv] = max(get(pmax_by_zone, zv, 0.0) - max(pmax[gi],0.0), 0.0)
                    push!(gens_in_zone[z], gi)
                    moved = true
                    break
                end
            end

            if !moved
                for zz in 1:k
                    zz == z && continue
                    if !isempty(get(gens_in_zone, zz, Int[]))
                        gi = pop!(gens_in_zone[zz])
                        b  = genbus[gi]
                        zone[pos[b]] = z
                        pmax_by_zone[z]  += max(pmax[gi],0.0)
                        pmax_by_zone[zz]  = max(pmax_by_zone[zz] - max(pmax[gi],0.0), 0.0)
                        push!(gens_in_zone[z], gi)
                        moved = true
                        break
                    end
                end
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

# --------------------------------------------
# Batch runner (função pública)
# --------------------------------------------
function partition_cases(indir::AbstractString=DEFAULT_DIR;
    ks::Vector{Int} = [2,3,4,5],
    seed::Integer = 42,
    glob::AbstractString = "*.m",
    overwrite::Bool = false,
    balance_tol::Float64 = 0.10,
    min_zone_ratio::Float64 = 0.15,
)
    require_poweroptimization()

    files = sort(Glob.glob(glob, indir))
    isempty(files) && @info "Nenhum arquivo $glob encontrado em $indir."

    for k in ks
        outdir = joinpath(indir, "case $(k) zones")
        isdir(outdir) || mkpath(outdir)
        @info "==== k=$k → salvando em $(outdir)"
        for f in files
            base = basename(f)
            outpath = joinpath(outdir, base)

            case = try
                load_case(f)
            catch e
                @warn "Falha ao carregar $base: $e"
                continue
            end

            # pular fontes que já têm >1 zona
            if has_multiple_zones(case)
                @info "PULANDO (fonte já tem >1 zona): $(base)"
                continue
            end

            if isfile(outpath) && !overwrite
                @info "PULANDO (destino já existe): $(outpath)"
                continue
            end

            Random.seed!(seed)
            g, _, _, present, posmap = build_bus_graph(case)
            z = partition_spectral(case, g, k, present, posmap;
                                   min_zone_ratio=min_zone_ratio,
                                   balance_tol=balance_tol)
            apply_zones!(case, z)

            # resumo de tamanhos
            counts = Dict{Int,Int}()
            for zz in z
                counts[zz] = get(counts, zz, 0) + 1
            end
            total_buses = length(z)
            @info "Resumo zonas ($base): total_buses = $total_buses"
            for zz in sort(collect(keys(counts)))
                @info "  zone $zz: $(counts[zz]) buses"
            end

            write_matpower_with_zones(outpath, case)
            @info "OK: $(base) → $(outpath)"
        end
    end
    @info "Concluído."
    return nothing
end

# --------------------------------------------
# CLI opcional
# --------------------------------------------
function parse_args()
    args = Dict{String,Any}()
    args["--in"] = DEFAULT_DIR
    args["--zones"] = Int[2]
    args["--seed"] = 42
    args["--glob"] = "*.m"
    args["--overwrite"] = false
    args["--balance_tol"] = 0.10
    args["--min_zone_ratio"] = 0.15

    i = 1
    while i <= length(ARGS)
        a = ARGS[i]
        if a == "--in";       i+=1; args["--in"] = ARGS[i]
        elseif a == "--zones"
            z = Int[]; i+=1
            while i <= length(ARGS) && !startswith(ARGS[i], "--")
                push!(z, parse(Int, ARGS[i])); i+=1
            end
            i-=1; args["--zones"] = z
        elseif a == "--seed";   i+=1; args["--seed"] = parse(Int, ARGS[i])
        elseif a == "--glob";   i+=1; args["--glob"] = ARGS[i]
        elseif a == "--overwrite"; args["--overwrite"] = true
        elseif a == "--balance_tol"; i+=1; args["--balance_tol"] = parse(Float64, ARGS[i])
        elseif a == "--min_zone_ratio"; i+=1; args["--min_zone_ratio"] = parse(Float64, ARGS[i])
        else
            @warn "Argumento desconhecido: $a"
        end
        i+=1
    end
    return args
end

function main()
    a = parse_args()
    partition_cases(
        String(a["--in"]);
        ks = Vector{Int}(a["--zones"]),
        seed = Int(a["--seed"]),
        glob = String(a["--glob"]),
        overwrite = Bool(a["--overwrite"]),
        balance_tol = Float64(a["--balance_tol"]),
        min_zone_ratio = Float64(a["--min_zone_ratio"]),
    )
end

# Para rodar direto via include:
partition_cases(
    raw"C:\Users\franc\PowerOptimization.jl\test\zones";
    ks=[2,3,4,5,25,50,100,150,200,300],
    seed=42,
    glob="*.m",
    overwrite=false,
    balance_tol=0.10,
    min_zone_ratio=0.15,   # nenhuma zona < 15% do total
)
