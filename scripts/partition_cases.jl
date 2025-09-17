############################################################
# Particionador por METIS (edge-cut balanceado) — drop-in
############################################################
using Printf, LinearAlgebra, SparseArrays, Graphs, Random, Glob
using PowerOptimization, Statistics
using Metis   # pacote Metis.jl (módulo Metis)

using SimpleWeightedGraphs  # <- garanta que está no topo do arquivo
const DEFAULT_DIR = normpath(@__DIR__, "..", "test", "matpower_cases")

# ---------------------------
# utilidades já existentes
# ---------------------------
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

# ------------------------------------------------------------------
# 1) Grafo + pesos (iguais aos seus) em estrutura para o METIS
# ------------------------------------------------------------------
function build_metis_inputs(case::Dict{String,Any})
    f = case["branch"]["f_bus"]; t = case["branch"]["t_bus"]
    x = case["branch"]["br_x"];  r = case["branch"]["br_r"]
    on = case["branch"]["br_status"]
    tap   = get(case["branch"], "tap",   fill(1.0, length(f)))
    shift = get(case["branch"], "shift", fill(0.0, length(f)))
    ratea = get(case["branch"], "rate_a", fill(0.0, length(f)))

    buses = case["bus"]["bus_i"]
    maxbus = maximum(buses)
    present = falses(maxbus)
    for b in buses; present[b] = true; end
    posmap = Dict(buses[i]=>i for i in eachindex(buses))

    g = SimpleGraph(maxbus)
    adj_dict = Dict{Tuple{Int,Int}, Int32}()

    raw = Float64[]
    for i in eachindex(f)
        on[i]==1 || continue
        u, v = f[i], t[i]
        present[u] && present[v] || continue
        add_edge!(g, u, v)

        w = (abs(x[i]) > 1e-10 ? 1/abs(x[i]) : 1e6)
        if abs(tap[i]-1.0) > 1e-6 || abs(shift[i]) > 1e-6
            w *= 50.0
        end
        if ratea[i] > 0
            w *= (1 + ratea[i])
        end
        push!(raw, w)
    end
    isempty(raw) && (raw = [1.0])
    maxw = maximum(raw)

    for i in eachindex(f)
        on[i]==1 || continue
        u, v = f[i], t[i]
        present[u] && present[v] || continue
        wfloat = (abs(x[i]) > 1e-10 ? 1/abs(x[i]) : 1e6)
        if abs(tap[i]-1.0) > 1e-6 || abs(shift[i]) > 1e-6
            wfloat *= 50.0
        end
        if ratea[i] > 0
            wfloat *= (1 + ratea[i])
        end
        w = Int32(round(1 + 10_000 * (wfloat / maxw)))
        key = u < v ? (u,v) : (v,u)
        adj_dict[key] = get(adj_dict, key, 0) + w
    end

    vw = zeros(Float64, maxbus)
    if haskey(case["bus"], "pd")
        for i in eachindex(buses)
            vw[buses[i]] += max(case["bus"]["pd"][i], 0.0)
        end
    end
    if haskey(case, "gen") && !isempty(case["gen"])
        gtbl = case["gen"]
        for i in eachindex(gtbl["gen_bus"])
            vw[gtbl["gen_bus"][i]] += max(gtbl["pmax"][i], 0.0)
        end
    end
    for b in 1:maxbus
        present[b] || (vw[b] = 0)
        vw[b] = vw[b] <= 0 ? 1.0 : vw[b]
    end
    vmax = maximum(vw[present])
    vwgt = Int32.(round.(1 .+ 1000 .* vw ./ vmax))

    return g, vwgt, adj_dict, present, posmap
end

# ------------------------------------------------------------------
# 3) Métrica de qualidade: edge-cut ratio e tamanhos
# ------------------------------------------------------------------
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
        zu = zone[pos[u]]; zv = zone[pos[v]]
        w = (abs(x[i]) > 1e-10 ? 1/abs(x[i]) : 1e6)
        if abs(tap[i]-1.0) > 1e-6 || abs(shift[i]) > 1e-6; w *= 50.0; end
        if ratea[i] > 0; w *= (1 + ratea[i]); end
        w_tot += w
        if zu != zv; w_cut += w; end
    end
    return w_cut / max(w_tot, eps())
end

# ------------------------------------------------------------------
# 4) Política de escolha de k (auto)
# ------------------------------------------------------------------
function choose_k(case; max_zone::Int=5000, kmin::Int=2, kmax::Int=2000)
    n = length(case["bus"]["bus_i"])
    k = max(kmin, ceil(Int, n / max_zone))
    return min(max(k, kmin), kmax)
end

# --------------------------------------------
# Pós-processos: conectividade + geração/slack
# --------------------------------------------
function fix_connectivity_and_slack!(case::Dict{String,Any}, zone::Vector{Int}, k::Int)
    buses = case["bus"]["bus_i"]
    f = case["branch"]["f_bus"]; t = case["branch"]["t_bus"]; on = case["branch"]["br_status"]

    # vizinhança por bus (apenas ramos ativos)
    neigh = Dict{Int,Vector{Int}}()
    for i in eachindex(f)
        on[i]==1 || continue
        push!(get!(neigh, f[i], Int[]), t[i])
        push!(get!(neigh, t[i], Int[]), f[i])
    end
    pos = Dict(buses[i]=>i for i in eachindex(buses))

    # 1) garantir conectividade por zona: se houver componentes desconexas
    #    numa mesma zona, "puxa" cada componente menor para a zona mais
    #    incidente na borda dessa componente.
    for z in 1:k
        nodes = [buses[i] for i in eachindex(buses) if zone[i]==z]
        isempty(nodes) && continue
        seen = Set{Int}()
        comps = Vector{Vector{Int}}()
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
        sort!(comps, by=length, rev=true)  # primeira é a maior (fica)
        for comp in comps[2:end]
            # zona mais "conectada" à fronteira dessa componente
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

    # 2) garantir pelo menos uma barra com geração (>0 pmax) por zona.
    gen = get(case, "gen", Dict{String,Any}())
    isempty(gen) && return zone

    genbus = gen["gen_bus"]; status = gen["gen_status"]; pmax = gen["pmax"]
    pos = Dict(case["bus"]["bus_i"][i]=>i for i in eachindex(case["bus"]["bus_i"]))

    # pmax por zona
    pmax_by_zone = Dict(z=>0.0 for z in 1:k)
    for i in eachindex(genbus)
        status[i]==1 || continue
        b  = genbus[i]
        z  = zone[pos[b]]
        pM = max(pmax[i], 0.0)
        pmax_by_zone[z] = get(pmax_by_zone, z, 0.0) + pM
    end

    # grafo entre zonas (arestas de corte)
    neighZ = Dict{Int,Set{Int}}()
    for i in eachindex(f)
        on[i]==1 || continue
        zf = zone[pos[f[i]]]; zt = zone[pos[t[i]]]
        if zf != zt
            push!(get!(neighZ, zf, Set{Int}()), zt)
            push!(get!(neighZ, zt, Set{Int}()), zf)
        end
    end

    # lista de geradores "movíveis" por zona (pmax>0)
    gens_in_zone = Dict(z=>Int[] for z in 1:k)
    for i in eachindex(genbus)
        status[i]==1 || continue
        if pmax[i] > 1e-6
            bz = zone[pos[genbus[i]]]
            push!(gens_in_zone[bz], i)
        end
    end

    # para cada zona sem pmax, "puxa" um gerador de uma zona vizinha (ou qualquer outra)
    for z in 1:k
        if get(pmax_by_zone, z, 0.0) > 1e-9
            continue
        end
        moved = false
        for zv in sort(collect(get(neighZ, z, Set{Int}())))  # preferir vizinhas
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
            # pega de qualquer outra zona
            for zz in 1:k
                zz == z && continue
                if !isempty(get(gens_in_zone, zz, Int[]))
                    gi = pop!(gens_in_zone[zz])
                    b  = genbus[gi]
                    zone[pos[b]] = z
                    pmax_by_zone[z]  += max(pmax[gi],0.0)
                    pmax_by_zone[zz]  = max(pmax_by_zone[zz] - max(pmax[gi],0.0), 0.0)
                    push!(gens_in_zone[z], gi)
                    break
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


# ------------------------------------------------------------------
# 5) Pipeline em lote usando METIS
# ------------------------------------------------------------------
using SimpleWeightedGraphs  # <- garanta que está no topo do arquivo

function metis_partition(case::Dict{String,Any}; k::Int, balance_tol::Float64=0.10)
    @assert k >= 1 "k deve ser ≥ 1"
    @assert haskey(case, "bus") && haskey(case, "branch") "case incompleto"

    # --- coleta dados básicos
    buses  = case["bus"]["bus_i"]
    present = Set(buses)
    posmap = Dict(buses[i] => i for i in eachindex(buses))

    f   = case["branch"]["f_bus"]
    t   = case["branch"]["t_bus"]
    x   = case["branch"]["br_x"]
    r   = case["branch"]["br_r"]
    on  = case["branch"]["br_status"]
    tap   = get(case["branch"], "tap",   fill(1.0, length(f)))
    shift = get(case["branch"], "shift", fill(0.0, length(f)))
    ratea = get(case["branch"], "rate_a", fill(0.0, length(f)))

    # --- mapeia bus_id (globais) -> 1..n (locais)
    glob2loc = Dict{Int,Int}()
    loc2glob = Vector{Int}(undef, length(buses))
    for (i,b) in enumerate(buses)
        glob2loc[b] = i
        loc2glob[i] = b
    end
    n = length(loc2glob)

    # --- 1) calcula pesos brutos para normalizar
    raw_w = Float64[]
    for i in eachindex(f)
        on[i]==1 || continue
        u, v = f[i], t[i]
        (u in present && v in present) || continue
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

    # --- 2) acumula pesos por aresta (sem duplicar/atualizar no grafo)
    wdict = Dict{Tuple{Int,Int}, Int}()  # pesos inteiros para o Metis
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

    # --- 3) constrói grafo ponderado uma única vez
    G = SimpleWeightedGraph(n)
    for ((u,v), w) in wdict
        add_edge!(G, u, v, w)
    end

    # --- 4) chama METIS (alto nível); retorna partições 1..k
    parts = Metis.partition(G, k; alg = :KWAY)

    # --- 5) monta vetor zone alinhado ao case + pós-processo
    Z = fill(1, length(buses))
    @inbounds for (i, bid) in enumerate(buses)
        Z[i] = parts[glob2loc[bid]]
    end

    fix_connectivity_and_slack!(case, Z, k)
    return Z
end



# --------------------------------------------
# Execução direta (edite os caminhos/flags se quiser)
# --------------------------------------------
partition_cases_metis(
    raw"C:\Users\franc\PowerOptimization.jl\test\zones";
    ks = :auto,                 # ou ks=[25,50,100]
    max_zone_size = 350,       # “3k–5k”
    target_edgecut = 0.02,      # 2%
    balance_tol = 0.10,         # ±10%
    glob = "*.m",
    overwrite = false,
)
