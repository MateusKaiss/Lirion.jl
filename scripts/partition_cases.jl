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
# Construção do grafo e particionamento (espectral)
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

# ---------- helpers espectrais seguros ----------
function fiedler_vector(L::AbstractMatrix{<:Real})
    Lf = SparseMatrixCSC{Float64,Int}(L)
    n = size(Lf,1)
    if n <= 2
        return nothing
    end
    if n <= 200
        vals, vecs = eigen(Symmetric(Matrix(Lf)))
        return length(vals) >= 2 ? vecs[:,2] : nothing
    else
        try
            vals, vecs = Arpack.eigs(Lf; nev=2, which=:SR)
            return size(vecs,2) >= 2 ? vecs[:,2] : nothing
        catch
            return nothing
        end
    end
end

function partition_spectral(case, g::SimpleGraph, k::Int, present::BitVector, posmap::Dict{Int,Int})
    buses = collect(keys(posmap)) |> sort
    zone = fill(1, length(buses))
    bucket = Dict(1 => collect(eachindex(buses)))
    nextlabel = 2

    while length(bucket) < k
        sizes = Dict(z => length(v) for (z, v) in bucket)
        zkey = argmax(sizes)    
        idxs = bucket[zkey]
        sub_nodes = buses[idxs]
        if length(sub_nodes) <= 2
            break
        end

        sg, mapping = induced_subgraph(g, sub_nodes)
        if nv(sg) <= 2
            break
        end

        L = SparseMatrixCSC{Float64,Int}(laplacian_matrix(sg))
        f = fiedler_vector(L)
        f === nothing && break

        thresh = median(f)
        left  = findall(<(thresh), f)
        right = findall(>=(thresh), f)
        if isempty(left) || isempty(right)
            break
        end

        for (side, list) in enumerate((left, right))
            for vh in list
                b = mapping[vh]   
                pos = posmap[b]      
                zone[pos] = (side == 1 ? zkey : nextlabel)
            end
        end

        bucket[zkey]      = [i for i in eachindex(zone) if zone[i]==zkey]
        bucket[nextlabel] = [i for i in eachindex(zone) if zone[i]==nextlabel]
        nextlabel += 1
        if nextlabel > k
            break
        end
    end
    fix_connectivity_and_slack!(case, zone, k)
    return zone
end

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
function partition_cases(indir::AbstractString=DEFAULT_DIR; ks::Vector{Int}=[2,3,4,5],
    seed::Integer=42, glob::AbstractString="*.m", overwrite::Bool=false)

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
            z = partition_spectral(case, g, k, present, posmap)
            apply_zones!(case, z)
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
    )
end

# Para rodar direto via include:
partition_cases(
    raw"C:\Users\franc\PowerOptimization.jl\test\zones";
    ks=[2,3,4,5],
    seed=42,
    glob="*.m",
    overwrite=false
)
