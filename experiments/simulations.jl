# =========================
# batch_experiments.jl
# =========================
using Dates, Printf, CSV, DataFrames

# --- tente carregar Lirion no ambiente atual ---
try
    @eval using Lirion
catch
    @warn "Não foi possível carregar Lirion; certifique-se de estar com --project correto."
end

# --------- Descoberta de casos na pasta atual ---------
"""
    discover_cases()

Retorna vetor de (label, path) para os 3 casos: :std, :api, :sad (se existirem).
"""
function discover_cases()
    files = readdir(pwd(); join=true)
    pick(pat) = begin
        xs = filter(f -> occursin(pat, basename(f)) && endswith(f, ".m"), files)
        isempty(xs) ? nothing : xs[1]
    end
    std = pick(r"_lmbd\.m$")           # std tem 1 underscore antes de lmbd
    api = pick(r"_lmbd__api\.m$")
    sad = pick(r"_lmbd__sad\.m$")
    out = Tuple{Symbol,String}[]
    std !== nothing && push!(out, (:std, std))
    api !== nothing && push!(out, (:api, api))
    sad !== nothing && push!(out, (:sad, sad))
    return out
end

# --------- Mapeamento de variantes por família ---------
# --------- Mapeamento de variantes por família ---------
const VARIANTS = Dict(
    :dc   => [(:admm_bus,    "distributed-admm-dc"),          # <- sem "bus" na string
              (:admm_branch, "distributed-admm-branch-dc"),
              (:admm_hybrid, "distributed-admm-hybrid-dc")],
    :ac   => [(:admm_bus,    "distributed-admm-ac"),          # <- sem "bus" na string
              (:admm_branch, "distributed-admm-branch-ac"),
              (:admm_hybrid, "distributed-admm-hybrid-ac")],
    :socp => [(:admm_bus,    "distributed-admm-socp"),        # <- sem "bus" na string
              (:admm_branch, "distributed-admm-branch-socp"),
              (:admm_hybrid, "distributed-admm-hybrid-socp")]
)

# --------- Helpers simples ---------
getflag(args::Vector{String}, key::String) = begin
    i = findfirst(==(key), args)
    (i === nothing || i == length(args)) ? missing : args[i+1]
end

# --------- Builder de ARGS p/ uma execução ---------
"""
    build_args(case_path; model_key, family, extra_args=String[])

Monta: [case_path, "--model", <modelo>, ...extra_args]
Use `extra_args` para passar tudo (ex.: solver, ρ, tolerâncias, warm, etc).
"""
function build_args(case_path::String; model_key::Symbol, family::Symbol, extra_args::Vector{String}=String[])
    model_str = only(filter(p -> p[1] == model_key, VARIANTS[family]))[2]
    args = String[case_path, "--model", model_str]
    append!(args, extra_args)
    return args
end

# --------- Execução de um job via sua entrypoint de ARGS ---------
"""
    run_once(args::Vector{String}) -> (status, obj, iters, time_s)

Chama sua entrada ARGS. Esperado: (total_real_cost, partitions, iter_done, elapsed_s).
"""
function run_once(args::Vector{String})
    try
        if isdefined(Lirion, :run_with_args)
            res = Lirion.run_with_args(args)
        elseif isdefined(Lirion, :main)
            res = Lirion.main(args)
        elseif isdefined(Main, :main)
            res = Main.main(args)
        else
            error("Não encontrei entrypoint por ARGS (ex.: Lirion.run_with_args / Lirion.main).")
        end
        obj, _, iters, time_s = res
        return ("ok", float(obj), Int(iters), float(time_s))
    catch err
        @warn "Falhou: $err" args=args
        return ("error", NaN, 0, NaN)
    end
end

# --------- Runner de lote ---------
"""
    run_batch(dir=pwd(); family=:socp, which=:all, extra_args=String[], outfile=nothing)

- `family`: :dc | :ac | :socp
- `which`:  :all | vetor de Symbols em [:admm_bus, :admm_branch, :admm_hybrid]
- `extra_args`: será anexado a *todas* as execuções (ex.: ["--optimizer","Clarabel.Optimizer", ...])
- `outfile`: caminho do CSV; se `nothing`, gera automático.
"""
function run_batch(dir::AbstractString=pwd(); family::Symbol=:socp, which=:all,
                   extra_args::Vector{String}=String[], outfile=nothing)
    haskey(VARIANTS, family) || error("family deve ser :dc, :ac ou :socp")
    cd(dir) do
        cases = discover_cases()
        isempty(cases) && error("Nenhum caso encontrado nesta pasta.")
        variants = which === :all ? VARIANTS[family] : filter(v -> v[1] in which, VARIANTS[family])

        rows = NamedTuple[]
        for (label, cpath) in cases
            for (vkey, model_str) in variants
                args = build_args(cpath; model_key=vkey, family=family, extra_args=extra_args)
                status, obj, iters, time_s = run_once(args)

                push!(rows, (
                    timestamp = Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"),
                    case_label = String(label),
                    case_name  = basename(cpath),
                    base_model = String(family),
                    variant    = String(vkey),
                    model_str  = model_str,
                    optimizer  = getflag(args, "--optimizer"),
                    status     = status,
                    objective  = obj,
                    iterations = iters,
                    time_s     = time_s,
                    args_line  = join(args, " ")
                ))
                @info "Concluído" case=basename(cpath) model=model_str variant=vkey status iters time_s obj=@sprintf("%.8f", obj)
            end
        end

        df = DataFrame(rows)
        if isnothing(outfile)
            ts = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
            outfile = joinpath(dir, "results_$(String(family))_$ts.csv")
        end
        CSV.write(outfile, df)
        @info "CSV salvo" outfile rows=length(rows)
        return outfile, df
    end
end

# =========================
# ARGSETS DE EXEMPLO (edite livremente)
# =========================
# Você pode usar estes vetores diretamente no `extra_args` do run_batch.

extra_dc = String[
    "--optimizer","OSQP.Optimizer",
    "--admm-max-iter","700",
    "--admm-rho-v","1e3",
    "--admm-rho-flow-p","1e3",
    "--admm-rho-flow-q","1e3",
    "--admm-eps-primal","1e-4",
    "--admm-eps-dual","1e-4",
    "--admm-adaptive","on",
    "--admm-mu-incr","3.0",
    "--admm-mu-decr","3.0",
    "--admm-tau-incr","1.1",
    "--admm-tau-decr","0.9",
    # "--admm-warm-start","on",
    # "--admm-acceleration","on"
]

extra_socp = String[
    "--optimizer","Clarabel.Optimizer",   # ou COSMO.Optimizer / SCS.Optimizer
    "--admm-max-iter","700",
    "--admm-rho-v","1e3",
    "--admm-rho-flow-p","1e3",
    "--admm-rho-flow-q","1e3",
    "--admm-eps-primal","1e-4",
    "--admm-eps-dual","1e-4",
    "--admm-adaptive","on",
    "--admm-mu-incr","3.0",
    "--admm-mu-decr","3.0",
    "--admm-tau-incr","1.1",
    "--admm-tau-decr","0.9",
    "--admm-warm-start","on",
    "--admm-acceleration","on"
]

extra_ac = String[
    "--optimizer","Ipopt.Optimizer",
    "--admm-max-iter","700",
    "--admm-rho-v","1e3",
    "--admm-rho-flow-p","1e3",
    "--admm-rho-flow-q","1e3",
    "--admm-eps-primal","1e-4",
    "--admm-eps-dual","1e-4",
    "--admm-adaptive","on",
    "--admm-mu-incr","3.0",
    "--admm-mu-decr","3.0",
    "--admm-tau-incr","1.1",
    "--admm-tau-decr","0.9",
    "--admm-warm-start","on",
    "--admm-acceleration","on"
]

run_batch("/home/mateus-kaiss/Projects/experiments/03-bus";
          family=:dc, extra_args=extra_dc)