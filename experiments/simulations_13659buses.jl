# experiments/simulations.jl
# Batch para Lirion: variações ADMM (bus/branch/hybrid) + centralizado — SEM MÓDULO

using Glob
using DataFrames
using CSV
using Dates
using Printf
using JuMP
import MathOptInterface as MOI
# ----------------------- Utils -----------------------
const _FMT_TIME = dateformat"yyyy-mm-dd HH:MM:SS"

_family_to_modelkey(fam::Symbol) = begin
    fam === :dc   && return "dc"
    fam === :ac   && return "ac"
    fam === :socp && return "socp"
    error("family inválido: $fam (use :dc, :ac ou :socp)")
end

_variants() = (:bus, :branch, :hybrid)

_modelkey_distributed(family::Symbol, variant::Symbol) = begin
    base = _family_to_modelkey(family)
    variant === :bus    && return "distributed-admm-$base"
    variant === :branch && return "distributed-admm-branch-$base"
    variant === :hybrid && return "distributed-admm-hybrid-$base"
    error("variant inválido: $variant")
end

function _list_cases(dir::AbstractString)
    isdir(dir) || error("Diretório não existe: $dir")
    return sort(glob("*.m", dir))
end

_now_str() = Dates.format(Dates.now(), _FMT_TIME)

# Extrai objetivo do retorno centralizado
function _objective_from_centralized(result_tuple)::Float64
    # Esperado: (result_dict, result_model, data)
    length(result_tuple) >= 2 || error("Retorno centralizado inesperado (tuple curto).")
    result_model = result_tuple[2]
    try
        return JuMP.objective_value(result_model)
    catch
        try
            return value(result_model[:generation_cost])
        catch
            error("Não foi possível extrair objetivo do modelo centralizado.")
        end
    end
end

# Extrai (obj, iters, tempo_s) do retorno distribuído
function _triple_from_distributed(result_any)
    # runners retornam: (total_real_cost, partitions, iter_done, elapsed_s)
    if result_any isa Tuple && length(result_any) >= 4
        obj = Float64(result_any[1])
        iters = Int(result_any[3])
        time_s = Float64(result_any[4])
        return obj, iters, time_s
    end
    error("Retorno distribuído inesperado; esperado tuple (obj, ..., iters, time_s).")
end

# Executa Lirion.main e captura erro/resultado
function _call_lirion_main(argsvec::Vector{String})
    try
        res = Lirion.main(argsvec)
        return (status = "ok", result = res, error = "")
    catch err
        return (status = "error", result = nothing, error = sprint(showerror, err))
    end
end

# ----------------------- Batch: Distribuído -----------------------
"""
run_batch(dir; family=:dc, extra_args=String[], outpath=nothing)

- Roda as *três variantes distribuídas* (bus/branch/hybrid) do `family` (dc/ac/socp)
  para todos os `.m` dentro de `dir`.
- CSV: [:timestamp, :case, :family, :variant, :objective, :iters, :time_s, :status, :error]
"""
function run_batch(dir::AbstractString; family::Symbol=:dc,
                   extra_args::Vector{String}=String[],
                   outpath::Union{Nothing,String}=nothing)

    famkey = _family_to_modelkey(family)
    cases = _list_cases(dir)
    isempty(cases) && @warn "Nenhum .m encontrado em $dir"

    rows = Vector{NamedTuple}()
    for f in cases
        for variant in _variants()
            modelkey = _modelkey_distributed(family, variant)
            args = String[f, "--model", modelkey]
            append!(args, extra_args)

            @info "[ADMM] $(basename(f)) | $modelkey"
            call = _call_lirion_main(args)

            if call.status == "ok"
                try
                    obj, iters, time_s = _triple_from_distributed(call.result)
                    push!(rows, (
                        timestamp = _now_str(),
                        case = f, family = famkey, variant = String(variant),
                        objective = obj, iters = iters, time_s = time_s,
                        status = "ok", error = "",
                    ))
                catch err
                    push!(rows, (
                        timestamp = _now_str(),
                        case = f, family = famkey, variant = String(variant),
                        objective = NaN, iters = -1, time_s = NaN,
                        status = "error", error = sprint(showerror, err),
                    ))
                end
            else
                push!(rows, (
                    timestamp = _now_str(),
                    case = f, family = famkey, variant = String(variant),
                    objective = NaN, iters = -1, time_s = NaN,
                    status = "error", error = call.error,
                ))
            end
        end
    end

    df = DataFrame(rows)
    outcsv = isnothing(outpath) ? joinpath(dir, "admm_results_$(famkey).csv") : outpath
    CSV.write(outcsv, df)
    @info "CSV (ADMM) salvo em: $outcsv"
    show(df)
    println()  # pular linha
    return df
end

# ----------------------- Batch: Centralizado -----------------------
"""
run_batch_centralized(dir; family=:dc, extra_args=String[], outpath=nothing)

- Roda o **centralizado** (dc/ac/socp) para todos os `.m` em `dir`.
- CSV: [:timestamp, :case, :family, :objective, :status, :error]
"""
function run_batch_centralized(dir::AbstractString; family::Symbol=:dc,
                               extra_args::Vector{String}=String[],
                               outpath::Union{Nothing,String}=nothing)

    famkey = _family_to_modelkey(family)
    cases = _list_cases(dir)
    isempty(cases) && @warn "Nenhum .m encontrado em $dir"

    rows = Vector{NamedTuple}()
    for f in cases
        args = String[f, "--model", famkey]
        append!(args, extra_args)

        @info "[CENT] $(basename(f)) | $famkey"
        call = _call_lirion_main(args)

        if call.status == "ok"
            try
                obj = _objective_from_centralized(call.result)
                push!(rows, (
                    timestamp = _now_str(),
                    case = f, family = famkey, objective = obj,
                    status = "ok", error = "",
                ))
            catch err
                push!(rows, (
                    timestamp = _now_str(),
                    case = f, family = famkey, objective = NaN,
                    status = "error", error = sprint(showerror, err),
                ))
            end
        else
            push!(rows, (
                timestamp = _now_str(),
                case = f, family = famkey, objective = NaN,
                status = "error", error = call.error,
            ))
        end
    end

    df = DataFrame(rows)
    outcsv = isnothing(outpath) ? joinpath(dir, "centralized_results_$(famkey).csv") : outpath
    CSV.write(outcsv, df)
    @info "CSV (centralizado) salvo em: $outcsv"
    show(df)
    println()  # pular linha
    return df
end

# ----------------------- EXEMPLOS (opcional) -----------------------
# Deixe comentado e rode a partir do REPL.
extra_dc = String[
    "--optimizer","Gurobi.Optimizer",
    "--admm-max-iter","2000",
    "--admm-rho","5e4",
    "--admm-rho-v","5e4",
    "--admm-rho-flow-p","1e4",
    "--admm-rho-flow-q","1e4",
    "--admm-eps-primal","1e-4",
    "--admm-eps-dual","1e-4",
    "--admm-adaptive","on",
    "--admm-mu-incr","3.0",
    "--admm-mu-decr","3.0",
    "--admm-tau-incr","1.05",
    "--admm-tau-decr","0.95",
]


# run_batch("/home/mateus-kaiss/Projects/experiments/13659-bus"; family=:dc, extra_args=extra_dc)
# run_batch_centralized("/home/mateus-kaiss/Projects/experiments/13659-bus";
#     family=:dc, extra_args=["--optimizer","Gurobi.Optimizer"])


extra_ac = String[
    "--admm-max-iter","2000",
    "--admm-rho","5e7",          # base alto
    "--admm-rho-v","5e7",        # |V|
    "--admm-rho-theta","5e7",        # |V|
    "--admm-rho-flow-p","1e5",   # P_flow
    "--admm-rho-flow-q","1e5",   # Q_flow
    "--admm-eps-primal","1e-4",
    "--admm-eps-dual","1e-4",
    "--admm-adaptive","on",
    "--admm-mu-incr","2.0",
    "--admm-mu-decr","2.0",
    "--admm-tau-incr","1.2",
    "--admm-tau-decr","0.8",
    # "--admm-warm-start","on", 
    # "--admm-accel",             # opcional
]

# # ----- Rodar AC (distribuído: bus/branch/hybrid) -----
run_batch("/home/mateus-kaiss/Projects/experiments/13659-bus";
    family=:ac, extra_args=extra_ac)

# # # # ----- Rodar AC (centralizado) -----
run_batch_centralized("/home/mateus-kaiss/Projects/experiments/13659-bus";
    family=:ac, extra_args=["--optimizer","Ipopt.Optimizer"])


extra_socp = String[
    "--optimizer",
    """optimizer_with_attributes(Clarabel.Optimizer,
        # mesmas tolerâncias/iters do centralizado (preset antichoque)
        "tol_feas" => 3e-4,
        "tol_gap_abs" => 3e-4,
        "tol_gap_rel" => 3e-4,
        "max_iter" => 3000,

        # fatoração direta estável para KKT de SOCP
        "direct_solve_method" => :qdldl,

        # regularização mais forte
        "static_regularization_enable" => true,
        "dynamic_regularization_enable" => true,
        "dynamic_regularization_eps"   => 1e-8,
        "dynamic_regularization_delta" => 1e-3,

        # busca de linha mais conservadora
        "max_step_fraction" => 0.5,
        "linesearch_backtrack_step" => 0.5,
        "min_terminate_step_length" => 1e-7,
        "min_switch_step_length" => 5e-2,

        # equilibrate menos agressivo
        "equilibrate_enable" => true,
        "equilibrate_min_scaling" => 1e-3,
        "equilibrate_max_scaling" => 1e3,

        # refinamento iterativo um pouco mais forte
        "iterative_refinement_reltol" => 5e-10,
        "iterative_refinement_abstol" => 5e-10
    )""",

    "--admm-max-iter","2000",
    "--admm-rho","1e5",
    "--admm-rho-v","1e5",
    "--admm-rho-flow-p","1e5",
    "--admm-rho-flow-q","1e5",
    "--admm-eps-primal","1e-4",
    "--admm-eps-dual","1e-4",
    "--admm-adaptive","on",
    "--admm-mu-incr","3.0",
    "--admm-mu-decr","3.0",
    "--admm-tau-incr","1.1",
    "--admm-tau-decr","0.9",
]


# # # ----- Rodar SOCP (distribuído: bus/branch/hybrid) -----
# run_batch("/home/mateus-kaiss/Projects/experiments/13659-bus";
#     family=:socp, extra_args=extra_socp)

# # # # ----- SOCP (centralizado) -----
# run_batch_centralized("/home/mateus-kaiss/Projects/experiments/13659-bus";
#     family=:socp, extra_args=[
#         "--optimizer",
#         """optimizer_with_attributes(Clarabel.Optimizer,
#             # tolerâncias moderadas e mais iterações
#             "tol_feas" => 3e-4,
#             "tol_gap_abs" => 3e-4,
#             "tol_gap_rel" => 3e-4,
#             "max_iter" => 3000,

#             # deixe QDLDL (evita ZeroPivot do CHOLMOD)
#             "direct_solve_method" => :qdldl,

#             # regularização mais forte
#             "static_regularization_enable" => true,
#             "dynamic_regularization_enable" => true,
#             "dynamic_regularization_eps"   => 1e-8,
#             "dynamic_regularization_delta" => 1e-3,

#             # busca de linha mais conservadora
#             "max_step_fraction" => 0.5,
#             "linesearch_backtrack_step" => 0.5,
#             "min_terminate_step_length" => 1e-7,
#             "min_switch_step_length" => 5e-2,

#             # equilíbrio menos agressivo
#             "equilibrate_enable" => true,
#             "equilibrate_min_scaling" => 1e-3,
#             "equilibrate_max_scaling" => 1e3,

#             # refinamento iterativo mais forte
#             "iterative_refinement_reltol" => 5e-10,
#             "iterative_refinement_abstol" => 5e-10
#         )"""
#     ])