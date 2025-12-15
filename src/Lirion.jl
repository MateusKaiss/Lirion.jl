module Lirion

# ------------------------------------------------------------------
# Configuração de threads (Julia + BLAS)
# ------------------------------------------------------------------
for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
    if !haskey(ENV, k)
        ENV[k] = "1"
    end
end

using LinearAlgebra
try
    LinearAlgebra.BLAS.set_num_threads(1)
catch
end

# ------------------------------------------------------------------
# Pacotes
# ------------------------------------------------------------------
using PlotlyLight
using JuMP
using Ipopt
using Gurobi
using DataFrames
using CSV
using Dates
using Statistics
using InfrastructureModels
using Logging
using ArgParse
using SCS
using OSQP
using Base.Threads
using MathOptInterface
using COSMO
using Clarabel

const MOI = MathOptInterface
const _FAST_USE_LOGS = true

function __init__()
    if _FAST_USE_LOGS
        try
            @info "[Lirion] Threads: Julia=$(nthreads()) | BLAS=$(LinearAlgebra.BLAS.get_num_threads())"
        catch
        end
    end
end

# ------------------------------------------------------------------
# Tipos básicos
# ------------------------------------------------------------------
include("types.jl")

# ------------------------------------------------------------------
# Collections: dados estáticos + timeseries + SystemData
# ------------------------------------------------------------------

# Timeseries utils (opera em Dict) + structs de TS
include("collections/timeseries/timeseries_utils.jl")
include("collections/timeseries/bus_timeseries.jl")
include("collections/timeseries/prosumers_timeseries.jl")
include("collections/timeseries/consumer_timeseries.jl")
include("collections/timeseries/gencost_timeseries.jl")

# Coleções estáticas
include("collections/bus.jl")
include("collections/branch.jl")
include("collections/generator.jl")
include("collections/battery.jl")
include("collections/prosumers.jl")
include("collections/objective.jl")
include("collections/power_balance.jl")

# SystemData + ADMM/estruturas de sistema distribuído
include("collections/system_data.jl")

# Pré-processamento (preprocess_case_data!, convert_to_pu!, etc.)
include("collections/preprocess.jl")

# Args e inputs (build_inputs usa preprocess_case_data!, build_system_data, etc.)
include("args/args.jl")
include("collections/inputs.jl")
include("collections/admm_structs.jl")

# ------------------------------------------------------------------
# Parsing de casos (.m / .csv)
# ------------------------------------------------------------------
include("parsing/parser.jl")
include("parsing/csv_parser.jl")

# ------------------------------------------------------------------
# Construção do modelo JuMP
# ------------------------------------------------------------------
include("build_model/bus.jl")
include("build_model/branch.jl")
include("build_model/generator.jl")
include("build_model/battery.jl")
include("build_model/prosumer.jl")
include("build_model/active_power_balance.jl")
include("build_model/reactive_power_balance.jl")
include("build_model/objective_function.jl")
include("build_model/build_model.jl")

# ------------------------------------------------------------------
# Pós-processamento / resultados
# ------------------------------------------------------------------
include("build_results/build_results.jl")

# ------------------------------------------------------------------
# Métodos distribuídos (ADMM, etc.)
# ------------------------------------------------------------------
include("distributed/distributed_opf.jl")

include("distributed/admm/utils.jl")
include("distributed/admm/preparation.jl")
include("distributed/admm/admm_core.jl")

include("distributed/admm/terms_bus.jl")
include("distributed/admm/terms_branch.jl")

include("distributed/admm/admm_dc.jl")
include("distributed/admm/admm_ac.jl")
include("distributed/admm/admm_socp.jl")

# ------------------------------------------------------------------
# Entrada principal (run_lirion / main / julia_main)
# ------------------------------------------------------------------
include("main.jl")

end
