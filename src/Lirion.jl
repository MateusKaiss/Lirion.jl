module Lirion
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
using HSL_jll

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

include("types.jl")
include("args/args.jl")

include("parsing/parser.jl")
include("parsing/csv_parser.jl")

include("preprocessing/timeseries_utils.jl")
include("preprocessing/pre_processing.jl")
include("preprocessing/parameter_extractors.jl")

include("build_model/helpers.jl")
include("build_model/bus.jl")
include("build_model/branch.jl")
include("build_model/generator.jl")
include("build_model/battery.jl")
include("build_model/prosumer.jl")
include("build_model/active_power_balance.jl")
include("build_model/reactive_power_balance.jl")
include("build_model/objective_function.jl")
include("build_model/build_model.jl")

include("build_results/build_results.jl")

include("distributed/distributed_opf.jl")
include("distributed/admm/partitions.jl")
include("distributed/admm/admm_core.jl")
include("distributed/admm/admm_dc.jl")
include("distributed/admm/admm_ac.jl")
include("distributed/admm/admm_socp.jl")

include("main.jl")

end
