module Lirion

# ------------------------------------------------------------------
# Thread configuration (Julia + BLAS)
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
# Packages
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
# Basic types
# ------------------------------------------------------------------
include("types.jl")

# ------------------------------------------------------------------
# Collections: static data + time series + SystemData
# ------------------------------------------------------------------
include("collections/timeseries/timeseries_utils.jl")
include("collections/timeseries/bus_timeseries.jl")
include("collections/timeseries/prosumers_timeseries.jl")
include("collections/timeseries/consumer_timeseries.jl")
include("collections/timeseries/gencost_timeseries.jl")

include("collections/bus.jl")
include("collections/branch.jl")
include("collections/generator.jl")
include("collections/battery.jl")
include("collections/prosumers.jl")
include("collections/objective.jl")
include("collections/power_balance.jl")
include("collections/system_data.jl")
include("collections/preprocess.jl")

# ------------------------------------------------------------------
# Inputs structs (legacy CLI + shared input pipeline)
# ------------------------------------------------------------------
include("collections/inputs.jl")
include("collections/admm_structs.jl")

# ------------------------------------------------------------------
# Case parsing (.m / .csv)
# ------------------------------------------------------------------
include("parsing/parser.jl")
include("parsing/csv_parser.jl")

# ------------------------------------------------------------------
# Model construction and centralized solution
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
# Post-processing / results
# ------------------------------------------------------------------
include("build_results/build_results.jl")
include("build_results/build_plots.jl")

# ------------------------------------------------------------------
# Distributed methods
# ------------------------------------------------------------------
include("distributed/distributed_opf.jl")

include("distributed/admm/admm_support.jl")
include("distributed/admm/buffers.jl")
include("distributed/admm/objective_terms.jl")
include("distributed/admm/consensus.jl")
include("distributed/admm/duals.jl")
include("distributed/admm/residuals.jl")
include("distributed/admm/rho.jl")

include("distributed/admm/admm_dc.jl")
include("distributed/admm/admm_ac.jl")
include("distributed/admm/admm_socp.jl")

# ------------------------------------------------------------------
# Legacy runner / CLI entrypoints
# (run_lirion is defined here today, so it must be included
#  BEFORE solve.jl which calls run_lirion)
# ------------------------------------------------------------------
include("api.jl")

# ------------------------------------------------------------------
# New function-based API (no Args required)
# ------------------------------------------------------------------
include("solve.jl")

# ------------------------------------------------------------------
# Public exports (API)
# ------------------------------------------------------------------
export solve, dc, ac, socp,
    centralized, admm_bus, admm_branch, admm_hybrid,
    partition_from_case, partition_metis

end
