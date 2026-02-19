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
using Logging
using SCS
using OSQP
using Base.Threads
using MathOptInterface
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
include("parsing/raw_structs.jl")

# ------------------------------------------------------------------
# Collections: static data + time series + SystemData
# ------------------------------------------------------------------
include("collections/timeseries_utils.jl")
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
include("collections/preprocess_case.jl")
include("collections/distributed_partitioning.jl")

# ------------------------------------------------------------------
# Inputs structs (legacy CLI + shared input pipeline)
# ------------------------------------------------------------------
include("collections/build_inputs.jl")
include("collections/admm_structs.jl")

# ------------------------------------------------------------------
# Case parsing (.m / .csv)
# ------------------------------------------------------------------
include("parsing/parser_utils.jl")
include("parsing/matpower_parser.jl")
include("parsing/csv_parser.jl")

# ------------------------------------------------------------------
# Model construction and centralized solution
# ------------------------------------------------------------------
# --- BUS
include("build_model/bus/variables.jl")
include("build_model/bus/dispatch.jl")

# --- BRANCH
include("build_model/branch/dc.jl")
include("build_model/branch/ac.jl")
include("build_model/branch/socp_rt.jl")
include("build_model/branch/socp.jl")
include("build_model/branch/dispatch.jl")

# --- GENERATOR
include("build_model/generator/variables.jl")
include("build_model/generator/constraints.jl")
include("build_model/generator/dispatch.jl")

# --- BATTERY
include("build_model/battery/variables.jl")
include("build_model/battery/constraints_common.jl")
include("build_model/battery/constraints_dc.jl")
include("build_model/battery/constraints_ac.jl")
include("build_model/battery/constraints_socp.jl")
include("build_model/battery/dispatch.jl")

# --- PROSUMER
include("build_model/prosumer/variables.jl")
include("build_model/prosumer/constraints.jl")
include("build_model/prosumer/dispatch.jl")

# --- BALANCE
include("build_model/power_balance/active.jl")
include("build_model/power_balance/reactive.jl")
include("build_model/power_balance/dispatch.jl")

# --- OBJECTIVE
include("build_model/objective/generation_cost.jl")
include("build_model/objective/dispatch.jl")

# --- builder/solver
include("build_model/build_model.jl")
include("build_model/solve_model.jl")

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
# ------------------------------------------------------------------
include("api.jl")

# ------------------------------------------------------------------
# New function-based 
# ------------------------------------------------------------------
include("solve.jl")

# ------------------------------------------------------------------
# Public exports (API)
# ------------------------------------------------------------------
export solve, dc, ac, socp,
    centralized, admm,
    partition_from_case, partition_metis

end
