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
# data: static data + time series + SystemData
# ------------------------------------------------------------------
include("data/timeseries_utils.jl")
include("data/timeseries/bus_timeseries.jl")
include("data/timeseries/prosumers_timeseries.jl")
include("data/timeseries/consumer_timeseries.jl")
include("data/timeseries/gencost_timeseries.jl")

include("data/bus.jl")
include("data/branch.jl")
include("data/generator.jl")
include("data/battery.jl")
include("data/prosumers.jl")
include("data/objective.jl")
include("data/power_balance.jl")
include("data/system_data.jl")
include("data/preprocess_case.jl")

# ------------------------------------------------------------------
# Partitioning / Sectioning
# ------------------------------------------------------------------
include("partitioning/partitioning_utils.jl")
include("partitioning/tie_line.jl")

# ------------------------------------------------------------------
# Inputs structs (legacy CLI + shared input pipeline)
# ------------------------------------------------------------------
include("data/build_inputs.jl")
include("data/admm_structs.jl")

# ------------------------------------------------------------------
# Case parsing (.m / .csv)
# ------------------------------------------------------------------
include("parsing/parser_utils.jl")
include("parsing/matpower_parser.jl")
include("parsing/csv_parser.jl")
include("parsing/parse_case.jl")

# ------------------------------------------------------------------
# Model construction and centralized solution
# ------------------------------------------------------------------
# --- BUS
include("modeling/bus/variables.jl")
include("modeling/bus/dispatch.jl")

# --- BRANCH
include("modeling/branch/dc.jl")
include("modeling/branch/ac.jl")
include("modeling/branch/socp_rt.jl")
include("modeling/branch/socp.jl")
include("modeling/branch/dispatch.jl")

# --- GENERATOR
include("modeling/generator/variables.jl")
include("modeling/generator/constraints.jl")
include("modeling/generator/dispatch.jl")

# --- BATTERY
include("modeling/battery/variables.jl")
include("modeling/battery/constraints_common.jl")
include("modeling/battery/constraints_dc.jl")
include("modeling/battery/constraints_ac.jl")
include("modeling/battery/constraints_socp.jl")
include("modeling/battery/dispatch.jl")

# --- PROSUMER
include("modeling/prosumer/variables.jl")
include("modeling/prosumer/constraints.jl")
include("modeling/prosumer/dispatch.jl")

# --- BALANCE
include("modeling/power_balance/active.jl")
include("modeling/power_balance/reactive.jl")
include("modeling/power_balance/dispatch.jl")

# --- OBJECTIVE
include("modeling/objective/generation_cost.jl")
include("modeling/objective/dispatch.jl")

# --- builder/solver
include("modeling/build_model.jl")
include("modeling/solve_model.jl")

# ------------------------------------------------------------------
# Post-processing / results
# ------------------------------------------------------------------
include("postprocessing/build_results.jl")
include("postprocessing/build_plots.jl")

# ------------------------------------------------------------------
# Distributed methods
# ------------------------------------------------------------------
include("distributed_algorithms/distributed_opf.jl")

include("distributed_algorithms/admm/admm_support.jl")
include("distributed_algorithms/admm/buffers.jl")
include("distributed_algorithms/admm/objective_terms.jl")
include("distributed_algorithms/admm/consensus.jl")
include("distributed_algorithms/admm/duals.jl")
include("distributed_algorithms/admm/residuals.jl")
include("distributed_algorithms/admm/rho.jl")

include("distributed_algorithms/admm/run_admm/admm_dc.jl")
include("distributed_algorithms/admm/run_admm/admm_ac.jl")
include("distributed_algorithms/admm/run_admm/admm_socp.jl")

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
