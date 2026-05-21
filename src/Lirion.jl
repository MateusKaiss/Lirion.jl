module Lirion

using LinearAlgebra
using SparseArrays

# ------------------------------------------------------------------
# Packages
# ------------------------------------------------------------------
using JuMP
using Ipopt
using Dates
using Statistics
using Logging
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

"""
    Lirion.threadinfo()

Print current threading configuration — Julia threads, BLAS threads,
and detected CPU cores. Useful before a distributed run to sanity-check
that Julia has enough threads for the number of zones.
"""
function threadinfo()
    println("Julia threads : $(Threads.nthreads())")
    println("BLAS threads  : $(LinearAlgebra.BLAS.get_num_threads())")
    println("CPU threads   : $(Sys.CPU_THREADS)")
    return nothing
end

# ------------------------------------------------------------------
# Basic types
# ------------------------------------------------------------------
include("types.jl")

# ------------------------------------------------------------------
# Raw case representation (returned by parsers)
# ------------------------------------------------------------------
include("data/raw_structs.jl")

# ------------------------------------------------------------------
# Case parsing (.m / .csv)
# ------------------------------------------------------------------
include("parsing/matlab_parser_core.jl")
include("parsing/matpower_tables.jl")
include("parsing/case_parser.jl")

# ------------------------------------------------------------------
# Data layer: static data + time series + SystemData
# ------------------------------------------------------------------
include("data/timeseries_utils.jl")
include("data/timeseries/bus_timeseries.jl")
include("data/timeseries/gencost_timeseries.jl")

include("data/bus.jl")
include("data/branch.jl")
include("data/generator.jl")
include("data/storage.jl")
include("data/objective.jl")
include("data/power_balance.jl")
include("data/system_data.jl")
include("data/preprocess_case.jl")

# ------------------------------------------------------------------
# Partitioning / Sectioning
# ------------------------------------------------------------------
include("partitioning/partitioning_utils.jl")
include("partitioning/tie_line.jl")
include("partitioning/dummy_bus.jl")
include("partitioning/split_bus.jl")
include("partitioning/create_partition.jl")

# ------------------------------------------------------------------
# Inputs structs (legacy CLI + shared input pipeline)
# ------------------------------------------------------------------
include("data/build_inputs.jl")
include("data/partition_structs.jl")

# ------------------------------------------------------------------
# Model construction and centralized solution
# ------------------------------------------------------------------

# --- BUS
include("modeling/bus/variables.jl")
include("modeling/bus/dispatch.jl")

# --- BRANCH
include("modeling/branch/dc.jl")
include("modeling/branch/ac.jl")
include("modeling/branch/acr.jl")
include("modeling/branch/ivr.jl")
include("modeling/branch/lpac.jl")
include("modeling/branch/socwr_rt.jl")
include("modeling/branch/socwr.jl")
include("modeling/branch/qc.jl")
include("modeling/branch/sdpwrm.jl")
include("modeling/branch/socbf.jl")
include("modeling/branch/bfa.jl")
include("modeling/branch/nfa.jl")
include("modeling/branch/dispatch.jl")

# --- GENERATOR
include("modeling/generator/variables.jl")
include("modeling/generator/constraints.jl")
include("modeling/generator/dispatch.jl")

# --- BATTERY
include("modeling/storage/variables.jl")
include("modeling/storage/constraints_common.jl")
include("modeling/storage/constraints_dc.jl")
include("modeling/storage/constraints_ac.jl")
include("modeling/storage/constraints_acr.jl")
include("modeling/storage/constraints_socwr.jl")
include("modeling/storage/dispatch.jl")

# --- BALANCE
include("modeling/power_balance/active.jl")
include("modeling/power_balance/reactive.jl")
include("modeling/power_balance/ivr_balance.jl")
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
include("postprocessing/lirion_result.jl")
include("postprocessing/build_results.jl")
include("postprocessing/build_plots.jl")

# ------------------------------------------------------------------
# Distributed methods
# ------------------------------------------------------------------
include("distributed_algorithms/optimizer.jl")
include("distributed_algorithms/local_models.jl")
include("distributed_algorithms/shared/algorithms.jl")
include("distributed_algorithms/shared/traits.jl")
include("distributed_algorithms/shared/totals.jl")
include("distributed_algorithms/shared/sectioning_helpers.jl")

include("distributed_algorithms/aladin/params.jl")

include("distributed_algorithms/shared/history.jl")

include("distributed_algorithms/admm/core/var_state.jl")
include("distributed_algorithms/admm/core/snapshots.jl")
include("distributed_algorithms/admm/core/flow_accessors.jl")
include("distributed_algorithms/admm/params.jl")

include("distributed_algorithms/admm/dc/workspace.jl")
include("distributed_algorithms/admm/dc/partition.jl")
include("distributed_algorithms/admm/dc/objective.jl")

include("distributed_algorithms/admm/ac/workspace.jl")
include("distributed_algorithms/admm/ac/partition.jl")
include("distributed_algorithms/admm/ac/objective.jl")

include("distributed_algorithms/admm/socwr/workspace.jl")
include("distributed_algorithms/admm/socwr/partition.jl")
include("distributed_algorithms/admm/socwr/objective.jl")

include("distributed_algorithms/admm/core/save.jl")
include("distributed_algorithms/admm/core/consensus.jl")
include("distributed_algorithms/admm/core/dual_update.jl")
include("distributed_algorithms/admm/core/residuals.jl")
include("distributed_algorithms/admm/core/tolerances.jl")
include("distributed_algorithms/admm/core/rho_update.jl")

include("distributed_algorithms/aladin/core/zone_helpers.jl")
include("distributed_algorithms/aladin/core/coordinator_qp.jl")

include("distributed_algorithms/aladin/dc/zone_data.jl")
include("distributed_algorithms/aladin/dc/sensitivities.jl")
include("distributed_algorithms/aladin/dc/active_set.jl")
include("distributed_algorithms/aladin/dc/pairs.jl")
include("distributed_algorithms/aladin/dc/partition.jl")
include("distributed_algorithms/aladin/dc/workspace.jl")
include("distributed_algorithms/aladin/dc/objective.jl")
include("distributed_algorithms/aladin/dc/coordinator_qp.jl")
include("distributed_algorithms/aladin/dc/line_search.jl")

include("distributed_algorithms/aladin/ac/zone_data.jl")
include("distributed_algorithms/aladin/ac/sensitivities.jl")
include("distributed_algorithms/aladin/ac/active_set.jl")
include("distributed_algorithms/aladin/ac/pairs.jl")
include("distributed_algorithms/aladin/ac/partition.jl")
include("distributed_algorithms/aladin/ac/workspace.jl")
include("distributed_algorithms/aladin/ac/objective.jl")
include("distributed_algorithms/aladin/ac/coordinator_qp.jl")
include("distributed_algorithms/aladin/ac/line_search.jl")

include("distributed_algorithms/aladin/socwr/zone_data.jl")
include("distributed_algorithms/aladin/socwr/sensitivities.jl")
include("distributed_algorithms/aladin/socwr/active_set.jl")
include("distributed_algorithms/aladin/socwr/pairs.jl")
include("distributed_algorithms/aladin/socwr/partition.jl")
include("distributed_algorithms/aladin/socwr/workspace.jl")
include("distributed_algorithms/aladin/socwr/objective.jl")
include("distributed_algorithms/aladin/socwr/coordinator_qp.jl")
include("distributed_algorithms/aladin/socwr/line_search.jl")

include("distributed_algorithms/aladin/core/sensitivities.jl")
include("distributed_algorithms/aladin/core/flow_accessors.jl")

include("distributed_algorithms/aladin/core/primal_update.jl")
include("distributed_algorithms/aladin/core/dual_update.jl")
include("distributed_algorithms/aladin/core/residuals.jl")
include("distributed_algorithms/aladin/core/tolerances.jl")
include("distributed_algorithms/aladin/core/adapt_rho_mu.jl")

include("distributed_algorithms/shared/workspace_traits.jl")
include("distributed_algorithms/shared/local_solve.jl")
include("distributed_algorithms/shared/logging.jl")

include("distributed_algorithms/admm/run.jl")
include("distributed_algorithms/aladin/run.jl")

# ------------------------------------------------------------------
# Distributed OPF top-level runner (dispatches ADMM vs ALADIN)
# ------------------------------------------------------------------
include("distributed_algorithms/distributed_opf.jl")

# ------------------------------------------------------------------
# Legacy runner / CLI entrypoints
# ------------------------------------------------------------------
include("api.jl")

# ------------------------------------------------------------------
# New function-based API
# ------------------------------------------------------------------
include("solve.jl")

for T in (
    DCStateWorkspace, DCPowerWorkspace, DCHybridWorkspace,
    ACStateWorkspace, ACPowerWorkspace, ACHybridWorkspace,
    SOCWRStateWorkspace, SOCWRPowerWorkspace, SOCWRHybridWorkspace,
)
    precompile(run_algorithm!, (T, ADMMParams))
end

for T in (
    DCStateWorkspaceALADIN, DCPowerWorkspaceALADIN, DCHybridWorkspaceALADIN,
    ACStateWorkspaceALADIN, ACPowerWorkspaceALADIN, ACHybridWorkspaceALADIN,
    SOCWRStateWorkspaceALADIN, SOCWRPowerWorkspaceALADIN, SOCWRHybridWorkspaceALADIN,
)
    precompile(run_algorithm!, (T, ALADINParams))
end

# ------------------------------------------------------------------
# Public exports
# ------------------------------------------------------------------
export solve, load,
    dc, ac, acr, ivr, lpac, socwr, qc, sdpwrm, socbf, bfa, nfa,
    centralized, admm, aladin,
    create_partition, create_partition_split_bus, write_matpower,
    LirionResult, to_csvs

end
