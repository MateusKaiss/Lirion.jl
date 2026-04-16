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
using SparseArrays
try
    LinearAlgebra.BLAS.set_num_threads(1)
catch e
    @debug "BLAS.set_num_threads(1) not supported by this backend" exception = e
end

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
include("data/battery.jl")
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
include("modeling/branch/socbf.jl")
include("modeling/branch/bfa.jl")
include("modeling/branch/nfa.jl")
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
include("modeling/battery/constraints_acr.jl")
include("modeling/battery/constraints_socwr.jl")
include("modeling/battery/dispatch.jl")

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
include("postprocessing/build_results.jl")
include("postprocessing/build_plots.jl")

# ------------------------------------------------------------------
# Distributed methods
# ------------------------------------------------------------------
include("distributed_algorithms/optimizer.jl")
include("distributed_algorithms/partitioning.jl")
include("distributed_algorithms/shared/algorithms.jl")
include("distributed_algorithms/shared/traits.jl")
include("distributed_algorithms/shared/totals.jl")
include("distributed_algorithms/shared/logging.jl")
include("distributed_algorithms/shared/local_solve.jl")

include("distributed_algorithms/admm/core/support.jl")
include("distributed_algorithms/admm/core/components.jl")
include("distributed_algorithms/admm/core/history.jl")
include("distributed_algorithms/admm/core/snapshots.jl")
include("distributed_algorithms/admm/core/flow_accessors.jl")
include("distributed_algorithms/admm/params.jl")
include("distributed_algorithms/admm/core/local_solve.jl")

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
include("distributed_algorithms/admm/core/logging.jl")
include("distributed_algorithms/admm/run.jl")

include("distributed_algorithms/aladin/core/support.jl")
include("distributed_algorithms/aladin/core/history.jl")
include("distributed_algorithms/aladin/params.jl")
include("distributed_algorithms/aladin/core/local_solve.jl")

include("distributed_algorithms/aladin/dc/sensitivities.jl")
include("distributed_algorithms/aladin/dc/workspace.jl")
include("distributed_algorithms/aladin/dc/objective.jl")
include("distributed_algorithms/aladin/dc/coordinator_qp.jl")

include("distributed_algorithms/aladin/ac/sensitivities.jl")
include("distributed_algorithms/aladin/ac/workspace.jl")
include("distributed_algorithms/aladin/ac/objective.jl")
include("distributed_algorithms/aladin/ac/coordinator_qp.jl")

include("distributed_algorithms/aladin/socwr/sensitivities.jl")
include("distributed_algorithms/aladin/socwr/workspace.jl")
include("distributed_algorithms/aladin/socwr/objective.jl")
include("distributed_algorithms/aladin/socwr/coordinator_qp.jl")
include("distributed_algorithms/aladin/socwr/line_search.jl")

include("distributed_algorithms/aladin/core/flow_accessors.jl")

include("distributed_algorithms/aladin/core/primal_update.jl")
include("distributed_algorithms/aladin/core/dual_update.jl")
include("distributed_algorithms/aladin/core/residuals.jl")
include("distributed_algorithms/aladin/core/tolerances.jl")
include("distributed_algorithms/aladin/core/adapt_rho_mu.jl")
include("distributed_algorithms/aladin/core/logging.jl")

include("distributed_algorithms/shared/workspace_traits.jl")

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

# ------------------------------------------------------------------
# Precompile hints — compile the ADMM / ALADIN main loops for every
# workspace type at module-load time so the first call of each
# `{algorithm, formulation, coupling}` combination doesn't pay JIT
# overhead. Native Julia `precompile` directive, no extra dependency.
# ------------------------------------------------------------------
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
export solve,
    dc, ac, acr, ivr, lpac, socwr, socbf, bfa, nfa,
    centralized, admm, aladin,
    partition_from_case, partition_metis

end
