using Test
using Lirion
using PowerModels
using Glob
using Ipopt
using SCS
using OSQP
using COSMO
using Gurobi

include("common.jl")

const DEFAULT_SOCP_SOLVER = pick_socp_optimizer()

include(joinpath(@__DIR__, "test_objective.jl"))
include(joinpath(@__DIR__, "test_timeseries.jl"))

include(joinpath(@__DIR__, "test_admm_dc.jl"))
include(joinpath(@__DIR__, "test_admm_ac.jl"))
include(joinpath(@__DIR__, "test_admm_socp.jl"))

@testset "Lirion - All Tests" begin
    # Centralizados
    run_test_objective()
    run_test_timeseries()

    # Distributed DC
    run_admm_dc_state_tests(; rtol = 0.01)
    run_admm_dc_power_tests(; rtol = 0.01)
    run_admm_dc_hybrid_tests(; rtol = 0.01)

    # # Distributed AC
    run_admm_ac_state_tests(;
        rtol = 0.01,
        central_optimizer = Ipopt.Optimizer,
        distributed_optimizer = Ipopt.Optimizer,
    )
    run_admm_ac_power_tests(;
        rtol = 0.01,
        central_optimizer = Ipopt.Optimizer,
        distributed_optimizer = Ipopt.Optimizer,
    )
    run_admm_ac_hybrid_tests(;
        rtol = 0.01,
        central_optimizer = Ipopt.Optimizer,
        distributed_optimizer = Ipopt.Optimizer,
    )

    # Distributed SOCP
    run_admm_socp_state_tests(;
        rtol = 0.01,
        central_optimizer = DEFAULT_SOCP_SOLVER,
        distributed_optimizer = DEFAULT_SOCP_SOLVER,
    )
    run_admm_socp_power_tests(;
        rtol = 0.01,
        central_optimizer = DEFAULT_SOCP_SOLVER,
        distributed_optimizer = DEFAULT_SOCP_SOLVER,
    )
    run_admm_socp_hybrid_tests(;
        rtol = 0.01,
        central_optimizer = DEFAULT_SOCP_SOLVER,
        distributed_optimizer = DEFAULT_SOCP_SOLVER,
    )
end

print_failures()
