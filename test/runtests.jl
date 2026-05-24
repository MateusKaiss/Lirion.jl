using Test
using Lirion
using Glob
using Ipopt
using OSQP
using Clarabel
using JuMP: optimizer_with_attributes

# PowerModels is a test-only dependency (declared in [extras] / [targets], so
# `Pkg.test()` provides it). It is used by the multi-period tests (which
# cross-validate against PM's multi-network OPF), by the opt-in gold test
# (env var LIRION_TEST_VS_PM=1), and by the fixture regeneration script. The
# centralized fixture suite itself is frozen and PM-free.
const LIRION_TEST_VS_PM = get(ENV, "LIRION_TEST_VS_PM", "0") == "1"
if LIRION_TEST_VS_PM
    using PowerModels
end

include("common.jl")

const DEFAULT_SOCWR_SOLVER = optimizer_with_attributes(
    Clarabel.Optimizer,
    "verbose" => false,
    "tol_feas" => 1e-4,
    "tol_gap_abs" => 1e-4,
    "tol_gap_rel" => 1e-4,
    "max_iter" => 1200,
    "dynamic_regularization_enable" => true,
    "dynamic_regularization_eps" => 1e-10,
    "dynamic_regularization_delta" => 1e-4,
    "static_regularization_enable" => true,
    "equilibrate_enable" => true,
    "equilibrate_min_scaling" => 1e-4,
    "equilibrate_max_scaling" => 1e4,
    "direct_solve_method" => :qdldl,
    "iterative_refinement_reltol" => 1e-9,
    "iterative_refinement_abstol" => 1e-9,
    "max_step_fraction" => 0.7,
    "linesearch_backtrack_step" => 0.5,
    "tol_ktratio" => 1e-5,
    "min_terminate_step_length" => 1e-6,
    "min_switch_step_length" => 5e-2,
    "reduced_tol_gap_abs" => 5e-4,
    "reduced_tol_gap_rel" => 5e-4,
    "reduced_tol_feas" => 5e-4,
    "reduced_tol_infeas_abs" => 5e-11,
    "reduced_tol_infeas_rel" => 5e-4,
    "reduced_tol_ktratio" => 1e-4,
)

# Fast suite (default): frozen fixtures, no PM dependency.
include(joinpath(@__DIR__, "test_centralized_fixtures.jl"))

# Gold suite (opt-in): live PM roundtrip.
LIRION_TEST_VS_PM && include(joinpath(@__DIR__, "test_objective.jl"))

include(joinpath(@__DIR__, "test_multiperiod.jl"))
include(joinpath(@__DIR__, "test_renewables.jl"))

include(joinpath(@__DIR__, "distributed", "test_admm_dc.jl"))
include(joinpath(@__DIR__, "distributed", "test_admm_ac.jl"))
include(joinpath(@__DIR__, "distributed", "test_admm_socwr.jl"))
include(joinpath(@__DIR__, "distributed", "test_aladin_dc.jl"))
include(joinpath(@__DIR__, "distributed", "test_aladin_ac.jl"))
include(joinpath(@__DIR__, "distributed", "test_aladin_socwr.jl"))

@testset "Lirion - All Tests" begin
    # Centralized: fast fixtures (default) + optional PM roundtrip.
    run_centralized_fixture_tests(; rtol = 1e-3)
    LIRION_TEST_VS_PM && run_test_objective()

    # Multi-period invariants (internal consistency)
    run_multiperiod_invariants(; optimizer = Ipopt.Optimizer)

    # Renewables
    run_renewables_tests(; optimizer = Ipopt.Optimizer)

    # Distributed DC
    run_admm_dc_state_tests(; rtol = 0.01)
    run_admm_dc_power_tests(; rtol = 0.01)
    run_admm_dc_hybrid_tests(; rtol = 0.01)
    run_admm_dc_dummybus_hybrid_tests(; rtol = 0.05)
    run_admm_dc_dummybus_power_tests(; rtol = 0.05)
    run_admm_dc_splitbus_hybrid_tests(; rtol = 0.01)

    # Distributed AC
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
    run_admm_ac_splitbus_hybrid_tests(;
        rtol = 0.05,
        central_optimizer = Ipopt.Optimizer,
        distributed_optimizer = Ipopt.Optimizer,
    )
    run_admm_ac_dummybus_hybrid_tests(;
        rtol = 0.05,
        central_optimizer = Ipopt.Optimizer,
        distributed_optimizer = Ipopt.Optimizer,
    )
    run_admm_ac_dummybus_power_tests(;
        rtol = 0.05,
        central_optimizer = Ipopt.Optimizer,
        distributed_optimizer = Ipopt.Optimizer,
    )

    # Distributed SOCWR
    run_admm_socwr_state_tests(;
        rtol = 0.01,
        central_optimizer = DEFAULT_SOCWR_SOLVER,
        distributed_optimizer = DEFAULT_SOCWR_SOLVER,
    )
    run_admm_socwr_power_tests(;
        rtol = 0.01,
        central_optimizer = DEFAULT_SOCWR_SOLVER,
        distributed_optimizer = DEFAULT_SOCWR_SOLVER,
    )
    run_admm_socwr_hybrid_tests(;
        rtol = 0.01,
        central_optimizer = DEFAULT_SOCWR_SOLVER,
        distributed_optimizer = DEFAULT_SOCWR_SOLVER,
    )
    run_admm_socwr_splitbus_hybrid_tests(;
        rtol = 0.05,
        central_optimizer = DEFAULT_SOCWR_SOLVER,
        distributed_optimizer = DEFAULT_SOCWR_SOLVER,
    )
    run_admm_socwr_dummybus_hybrid_tests(;
        rtol = 0.05,
        central_optimizer = DEFAULT_SOCWR_SOLVER,
        distributed_optimizer = DEFAULT_SOCWR_SOLVER,
    )
    run_admm_socwr_dummybus_power_tests(;
        rtol = 0.05,
        central_optimizer = DEFAULT_SOCWR_SOLVER,
        distributed_optimizer = DEFAULT_SOCWR_SOLVER,
    )

    # Distributed ALADIN — DC
    run_aladin_dc_state_tests(; rtol = 0.01)
    run_aladin_dc_power_tests(; rtol = 0.01)
    run_aladin_dc_hybrid_tests(; rtol = 0.01)
    run_aladin_dc_hybrid_pflow_tests(; rtol = 0.01)
    run_aladin_dc_splitbus_hybrid_tests(; rtol = 0.05)
    run_aladin_dc_dummybus_hybrid_tests(; rtol = 0.05)
    run_aladin_dc_dummybus_power_tests(; rtol = 0.05)

    # Distributed ALADIN — AC
    run_aladin_ac_state_tests(;
        rtol = 0.01,
        central_optimizer = Ipopt.Optimizer,
        distributed_optimizer = Ipopt.Optimizer,
    )
    run_aladin_ac_power_tests(;
        rtol = 0.01,
        central_optimizer = Ipopt.Optimizer,
        distributed_optimizer = Ipopt.Optimizer,
    )
    run_aladin_ac_hybrid_tests(;
        rtol = 0.01,
        central_optimizer = Ipopt.Optimizer,
        distributed_optimizer = Ipopt.Optimizer,
    )
    run_aladin_ac_splitbus_hybrid_tests(;
        rtol = 0.05,
        central_optimizer = Ipopt.Optimizer,
        distributed_optimizer = Ipopt.Optimizer,
    )
    run_aladin_ac_dummybus_hybrid_tests(;
        rtol = 0.05,
        central_optimizer = Ipopt.Optimizer,
        distributed_optimizer = Ipopt.Optimizer,
    )
    run_aladin_ac_dummybus_power_tests(;
        rtol = 0.05,
        central_optimizer = Ipopt.Optimizer,
        distributed_optimizer = Ipopt.Optimizer,
    )

    # Distributed ALADIN — SOCWR
    run_aladin_socwr_state_tests(;
        rtol = 0.01,
        central_optimizer = DEFAULT_SOCWR_SOLVER,
        distributed_optimizer = DEFAULT_SOCWR_SOLVER,
    )
    run_aladin_socwr_power_tests(;
        rtol = 0.01,
        central_optimizer = DEFAULT_SOCWR_SOLVER,
        distributed_optimizer = DEFAULT_SOCWR_SOLVER,
    )
    run_aladin_socwr_hybrid_tests(;
        rtol = 0.01,
        central_optimizer = DEFAULT_SOCWR_SOLVER,
        distributed_optimizer = DEFAULT_SOCWR_SOLVER,
    )
    run_aladin_socwr_hybrid_pflow_tests(;
        rtol = 0.01,
        central_optimizer = DEFAULT_SOCWR_SOLVER,
        distributed_optimizer = DEFAULT_SOCWR_SOLVER,
    )
    run_aladin_socwr_splitbus_hybrid_tests(;
        rtol = 0.05,
        central_optimizer = DEFAULT_SOCWR_SOLVER,
        distributed_optimizer = DEFAULT_SOCWR_SOLVER,
    )
end

print_failures()
