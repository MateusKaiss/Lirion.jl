using Test
using PowerOptimization
using PowerModels
using Glob
using Gurobi
using Ipopt
using SCS

include("test_objective.jl")
include("test_timeseries.jl")
include("test_distributed_objective.jl")
include("test_distributed_branch_objective.jl")
include("test_distributed_hybrid_objective.jl")

@testset "PowerOptimization - All Tests" begin
    run_test_objective()
    run_test_timeseries()
    run_distributed_objective_tests(; rtol = 0.05)
    run_distributed_branch_objective_tests(; rtol = 0.05)
    run_distributed_hybrid_objective_tests(; rtol = 0.05)
end

println("\n===== GLOBAL SUMMARY OF FAILURES BY MODEL =====")
for (label, failures) in GLOBAL_FAILED_CASES
    println("\n[$label]")
    for (source, cases) in failures
        println("  $source error in $(length(cases)) case(s):")
        for c in cases
            println("    - $c")
        end
    end
end
println("==============================================\n")
