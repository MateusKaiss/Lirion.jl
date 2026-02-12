using Test
using Lirion
using PowerModels
using Glob
using Gurobi
using Ipopt
using SCS
using OSQP
using COSMO   # <- adicionado para o fallback

# -----------------------------
# Helper: escolher SOCP solver
# -----------------------------
# Ordem padrão: Gurobi → COSMO → SCS. Você pode alterar via:
# ENV["LIRION_SOCP_PREFER"] = "auto" | "gurobi" | "cosmo" | "scs"
const _SOCP_CHAIN_MAP = Dict(
    :auto => (Gurobi.Optimizer, COSMO.Optimizer, SCS.Optimizer),
    :gurobi => (Gurobi.Optimizer, COSMO.Optimizer, SCS.Optimizer),
    :cosmo => (COSMO.Optimizer, Gurobi.Optimizer, SCS.Optimizer),
    :scs => (SCS.Optimizer, Gurobi.Optimizer, COSMO.Optimizer),
)

function pick_socp_optimizer(;
    prefer::Symbol = Symbol(
        lowercase(get(ENV, "LIRION_SOCP_PREFER", "auto")),
    ),
)
    chain = get(_SOCP_CHAIN_MAP, prefer, _SOCP_CHAIN_MAP[:auto])
    for factory in chain
        try
            _ = factory()             # testa instanciar (pega licença do Gurobi, etc.)
            return factory
        catch
            # tenta o próximo
        end
    end
    return SCS.Optimizer               # fallback final paranoico
end

# Use o mesmo solver para centralizado e distribuído (pode diferenciar se quiser)
const DEFAULT_SOCP_SOLVER = pick_socp_optimizer()

include("test_objective.jl")
include("test_timeseries.jl")
include("test_distributed_objective.jl")
include("test_distributed_branch_objective.jl")
include("test_distributed_hybrid_objective.jl")
include("test_distributed_objective_ac_admm.jl")
include("test_distributed_branch_objective_ac.jl")
include("test_distributed_hybrid_objective_ac.jl")
include("test_distributed_socp_objective.jl")
include("test_distributed_branch_objective_socp.jl")
include("test_distributed_hybrid_objective_socp.jl")

@testset "Lirion - All Tests" begin
    run_test_objective()
    run_test_timeseries()

    run_distributed_objective_tests(; rtol = 0.01)
    run_distributed_branch_objective_tests(; rtol = 0.01)
    run_distributed_hybrid_objective_tests(; rtol = 0.01)


    run_distributed_ac_objective_tests(; rtol = 0.01)
    run_distributed_branch_ac_objective_tests(; rtol = 0.01)
    run_distributed_hybrid_ac_objective_tests(; rtol = 0.01)


    run_distributed_socp_objective_tests(;
        rtol = 0.01,
        central_optimizer = DEFAULT_SOCP_SOLVER,
        distributed_optimizer = DEFAULT_SOCP_SOLVER,
    )

    run_distributed_branch_socp_objective_tests(;
        rtol = 0.01,
        central_optimizer = DEFAULT_SOCP_SOLVER,
        distributed_optimizer = DEFAULT_SOCP_SOLVER,
    )
    
    run_distributed_hybrid_socp_objective_tests(;
        rtol = 0.01,
        central_optimizer = DEFAULT_SOCP_SOLVER,
        distributed_optimizer = DEFAULT_SOCP_SOLVER,
    )
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
