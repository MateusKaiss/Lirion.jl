# ============================================================
# Threading benchmark — distributed ADMM
# ============================================================
#
# Compares wall-clock time of distributed ADMM runs with
# BLAS=1 vs BLAS=physical-cores, holding JULIA_NUM_THREADS
# constant.
#
# The production `_with_scoped_blas_threads` always forces
# BLAS=1 inside `run_distributed`. For the benchmark we
# redefine it as a no-op so the externally-set BLAS count
# actually reaches the solver.
#
# Usage:
#
#   JULIA_NUM_THREADS=28 julia --project benchmark/threading_benchmark.jl
#
# ============================================================

using Lirion
using LinearAlgebra
using Printf
using Statistics


# ============================================================
# Monkey-patch the BLAS scoping wrapper (benchmark-only)
# ============================================================

Lirion.eval(quote
    function _with_scoped_blas_threads(f, n::Integer = 1)
        return f()
    end
end)


# ============================================================
# Configuration
# ============================================================

const PHYSICAL_CORES = 14
const TRIALS         = 3

const BENCH_CASES = [
    ("case57_ieee",     (; rho = 1e5, max_iter = 500)),
    ("case73_ieee_rts", (; rho = 1e5, max_iter = 500)),
]


# ============================================================
# Helpers
# ============================================================

_case_path(case::String) = joinpath(
    pkgdir(Lirion), "paper", "cases", case, "std", "$(case)_std_tie.m",
)

function _bench_solve(case::String, params)
    return Lirion.solve(
        _case_path(case);
        model     = dc(),
        algorithm = admm(;
            coupling     = :state,
            sectioning   = :tie_line,
            max_iter     = params.max_iter,
            rho          = params.rho,
            eps_abs      = 1e-4,
            eps_rel      = 1e-4,
            adaptive_rho = true,
            mu_incr      = 30.0,
            mu_decr      = 10.0,
            tau_incr     = 1.2,
            tau_decr     = 0.5,
        ),
        optimizer = nothing,
        outputs   = false,
        plots     = false,
    )
end

function _trials(f)
    times = Float64[]
    result = nothing
    for _ in 1:TRIALS
        push!(times, @elapsed (result = f()))
    end
    return (; med = median(times), min = minimum(times), max = maximum(times), result)
end


# ============================================================
# Main
# ============================================================

function main()
    @printf("\nJulia threads = %d  |  Machine CPU threads = %d\n\n",
            Threads.nthreads(), Sys.CPU_THREADS)

    @info "Warmup..."
    BLAS.set_num_threads(1);              _bench_solve("case57_ieee", (; rho = 1e5, max_iter = 50))
    BLAS.set_num_threads(PHYSICAL_CORES); _bench_solve("case57_ieee", (; rho = 1e5, max_iter = 50))
    @info "Warmup done"

    println()
    @printf("%-15s | %5s | %10s | %10s | %10s | %s\n",
            "case", "BLAS", "median(s)", "min(s)", "max(s)", "iters/status")
    println("-"^85)

    for (case, params) in BENCH_CASES
        @info "Running $case..."
        for blas in (1, PHYSICAL_CORES)
            BLAS.set_num_threads(blas)
            tr = _trials(() -> _bench_solve(case, params))
            iters  = try Int(tr.result.iter_done) catch; -1 end
            status = try string(tr.result.status) catch; "unknown" end
            @printf("%-15s | %5d | %10.2f | %10.2f | %10.2f | %d / %s\n",
                    case, blas, tr.med, tr.min, tr.max, iters, status)
        end
        println()
    end
end

main()
