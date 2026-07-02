# ============================================================
# paper/scripts/experiments/run_c0b_threads.jl
# ============================================================
# Stage C0b — Julia threading alignment for distributed runs.
# DIAGNOSTIC ONLY (not tuning, not the main matrix). Measures whether
# JULIA_NUM_THREADS>1 speeds up distributed solves while BLAS is fixed
# at 1. Distributed ADMM/ALADIN solve zone subproblems via
# Threads.@threads, so the natural parallelism is ACROSS zones — capped
# by k. case118 k8 (8 zones) is included to reveal the benefit beyond
# the k2 cases.
#
# JULIA_NUM_THREADS is a startup setting, so this script is launched
# once per thread count (julia -t N) and writes a tagged per-N CSV; a
# separate merge step aggregates them.
#
# Official distributed solvers (C0): DC=OSQP, SOCWR=Clarabel(hardened),
# AC=Ipopt+MA57 via LIRION_HSLLIB. BLAS threads = 1. Fixed, safe rho
# (no tuning, no rho=0): DC rho=1e4, SOCWR/AC rho=1e3.
# ============================================================

include(joinpath(@__DIR__, "lib", "tuning_common.jl"))

using LinearAlgebra
BLAS.set_num_threads(1)

const NT = Threads.nthreads()
const RES_RTA = joinpath(RESULTS_DIR, "runtime_alignment")

rho_for(f) = f == "dc" ? 1e4 : 1e3

# (case, k, formulation, algorithm)
const CELLS = [
    ("case57_ieee", 2, "dc", "admm"),
    ("case57_ieee", 2, "socwr", "admm"),
    ("case57_ieee", 2, "ac", "admm"),
    ("case73_ieee_rts", 2, "dc", "admm"),
    ("case73_ieee_rts", 2, "ac", "admm"),
    ("case118_ieee", 2, "dc", "admm"),
    ("case118_ieee", 2, "ac", "admm"),
    ("case118_ieee", 8, "dc", "admm"),   # 8 zones -> reveals zone-parallelism
    ("case57_ieee", 2, "dc", "aladin"),
]

const C0B_COLUMNS = [
    :julia_threads, :blas_threads, :case, :k, :formulation, :algorithm,
    :sectioning, :coupling, :status, :iterations, :primal_residual, :gap_pct,
    :runtime_s, :objective, :warnings,
]

function run_one(case, k, f, algo; t1)
    path = partitioned_case_path(case, "std", k)
    cen = get(t1, (case, f), NaN)
    opt = official_dist_optimizer(f)
    status = "error"; obj = NaN; iters = 0; pr = NaN; rt = NaN; warn = ""
    alg = algo == "admm" ?
        Lirion.admm(; coupling = :hybrid, sectioning = :tie_line, max_iter = 200,
                    rho = rho_for(f), eps_abs = 1e-4, eps_rel = 1e-4, adaptive_rho = false) :
        Lirion.aladin(; coupling = :hybrid, sectioning = :tie_line, max_iter = 60,
                      rho = 1e3, mu = 1e3, rho_growth = 1.05)
    t0 = time()
    try
        r = Lirion.solve(path; model = model_ctor(f)(), algorithm = alg,
                         optimizer = opt, outputs = false, plots = false)
        rt = time() - t0
        obj = Float64(r.objective); status = string(r.status)
        iters = r.iter_done === nothing ? 0 : r.iter_done
        pr, _ = extract_residuals(r)
    catch e
        rt = time() - t0; warn = first(split(sprint(showerror, e), "\n"))
    end
    gap = (isnan(obj) || isnan(cen) || cen == 0) ? NaN : 100 * (obj - cen) / abs(cen)
    return (julia_threads = NT, blas_threads = BLAS.get_num_threads(),
            case = case, k = k, formulation = f, algorithm = algo,
            sectioning = "tie_line", coupling = "hybrid", status = status,
            iterations = iters, primal_residual = pr, gap_pct = gap,
            runtime_s = rt, objective = obj, warnings = warn)
end

function run_c0b()
    mkpath(RES_RTA)
    t1 = load_t1_objectives()
    @info "C0b thread run" julia_threads = NT blas_threads = BLAS.get_num_threads() ma57 = have_hsl()
    # discarded warm-up (compile the distributed paths)
    for (case, k, f, algo) in CELLS
        try; run_one(case, k, f, algo; t1 = t1); catch; end
    end
    rows = NamedTuple[]
    for (case, k, f, algo) in CELLS
        row = run_one(case, k, f, algo; t1 = t1)
        push!(rows, row)
        @info @sprintf("  T=%-2d %-14s k%d %-6s %-6s -> %-16s it=%-4d rt=%.3fs gap=%s%s",
            NT, case, k, f, algo, row.status, row.iterations, row.runtime_s,
            isnan(row.gap_pct) ? "NaN" : @sprintf("%.3f%%", row.gap_pct),
            isempty(row.warnings) ? "" : "  WARN")
    end
    out = write_csv(joinpath(RES_RTA, "c0b_t$(NT).csv"), rows; columns = C0B_COLUMNS)
    @info "Wrote" out = relpath(out, PAPER_DIR)
    return rows
end

run_c0b()
