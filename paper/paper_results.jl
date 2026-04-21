# ============================================================
# paper_results.jl
# ============================================================
# Reproduction script for the Lirion paper benchmarks.
#
# Usage (after `using Lirion`):
#
#     julia> include(joinpath(pkgdir(Lirion), "paper", "paper_results.jl"))
#
# The script:
#   (1) solves the centralised DC / AC / SOCWR OPF for each of the three
#       benchmark cases (case14_ieee, case30_ieee, case57_ieee — std
#       variant), using Gurobi when available, else Ipopt / Clarabel.
#   (2) solves the ADMM distributed OPF over the three couplings
#       (state / power / hybrid) combined with the three sectioning
#       methods (tie_line, dummy_bus, split_bus) — 6 algorithmic
#       variants per (case × base model) = 54 distributed runs total.
#   (3) reports iteration count, wall time, and gap vs. centralised.
#
# The partition files in `paper/cases/<case>/std/` were pre-generated
# with Metis (k = 2 zones) and committed so that reviewers reproduce
# the identical partition without re-running the MetisExt stage.
#
# This file lives under `paper/` which is NOT part of Lirion's loaded
# source tree — installing the package with `Pkg.add("Lirion")` does
# not pull this script in. Users who clone the repository get it
# alongside the partition files.
# ============================================================

using Lirion
using JuMP: optimizer_with_attributes
using Printf
using Dates

# ------------------------------------------------------------
# Solver detection — prefer commercial when licensed, otherwise
# open-source fallbacks. Determines which benchmark timings are
# comparable to the numbers reported in the paper.
# ------------------------------------------------------------
const _HAVE_GUROBI   = try; eval(:(using Gurobi));   true; catch; false; end
const _HAVE_IPOPT    = try; eval(:(using Ipopt));    true; catch; false; end
const _HAVE_OSQP     = try; eval(:(using OSQP));     true; catch; false; end
const _HAVE_CLARABEL = try; eval(:(using Clarabel)); true; catch; false; end

_HAVE_IPOPT || @warn "Ipopt not loaded — AC results will be unavailable."

const TIME_LIMIT_S = 300.0

_gurobi_opt() = optimizer_with_attributes(
    () -> Gurobi.Optimizer(Gurobi.Env()),
    "OutputFlag"     => 0,
    "TimeLimit"      => TIME_LIMIT_S,
    "FeasibilityTol" => 1e-7,
    "OptimalityTol"  => 1e-7,
)
_ipopt_opt() = optimizer_with_attributes(
    Ipopt.Optimizer, "print_level" => 0, "sb" => "yes",
    "max_cpu_time" => TIME_LIMIT_S,
)
_osqp_opt() = optimizer_with_attributes(
    OSQP.Optimizer, "verbose" => false,
    "time_limit" => TIME_LIMIT_S, "max_iter" => 50_000,
)
_clarabel_opt() = optimizer_with_attributes(
    Clarabel.Optimizer, "verbose" => false,
    "tol_feas" => 1e-5, "tol_gap_abs" => 1e-5, "tol_gap_rel" => 1e-5,
    "max_iter" => 4_000, "time_limit" => TIME_LIMIT_S,
)

# Choose best available solver for each problem class. LP/QP: Gurobi then
# OSQP; NLP: Ipopt only (no open alternative for the AC NLP here); SOCP:
# Gurobi (SOC barrier) then Clarabel.
_lp_solver()   = _HAVE_GUROBI ? _gurobi_opt() : (_HAVE_OSQP     ? _osqp_opt()     : nothing)
_nlp_solver()  = _HAVE_IPOPT  ? _ipopt_opt() : nothing
_socp_solver() = _HAVE_GUROBI ? _gurobi_opt() : (_HAVE_CLARABEL ? _clarabel_opt() : nothing)

_solver_for_base(base::Symbol) =
    base === :dc    ? _lp_solver()   :
    base === :ac    ? _nlp_solver()  :
    base === :socwr ? _socp_solver() :
    error("unknown base $base")

# ------------------------------------------------------------
# Per-(case, base, coupling) ρ / max_iter table.
#
# Values selected empirically on the paper benchmark — see the study's
# Section § on convergence behaviour for a detailed discussion. Briefly:
#
#   - `state` coupling propagates only θ (+ v_mag for AC) and requires
#     a larger ρ to enforce primal consensus. We use ρ = 1e5 on 14-bus
#     and ρ = 1e6 on 30-bus and 57-bus.
#   - `power` coupling without θ suffers from angle drift at the zone
#     interface. On case30 the gap plateaus at ≈ −1.2 % until ρ reaches
#     1e8 (DC only) with max_iter = 10 000. AC/SOCWR cannot absorb such
#     a large ρ (solver ill-conditioning) so they stay at the default.
#   - `hybrid` = state + power couplings combined → converges at the
#     default ρ = 1e3 on every case we tested.
# ------------------------------------------------------------
const DEFAULT_RHO = Dict(:state => 1e5, :power => 1e3, :hybrid => 1e3)
const DEFAULT_MAX_ITER = 5000

const CASE_OVERRIDES = Dict{Tuple{String, Symbol, Symbol}, Dict{Symbol, Any}}(
    # case14_ieee — state on any base needs ρ = 1e5
    ("case14_ieee", :dc,    :state)  => Dict(:rho => 1e5),
    ("case14_ieee", :ac,    :state)  => Dict(:rho => 1e5),
    ("case14_ieee", :socwr, :state)  => Dict(:rho => 1e5),

    # case30_ieee — state needs ρ = 1e6; DC power requires ρ = 1e8 + 10k iters
    ("case30_ieee", :dc,    :state)  => Dict(:rho => 1e6),
    ("case30_ieee", :ac,    :state)  => Dict(:rho => 1e6),
    ("case30_ieee", :socwr, :state)  => Dict(:rho => 1e6),
    ("case30_ieee", :dc,    :power)  => Dict(:rho => 1e8, :max_iter => 10_000),

    # case57_ieee — analogous to case30 for state
    ("case57_ieee", :dc,    :state)  => Dict(:rho => 1e6),
    ("case57_ieee", :ac,    :state)  => Dict(:rho => 1e6),
    ("case57_ieee", :socwr, :state)  => Dict(:rho => 1e6),
)

_rho_maxiter(case::String, base::Symbol, coupling::Symbol) = begin
    ov = get(CASE_OVERRIDES, (case, base, coupling), nothing)
    ov === nothing && return DEFAULT_RHO[coupling], DEFAULT_MAX_ITER
    get(ov, :rho, DEFAULT_RHO[coupling]), get(ov, :max_iter, DEFAULT_MAX_ITER)
end

# ------------------------------------------------------------
# Benchmark sweep
# ------------------------------------------------------------
const PAPER_DIR = @__DIR__
const CASES = ["case14_ieee", "case30_ieee", "case57_ieee"]

# (base, model_ctor). `ac` is skipped if Ipopt is not available since
# Lirion's AC formulation is a non-convex NLP.
function _base_specs()
    specs = [(:dc, dc), (:socwr, socwr)]
    _HAVE_IPOPT && insert!(specs, 2, (:ac, ac))
    return specs
end

const COUPLINGS = [:state, :power, :hybrid]

# TieLine and DummyBus share the Metis k-way partition file; SplitBus
# uses the bilateral-constrained partition (guarantees each bus borders
# ≤ 1 foreign zone).
const ALGO_VARIANTS = [
    (:tie_line,  :state),
    (:tie_line,  :power),
    (:tie_line,  :hybrid),
    (:dummy_bus, :power),
    (:dummy_bus, :hybrid),
    (:split_bus, :hybrid),
]

_partition_file(case::String, sectioning::Symbol) = begin
    dir = joinpath(PAPER_DIR, "cases", case, "std")
    if sectioning === :split_bus
        joinpath(dir, "$(case)_std_splitbus.m")
    else
        joinpath(dir, "$(case)_std_tie.m")
    end
end

function _classify(s::Symbol)
    s in (:optimal, :locally_solved, :local_optimal, :almost_optimal) && return "converged"
    s in (:iteration_limit, :time_limit, :slow_progress, :numerical_error) &&
        return "did_not_converge"
    s in (:infeasible, :locally_infeasible, :dual_infeasible) && return "infeasible"
    return String(s)
end

# A single distributed run. Returns a NamedTuple that the driver prints
# into the final table. Any exception is caught so one bad combo does
# not abort the whole reproduction run.
function _run_one(case::String, base::Symbol, base_ctor,
                  sectioning::Symbol, coupling::Symbol, ref_obj::Float64)
    path = _partition_file(case, sectioning)
    rho, max_iter = _rho_maxiter(case, base, coupling)
    t0 = time()
    status = "error"; obj = NaN; iters = 0
    try
        alg = admm(;
            coupling = coupling, sectioning = sectioning,
            max_iter = max_iter, rho = rho,
            eps_abs = 1e-4, eps_rel = 1e-4, adaptive_rho = false,
        )
        r = Lirion.solve(path;
            model = base_ctor(), algorithm = alg,
            optimizer = _solver_for_base(base),
            outputs = false, plots = false,
        )
        status = _classify(r.status)
        obj    = isfinite(r.objective) ? r.objective : NaN
        iters  = isa(r.iter_done, Integer) ? Int(r.iter_done) : 0
    catch err
        status = "error"
    end
    dt = time() - t0
    gap = (isfinite(obj) && isfinite(ref_obj) && ref_obj > 0) ?
          100.0 * (obj - ref_obj) / ref_obj : NaN
    (case = case, base = base, sectioning = sectioning, coupling = coupling,
     rho = rho, max_iter = max_iter, status = status,
     objective = obj, gap_pct = gap, iters = iters, time_s = dt)
end

function main()
    println("╔════════════════════════════════════════════════════════════╗")
    println("║  Lirion paper benchmark reproduction                       ║")
    println("║  Cases: ", join(CASES, ", "), " (std)")
    @printf("║  Solvers: Gurobi=%s  Ipopt=%s  Clarabel=%s  OSQP=%s\n",
            _HAVE_GUROBI, _HAVE_IPOPT, _HAVE_CLARABEL, _HAVE_OSQP)
    println("║  Started: ", now())
    println("╚════════════════════════════════════════════════════════════╝")

    rows = NamedTuple[]
    for case in CASES
        original = joinpath(PAPER_DIR, "cases", case, "std", "$(case)_std.m")

        # Centralised references — one per base model. ADMM gaps are
        # reported relative to these.
        refs = Dict{Symbol, Float64}()
        for (base, base_ctor) in _base_specs()
            opt = _solver_for_base(base)
            opt === nothing && continue
            r = Lirion.solve(original; model = base_ctor(), algorithm = centralized(),
                optimizer = opt, outputs = false, plots = false)
            refs[base] = r.objective
            @printf("[%s / centralised / %-5s]  obj = %12.2f   t = %.2fs\n",
                    case, base, r.objective, 0.0)
        end

        for (base, base_ctor) in _base_specs()
            haskey(refs, base) || continue
            for (sectioning, coupling) in ALGO_VARIANTS
                res = _run_one(case, base, base_ctor, sectioning, coupling, refs[base])
                push!(rows, res)
                gap_str = isfinite(res.gap_pct) ?
                          @sprintf("%+6.3f %%", res.gap_pct) : "   NaN  "
                @printf("[%s / %-5s / %-9s / %-6s]  ρ=%8.0e  iters=%5d  gap=%s  t=%6.2fs  %s\n",
                        case, base, sectioning, coupling, res.rho, res.iters,
                        gap_str, res.time_s, res.status)
            end
        end
    end

    # Flatten into a compact summary table.
    println("\n" * "="^80)
    println(" SUMMARY")
    println("="^80)
    @printf("%-18s %-5s %-9s %-6s  %8s %5s %8s %7s\n",
            "case", "base", "sect", "coup", "rho", "iters", "gap%", "t(s)")
    for r in rows
        @printf("%-18s %-5s %-9s %-6s  %8.0e %5d %+8.3f %7.2f\n",
                r.case, r.base, r.sectioning, r.coupling,
                r.rho, r.iters, r.gap_pct, r.time_s)
    end

    return rows
end

# Entering the script via `include` executes the benchmark. Users who
# want the driver (e.g. to customise solvers) can `import` this file
# as a module instead by wrapping with `module Paper … end`.
main()
