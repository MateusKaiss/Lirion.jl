# ============================================================
# paper_results_admm.jl
# ============================================================
# Distributed ADMM benchmark for the Lirion paper.
#
# Reads centralised reference objectives from `paper/cases/_centralized.toml`
# (produced by `paper_results_centralized.jl`). If absent, computes the
# reference on-the-fly per (case, variant, base).
#
# Sweep: CASES × VARIANTS × BASE_MODELS × ALGO_VARIANTS, with per-(case,
# base, coupling) ρ / max_iter overrides populated by `tune_case14.jl`
# and subsequent per-case tuning runs.
#
# Writes `paper/cases/_admm.toml` incrementally after every solve so a
# kill does not lose progress.
#
# Usage:
#     julia> include(joinpath(pkgdir(Lirion), "paper", "paper_results_admm.jl"))
# ============================================================

using Lirion
using JuMP: optimizer_with_attributes
using Printf
using Dates
using TOML

const _HAVE_GUROBI   = try; eval(:(using Gurobi));   true; catch; false; end
const _HAVE_IPOPT    = try; eval(:(using Ipopt));    true; catch; false; end
const _HAVE_OSQP     = try; eval(:(using OSQP));     true; catch; false; end
const _HAVE_CLARABEL = try; eval(:(using Clarabel)); true; catch; false; end

_HAVE_IPOPT || @warn "Ipopt not loaded — AC results will be unavailable."

const TIME_LIMIT_S = 300.0

# Route so Lirion's tuned silent factories kick in (see
# `src/distributed_algorithms/optimizer.jl`): bare types `Ipopt.Optimizer`
# / `Clarabel.Optimizer` trigger the `_silent_factory(::Type{...})`
# methods carrying the tolerance/adaptive settings. `nothing` for DC
# lets Lirion fall through to its hand-tuned OSQP attributes required
# to keep subproblem QPs bounded under large ρ.
_gurobi_opt() = optimizer_with_attributes(
    () -> Gurobi.Optimizer(Gurobi.Env()),
    "OutputFlag" => 0, "TimeLimit" => TIME_LIMIT_S,
    "FeasibilityTol" => 1e-7, "OptimalityTol" => 1e-7,
)

_solver_for_base(base::Symbol) = begin
    base === :dc    && return _HAVE_GUROBI ? _gurobi_opt() : nothing
    base === :ac    && return _HAVE_IPOPT  ? Ipopt.Optimizer : nothing
    base === :socwr && return _HAVE_GUROBI ? _gurobi_opt() :
                              (_HAVE_CLARABEL ? Clarabel.Optimizer : nothing)
    error("unknown base $base")
end

const PAPER_DIR   = @__DIR__
const CASES_DIR   = joinpath(PAPER_DIR, "cases")
const REPORT_CEN  = joinpath(CASES_DIR, "_centralized.toml")
const REPORT_ADMM = joinpath(CASES_DIR, "_admm.toml")

const CASES    = ["case14_ieee", "case30_ieee", "case57_ieee",
                  "case73_ieee_rts", "case118_ieee"]
const VARIANTS = ["std", "api", "sad"]

function _base_specs()
    specs = [(:dc, dc), (:socwr, socwr)]
    _HAVE_IPOPT && insert!(specs, 2, (:ac, ac))
    return specs
end

const ALGO_VARIANTS = [
    (:tie_line,  :state),
    (:tie_line,  :power),
    (:tie_line,  :hybrid),
    (:dummy_bus, :power),
    (:dummy_bus, :hybrid),
    (:split_bus, :hybrid),
]

# ------------------------------------------------------------
# Per-(case, variant, base, sectioning, coupling) ADMM overrides.
#
# Populated by `paper/tune_case14.jl` and equivalent per-case tune
# scripts. Each override records: ρ, max_iter, adaptive_rho, and
# (optionally) the residual-balancing knobs μ_incr/decr, τ_incr/decr.
#
# Falls back to `DEFAULT_RHO` / `DEFAULT_MAX_ITER` when no override is
# present. Keys use the tuple form for O(1) lookup.
# ------------------------------------------------------------
const DEFAULT_RHO::Dict{Symbol, Float64} = Dict(:state => 1e5, :power => 1e3, :hybrid => 1e3)
const DEFAULT_MAX_ITER::Int = 5000

const ADMM_OVERRIDES = Dict{Tuple{String, String, Symbol, Symbol, Symbol}, Dict{Symbol, Any}}(
    # --- case14_ieee (auto-generated from _tune_case14_admm.toml) -----
    ("case14_ieee", "std", :dc, :tie_line, :state) => Dict(:rho => 1e+05),  # gap=-0.304%, iters=1217
    ("case14_ieee", "std", :dc, :tie_line, :power) => Dict(:rho => 1e+03),  # gap=+0.018%, iters=34
    ("case14_ieee", "std", :dc, :tie_line, :hybrid) => Dict(:rho => 1e+03),  # gap=+0.003%, iters=83
    ("case14_ieee", "std", :dc, :dummy_bus, :power) => Dict(:rho => 1e+03),  # gap=-0.022%, iters=26
    ("case14_ieee", "std", :dc, :dummy_bus, :hybrid) => Dict(:rho => 1e+03),  # gap=-0.000%, iters=170
    ("case14_ieee", "std", :dc, :split_bus, :hybrid) => Dict(:rho => 1e+03),  # gap=-0.003%, iters=229
    ("case14_ieee", "std", :ac, :tie_line, :state) => Dict(:rho => 1e+05),  # gap=+0.133%, iters=434
    ("case14_ieee", "std", :ac, :tie_line, :power) => Dict(:rho => 1e+03),  # gap=-0.007%, iters=60
    ("case14_ieee", "std", :ac, :tie_line, :hybrid) => Dict(:rho => 1e+03),  # gap=-0.004%, iters=65
    ("case14_ieee", "std", :ac, :dummy_bus, :power) => Dict(:rho => 1e+03),  # gap=-0.045%, iters=222
    ("case14_ieee", "std", :ac, :dummy_bus, :hybrid) => Dict(:rho => 1e+03),  # gap=-0.006%, iters=177
    ("case14_ieee", "std", :ac, :split_bus, :hybrid) => Dict(:rho => 1e+03),  # gap=+0.003%, iters=194
    ("case14_ieee", "std", :socwr, :tie_line, :state) => Dict(:rho => 1e+05),  # gap=+0.019%, iters=2771
    ("case14_ieee", "std", :socwr, :tie_line, :power) => Dict(:rho => 1e+03),  # gap=-0.004%, iters=161
    ("case14_ieee", "std", :socwr, :tie_line, :hybrid) => Dict(:rho => 1e+03),  # gap=+0.007%, iters=208
    ("case14_ieee", "std", :socwr, :dummy_bus, :power) => Dict(:rho => 1e+03),  # gap=-0.021%, iters=659
    ("case14_ieee", "std", :socwr, :dummy_bus, :hybrid) => Dict(:rho => 1e+03),  # gap=+0.008%, iters=96
    ("case14_ieee", "std", :socwr, :split_bus, :hybrid) => Dict(:rho => 1e+03),  # gap=+0.001%, iters=125
    ("case14_ieee", "api", :dc, :tie_line, :state) => Dict(:rho => 1e+04),  # gap=-0.206%, iters=1234
    ("case14_ieee", "api", :dc, :tie_line, :hybrid) => Dict(:rho => 1e+03),  # gap=-0.918%, iters=186
    ("case14_ieee", "api", :dc, :dummy_bus, :hybrid) => Dict(:rho => 1e+03),  # gap=-0.922%, iters=255
    ("case14_ieee", "api", :dc, :split_bus, :hybrid) => Dict(:rho => 1e+03),  # gap=-0.913%, iters=276
    ("case14_ieee", "api", :ac, :tie_line, :state) => Dict(:rho => 1e+05),  # gap=-0.071%, iters=100
    ("case14_ieee", "api", :ac, :tie_line, :power) => Dict(:rho => 1e+03),  # gap=-0.180%, iters=1037
    ("case14_ieee", "api", :ac, :tie_line, :hybrid) => Dict(:rho => 1e+03),  # gap=-0.040%, iters=80
    ("case14_ieee", "api", :ac, :dummy_bus, :hybrid) => Dict(:rho => 1e+03),  # gap=-0.241%, iters=130
    ("case14_ieee", "api", :ac, :split_bus, :hybrid) => Dict(:rho => 1e+03),  # gap=-0.059%, iters=151
    ("case14_ieee", "api", :socwr, :tie_line, :state) => Dict(:rho => 1e+05),  # gap=+0.009%, iters=1448
    ("case14_ieee", "api", :socwr, :tie_line, :power) => Dict(:rho => 1e+03),  # gap=-0.052%, iters=99
    ("case14_ieee", "api", :socwr, :tie_line, :hybrid) => Dict(:rho => 1e+03),  # gap=+0.003%, iters=120
    ("case14_ieee", "api", :socwr, :dummy_bus, :power) => Dict(:rho => 1e+03),  # gap=-0.347%, iters=102
    ("case14_ieee", "api", :socwr, :dummy_bus, :hybrid) => Dict(:rho => 1e+03),  # gap=-0.092%, iters=63
    ("case14_ieee", "api", :socwr, :split_bus, :hybrid) => Dict(:rho => 1e+03),  # gap=-0.001%, iters=106
    ("case14_ieee", "sad", :ac, :tie_line, :state) => Dict(:rho => 1e+05),  # gap=+0.815%, iters=124
    ("case14_ieee", "sad", :ac, :tie_line, :hybrid) => Dict(:rho => 1e+03),  # gap=-0.016%, iters=113
    ("case14_ieee", "sad", :ac, :dummy_bus, :hybrid) => Dict(:rho => 1e+03),  # gap=+0.117%, iters=163
    ("case14_ieee", "sad", :ac, :split_bus, :hybrid) => Dict(:rho => 1e+03),  # gap=+0.019%, iters=225
    ("case14_ieee", "sad", :socwr, :tie_line, :state) => Dict(:rho => 1e+05),  # gap=+0.006%, iters=2729
    ("case14_ieee", "sad", :socwr, :tie_line, :power) => Dict(:rho => 1e+03),  # gap=-0.009%, iters=61
    ("case14_ieee", "sad", :socwr, :tie_line, :hybrid) => Dict(:rho => 1e+03),  # gap=+0.000%, iters=50
    ("case14_ieee", "sad", :socwr, :dummy_bus, :power) => Dict(:rho => 1e+03),  # gap=-0.033%, iters=94
    ("case14_ieee", "sad", :socwr, :dummy_bus, :hybrid) => Dict(:rho => 1e+03),  # gap=-0.012%, iters=88
    ("case14_ieee", "sad", :socwr, :split_bus, :hybrid) => Dict(:rho => 1e+03),  # gap=+0.001%, iters=117
)

function _admm_params(case::String, variant::String,
                      base::Symbol, sectioning::Symbol, coupling::Symbol)
    ov = get(ADMM_OVERRIDES, (case, variant, base, sectioning, coupling), nothing)
    rho          = ov === nothing ? DEFAULT_RHO[coupling]     : get(ov, :rho, DEFAULT_RHO[coupling])
    max_iter     = ov === nothing ? DEFAULT_MAX_ITER          : get(ov, :max_iter, DEFAULT_MAX_ITER)
    adaptive_rho = ov === nothing ? false                     : get(ov, :adaptive_rho, false)
    mu_incr      = ov === nothing ? 10.0                      : get(ov, :mu_incr, 10.0)
    mu_decr      = ov === nothing ? 10.0                      : get(ov, :mu_decr, 10.0)
    tau_incr     = ov === nothing ? 2.0                       : get(ov, :tau_incr, 2.0)
    tau_decr     = ov === nothing ? 0.5                       : get(ov, :tau_decr, 0.5)
    return (; rho, max_iter, adaptive_rho, mu_incr, mu_decr, tau_incr, tau_decr)
end

_partition_file(case::String, variant::String, sectioning::Symbol) = begin
    dir = joinpath(CASES_DIR, case, variant)
    sectioning === :split_bus ?
        joinpath(dir, "$(case)_$(variant)_splitbus.m") :
        joinpath(dir, "$(case)_$(variant)_tie.m")
end

_original_file(case::String, variant::String) =
    joinpath(CASES_DIR, case, variant, "$(case)_$(variant).m")

function _classify(s::Symbol)
    s in (:optimal, :locally_solved, :local_optimal, :almost_optimal) && return "converged"
    s in (:iteration_limit, :time_limit, :slow_progress, :numerical_error) &&
        return "did_not_converge"
    s in (:infeasible, :locally_infeasible, :dual_infeasible) && return "infeasible"
    return String(s)
end

function _save_toml!(report::Dict, path::AbstractString)
    tmp = path * ".tmp"
    open(tmp, "w") do io
        TOML.print(io, report; sorted = true)
    end
    mv(tmp, path; force = true)
    return nothing
end

# Try the pre-computed centralised report first; else compute inline.
function _ref_objective(cen::Dict, case::String, variant::String, base::Symbol, base_ctor)
    rec = get(get(get(cen, case, Dict()), variant, Dict()), string(base), Dict())
    if get(rec, "status", "") == "converged" &&
       isfinite(get(rec, "objective", NaN)) && get(rec, "objective", 0.0) > 0
        return float(rec["objective"])
    end
    # Fresh centralised solve. For DC, `opt === nothing` means "use
    # Lirion's tuned OSQP default" (not "no solver available"); for
    # AC/SOCWR it means the solver package isn't loaded and we must
    # skip.
    base === :ac    && !_HAVE_IPOPT                   && return NaN
    base === :socwr && !_HAVE_GUROBI && !_HAVE_CLARABEL && return NaN
    src = _original_file(case, variant)
    isfile(src) || return NaN
    opt = _solver_for_base(base)
    try
        r = opt === nothing ?
            Lirion.solve(src; model = base_ctor(), algorithm = centralized(),
                outputs = false, plots = false) :
            Lirion.solve(src; model = base_ctor(), algorithm = centralized(),
                optimizer = opt, outputs = false, plots = false)
        (_classify(r.status) == "converged" && isfinite(r.objective) && r.objective > 0) ||
            return NaN
        return float(r.objective)
    catch
        return NaN
    end
end

function _run_one(case::String, variant::String,
                  base::Symbol, base_ctor,
                  sectioning::Symbol, coupling::Symbol, ref_obj::Float64)
    path = _partition_file(case, variant, sectioning)
    p = _admm_params(case, variant, base, sectioning, coupling)

    t0 = time()
    status = "error"; obj = NaN; iters = 0
    try
        alg = admm(;
            coupling = coupling, sectioning = sectioning,
            max_iter = p.max_iter, rho = p.rho,
            eps_abs = 1e-4, eps_rel = 1e-4,
            adaptive_rho = p.adaptive_rho,
            mu_incr = p.mu_incr, mu_decr = p.mu_decr,
            tau_incr = p.tau_incr, tau_decr = p.tau_decr,
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
    return Dict{String, Any}(
        "rho"          => p.rho,
        "max_iter"     => p.max_iter,
        "adaptive_rho" => p.adaptive_rho,
        "status"       => status,
        "objective"    => isfinite(obj) ? obj : 0.0,
        "gap_pct"      => gap,
        "iters"        => iters,
        "time_s"       => dt,
    )
end

function main()
    println("╔════════════════════════════════════════════════════════════╗")
    println("║  Lirion paper — DISTRIBUTED ADMM                           ║")
    println("║  Cases:    ", join(CASES, ", "))
    println("║  Variants: ", join(VARIANTS, ", "))
    @printf("║  Solvers:  Gurobi=%s  Ipopt=%s  Clarabel=%s  OSQP=%s\n",
            _HAVE_GUROBI, _HAVE_IPOPT, _HAVE_CLARABEL, _HAVE_OSQP)
    println("║  Started:  ", now())
    println("║  Output:   ", REPORT_ADMM)
    println("╚════════════════════════════════════════════════════════════╝")
    flush(stdout)

    isdir(CASES_DIR) || mkpath(CASES_DIR)
    cen    = isfile(REPORT_CEN)  ? TOML.parsefile(REPORT_CEN)  : Dict{String, Any}()
    report = isfile(REPORT_ADMM) ? TOML.parsefile(REPORT_ADMM) : Dict{String, Any}()

    for case in CASES, variant in VARIANTS
        haskey(report, case) || (report[case] = Dict{String, Any}())
        haskey(report[case], variant) || (report[case][variant] = Dict{String, Any}())
        for (base, base_ctor) in _base_specs()
            bkey = string(base)
            haskey(report[case][variant], bkey) ||
                (report[case][variant][bkey] = Dict{String, Any}())
            ref_obj = _ref_objective(cen, case, variant, base, base_ctor)
            if !isfinite(ref_obj)
                @printf("[%s/%-3s/%-5s] SKIP (no centralised reference)\n",
                        case, variant, base)
                flush(stdout)
                continue
            end
            for (sectioning, coupling) in ALGO_VARIANTS
                key = "$(sectioning)_$(coupling)"
                if haskey(report[case][variant][bkey], key) &&
                   haskey(report[case][variant][bkey][key], "status")
                    rec = report[case][variant][bkey][key]
                    @printf("[%s/%-3s/%-5s/%-9s/%-6s] cached  gap=%s iters=%d\n",
                            case, variant, base, sectioning, coupling,
                            isfinite(get(rec, "gap_pct", NaN)) ?
                                @sprintf("%+6.3f%%", rec["gap_pct"]) : "NaN",
                            get(rec, "iters", 0))
                    flush(stdout)
                    continue
                end
                res = _run_one(case, variant, base, base_ctor,
                               sectioning, coupling, ref_obj)
                report[case][variant][bkey][key] = res
                _save_toml!(report, REPORT_ADMM)
                gap_str = isfinite(res["gap_pct"]) ?
                          @sprintf("%+6.3f%%", res["gap_pct"]) : "  NaN "
                @printf("[%s/%-3s/%-5s/%-9s/%-6s] ρ=%8.0e iters=%5d gap=%s t=%6.2fs %s%s\n",
                        case, variant, base, sectioning, coupling,
                        res["rho"], res["iters"], gap_str,
                        res["time_s"], res["status"],
                        res["adaptive_rho"] ? " (adaptive)" : "")
                flush(stdout)
            end
        end
    end

    report["_meta"] = Dict(
        "generated_at" => string(now()),
        "script"       => "paper_results_admm.jl",
        "default_rho"  => Dict(string(k) => v for (k, v) in DEFAULT_RHO),
        "default_max_iter" => DEFAULT_MAX_ITER,
        "eps_abs" => 1e-4, "eps_rel" => 1e-4,
    )
    _save_toml!(report, REPORT_ADMM)
    println("\nADMM results written to: ", REPORT_ADMM)
    return report
end

main()
