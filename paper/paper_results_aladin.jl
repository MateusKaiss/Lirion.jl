# ============================================================
# paper_results_aladin.jl
# ============================================================
# Distributed ALADIN benchmark for the Lirion paper.
#
# Same sweep shape as paper_results_admm.jl but with ALADIN. Reads the
# centralised reference objectives from `paper/cases/_centralized.toml`.
#
# ALADIN carries more knobs than ADMM (ρ, μ, σ_{theta,pg,other}, α1/2/3,
# rho_growth, mu_growth, rho_upper, mu_upper, adaptive). We start from
# the Houska-Frasch-Diehl canonical defaults and let per-(case, variant,
# base, sectioning, coupling) overrides in `ALADIN_OVERRIDES` supersede
# them. The overrides are populated by `paper/tune_case14.jl` and
# subsequent per-case tuning runs.
#
# Writes `paper/cases/_aladin.toml` incrementally.
#
# Usage:
#     julia> include(joinpath(pkgdir(Lirion), "paper", "paper_results_aladin.jl"))
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

# ALADIN's coordinator QP is ill-conditioned at large ρ; we need
# tight solver tolerances to keep it bounded. Lirion's tuned
# silent-factories carry those settings — pass bare types (or nothing
# for DC) so the resolver in `src/distributed_algorithms/optimizer.jl`
# picks them up instead of a loosely-configured user factory.
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

const PAPER_DIR    = @__DIR__
const CASES_DIR    = joinpath(PAPER_DIR, "cases")
const REPORT_CEN   = joinpath(CASES_DIR, "_centralized.toml")
const REPORT_ALAD  = joinpath(CASES_DIR, "_aladin.toml")

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
# ALADIN defaults (Houska-Frasch-Diehl canonical values, matching
# Lirion's ALADINConfig constructor in `src/types.jl`).
# ------------------------------------------------------------
const ALADIN_DEFAULTS = (
    max_iter           = 200,
    rho                = 1e3,
    mu                 = 1e4,
    alpha1             = 1.0,
    alpha2             = 1.0,
    alpha3             = 1.0,
    sigma_theta_shared = 10.0,
    sigma_theta_other  = 10.0,
    sigma_pg           = 1.0,
    sigma_other        = 1.0,
    active_set_tol     = 1e-6,
    adaptive           = true,
    rho_growth         = 1.3,
    mu_growth          = 2.0,
    rho_upper          = 1e6,
    mu_upper           = 1e6,
)

# Key: (case, variant, base, sectioning, coupling) → Dict of kwargs
# overrides. Any absent key falls through to ALADIN_DEFAULTS.
const ALADIN_OVERRIDES = Dict{Tuple{String, String, Symbol, Symbol, Symbol}, Dict{Symbol, Any}}(
    # --- case14_ieee (from _tune_case14_aladin_v2.toml) --------------
    ("case14_ieee", "std", :dc, :tie_line, :state) => Dict(:rho => 1e+03, :mu => 1e+03, :adaptive => true, :rho_growth => 1.05, :mu_growth => 1.1),  # gentle_adaptive gap=+0.001%, iters=9
    ("case14_ieee", "std", :dc, :tie_line, :power) => Dict(:rho => 1e+03, :mu => 1e+03, :adaptive => true, :rho_growth => 1.05, :mu_growth => 1.1),  # gentle_adaptive gap=+0.001%, iters=3
    ("case14_ieee", "std", :dc, :tie_line, :hybrid) => Dict(:rho => 1e+03, :mu => 1e+03, :adaptive => true, :rho_growth => 1.05, :mu_growth => 1.1),  # gentle_adaptive gap=+0.001%, iters=3
    ("case14_ieee", "std", :dc, :dummy_bus, :power) => Dict(:rho => 1e+03, :mu => 1e+03, :adaptive => true, :rho_growth => 1.05, :mu_growth => 1.1),  # gentle_adaptive gap=+0.001%, iters=3
    ("case14_ieee", "std", :dc, :dummy_bus, :hybrid) => Dict(:rho => 1e+03, :mu => 1e+03, :adaptive => true, :rho_growth => 1.05, :mu_growth => 1.1),  # gentle_adaptive gap=+0.001%, iters=7
    ("case14_ieee", "std", :dc, :split_bus, :hybrid) => Dict(:rho => 1e+03, :mu => 1e+03, :adaptive => true, :rho_growth => 1.05, :mu_growth => 1.1),  # gentle_adaptive gap=+0.001%, iters=3
    ("case14_ieee", "std", :ac, :tie_line, :state) => Dict(:rho => 1e+03, :mu => 1e+03, :adaptive => true, :rho_growth => 1.05, :mu_growth => 1.1),  # gentle_adaptive gap=+0.294%, iters=405
    ("case14_ieee", "std", :ac, :tie_line, :power) => Dict(:rho => 1e+03, :mu => 1e+03, :adaptive => true, :rho_growth => 1.05, :mu_growth => 1.1),  # gentle_adaptive gap=+0.019%, iters=43
    ("case14_ieee", "std", :ac, :tie_line, :hybrid) => Dict(:rho => 1e+03, :mu => 1e+03, :adaptive => true, :rho_growth => 1.05, :mu_growth => 1.1),  # gentle_adaptive gap=+0.044%, iters=55
    ("case14_ieee", "std", :ac, :dummy_bus, :power) => Dict(:rho => 1e+04, :mu => 1e+03, :adaptive => true, :rho_growth => 1.05, :mu_growth => 1.1),  # gentle_hi_rho gap=+0.584%, iters=75
    ("case14_ieee", "std", :ac, :dummy_bus, :hybrid) => Dict(:rho => 1e+04, :mu => 1e+03, :adaptive => true, :rho_growth => 1.05, :mu_growth => 1.1),  # gentle_hi_rho gap=+0.478%, iters=121
    ("case14_ieee", "std", :ac, :split_bus, :hybrid) => Dict(:rho => 1e+03, :mu => 1e+03, :adaptive => true, :rho_growth => 1.05, :mu_growth => 1.1),  # gentle_adaptive gap=+0.112%, iters=59
    ("case14_ieee", "std", :socwr, :tie_line, :state) => Dict(:rho => 100, :mu => 100, :adaptive => false),  # fixed_mid gap=-0.228%, iters=612
    ("case14_ieee", "std", :socwr, :tie_line, :power) => Dict(:rho => 100, :mu => 10, :adaptive => false),  # fixed_1e2_mu1 gap=-0.001%, iters=293
    ("case14_ieee", "std", :socwr, :tie_line, :hybrid) => Dict(:rho => 100, :mu => 10, :adaptive => false),  # fixed_1e2_mu1 gap=+0.008%, iters=347
    ("case14_ieee", "std", :socwr, :split_bus, :hybrid) => Dict(:rho => 100, :mu => 100, :adaptive => false),  # fixed_mid gap=-0.005%, iters=436
    ("case14_ieee", "api", :dc, :tie_line, :hybrid) => Dict(:rho => 1e+03, :mu => 1e+03, :adaptive => true, :rho_growth => 1.05, :mu_growth => 1.1),  # gentle_adaptive gap=-0.927%, iters=27
    ("case14_ieee", "api", :dc, :dummy_bus, :hybrid) => Dict(:rho => 1e+03, :mu => 1e+03, :adaptive => true, :rho_growth => 1.05, :mu_growth => 1.1),  # gentle_adaptive gap=-0.925%, iters=29
    ("case14_ieee", "api", :dc, :split_bus, :hybrid) => Dict(:rho => 1e+03, :mu => 1e+03, :adaptive => true, :rho_growth => 1.05, :mu_growth => 1.1),  # gentle_adaptive gap=-0.926%, iters=28
    ("case14_ieee", "api", :ac, :tie_line, :state) => Dict(:rho => 1e+03, :mu => 1e+03, :adaptive => true, :rho_growth => 1.05, :mu_growth => 1.1),  # gentle_adaptive gap=-0.040%, iters=86
    ("case14_ieee", "api", :ac, :tie_line, :power) => Dict(:rho => 1e+03, :mu => 1e+03, :adaptive => false),  # fixed_1e3 gap=-0.113%, iters=580
    ("case14_ieee", "api", :ac, :tie_line, :hybrid) => Dict(:rho => 1e+03, :mu => 1e+03, :adaptive => false),  # fixed_1e3 gap=-0.040%, iters=867
    ("case14_ieee", "api", :ac, :split_bus, :hybrid) => Dict(:rho => 1e+03, :mu => 1e+03, :adaptive => true, :rho_growth => 1.05, :mu_growth => 1.1),  # gentle_adaptive gap=-0.047%, iters=60
    ("case14_ieee", "api", :socwr, :tie_line, :state) => Dict(:rho => 1e+03, :mu => 1e+03, :adaptive => false),  # fixed_1e3 gap=+0.017%, iters=411
    ("case14_ieee", "api", :socwr, :tie_line, :power) => Dict(:rho => 100, :mu => 10, :adaptive => false),  # fixed_1e2_mu1 gap=-0.078%, iters=985
    ("case14_ieee", "api", :socwr, :tie_line, :hybrid) => Dict(:rho => 100, :mu => 100, :adaptive => false),  # fixed_mid gap=-0.007%, iters=130
    ("case14_ieee", "api", :socwr, :split_bus, :hybrid) => Dict(:rho => 100, :mu => 100, :adaptive => false),  # fixed_mid gap=-0.003%, iters=160
    ("case14_ieee", "sad", :ac, :tie_line, :state) => Dict(:rho => 1e+04, :mu => 1e+03, :adaptive => true, :rho_growth => 1.05, :mu_growth => 1.1),  # gentle_hi_rho gap=-0.426%, iters=534
    ("case14_ieee", "sad", :ac, :tie_line, :hybrid) => Dict(:rho => 1e+04, :mu => 1e+03, :adaptive => true, :rho_growth => 1.05, :mu_growth => 1.1),  # gentle_hi_rho gap=+0.198%, iters=126
    ("case14_ieee", "sad", :ac, :dummy_bus, :power) => Dict(:rho => 1e+04, :mu => 1e+03, :adaptive => true, :rho_growth => 1.05, :mu_growth => 1.1),  # gentle_hi_rho gap=-0.787%, iters=112
    ("case14_ieee", "sad", :ac, :dummy_bus, :hybrid) => Dict(:rho => 1e+04, :mu => 1e+03, :adaptive => true, :rho_growth => 1.05, :mu_growth => 1.1),  # gentle_hi_rho gap=+0.484%, iters=148
    ("case14_ieee", "sad", :ac, :split_bus, :hybrid) => Dict(:rho => 1e+03, :mu => 1e+03, :adaptive => true, :rho_growth => 1.05, :mu_growth => 1.1),  # gentle_adaptive gap=+0.023%, iters=61
    ("case14_ieee", "sad", :socwr, :tie_line, :state) => Dict(:rho => 100, :mu => 100, :adaptive => false),  # fixed_mid gap=-0.261%, iters=838
    ("case14_ieee", "sad", :socwr, :tie_line, :power) => Dict(:rho => 100, :mu => 10, :adaptive => false),  # fixed_1e2_mu1 gap=-0.013%, iters=357
    ("case14_ieee", "sad", :socwr, :tie_line, :hybrid) => Dict(:rho => 100, :mu => 10, :adaptive => false),  # fixed_1e2_mu1 gap=-0.002%, iters=479
    ("case14_ieee", "sad", :socwr, :split_bus, :hybrid) => Dict(:rho => 100, :mu => 100, :adaptive => false),  # fixed_mid gap=+0.007%, iters=381
)

function _aladin_kwargs(case::String, variant::String,
                        base::Symbol, sectioning::Symbol, coupling::Symbol)
    kw = Dict{Symbol, Any}(pairs(ALADIN_DEFAULTS))
    ov = get(ALADIN_OVERRIDES, (case, variant, base, sectioning, coupling), nothing)
    ov === nothing && return kw
    for (k, v) in ov
        kw[k] = v
    end
    return kw
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

function _ref_objective(cen::Dict, case::String, variant::String, base::Symbol, base_ctor)
    rec = get(get(get(cen, case, Dict()), variant, Dict()), string(base), Dict())
    if get(rec, "status", "") == "converged" &&
       isfinite(get(rec, "objective", NaN)) && get(rec, "objective", 0.0) > 0
        return float(rec["objective"])
    end
    base === :ac    && !_HAVE_IPOPT                     && return NaN
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
    kw   = _aladin_kwargs(case, variant, base, sectioning, coupling)

    t0 = time()
    status = "error"; obj = NaN; iters = 0
    try
        alg = aladin(;
            coupling = coupling, sectioning = sectioning,
            eps_abs = 1e-4, eps_rel = 1e-4,
            kw...,
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
        "rho"        => kw[:rho],
        "mu"         => kw[:mu],
        "max_iter"   => kw[:max_iter],
        "adaptive"   => kw[:adaptive],
        "rho_growth" => kw[:rho_growth],
        "mu_growth"  => kw[:mu_growth],
        "status"     => status,
        "objective"  => isfinite(obj) ? obj : 0.0,
        "gap_pct"    => gap,
        "iters"      => iters,
        "time_s"     => dt,
    )
end

function main()
    println("╔════════════════════════════════════════════════════════════╗")
    println("║  Lirion paper — DISTRIBUTED ALADIN                         ║")
    println("║  Cases:    ", join(CASES, ", "))
    println("║  Variants: ", join(VARIANTS, ", "))
    @printf("║  Solvers:  Gurobi=%s  Ipopt=%s  Clarabel=%s  OSQP=%s\n",
            _HAVE_GUROBI, _HAVE_IPOPT, _HAVE_CLARABEL, _HAVE_OSQP)
    println("║  Started:  ", now())
    println("║  Output:   ", REPORT_ALAD)
    println("╚════════════════════════════════════════════════════════════╝")
    flush(stdout)

    isdir(CASES_DIR) || mkpath(CASES_DIR)
    cen    = isfile(REPORT_CEN)  ? TOML.parsefile(REPORT_CEN)  : Dict{String, Any}()
    report = isfile(REPORT_ALAD) ? TOML.parsefile(REPORT_ALAD) : Dict{String, Any}()

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
                _save_toml!(report, REPORT_ALAD)
                gap_str = isfinite(res["gap_pct"]) ?
                          @sprintf("%+6.3f%%", res["gap_pct"]) : "  NaN "
                @printf("[%s/%-3s/%-5s/%-9s/%-6s] ρ=%8.0e μ=%8.0e iters=%5d gap=%s t=%6.2fs %s%s\n",
                        case, variant, base, sectioning, coupling,
                        res["rho"], res["mu"], res["iters"], gap_str,
                        res["time_s"], res["status"],
                        res["adaptive"] ? " (adaptive)" : "")
                flush(stdout)
            end
        end
    end

    report["_meta"] = Dict(
        "generated_at" => string(now()),
        "script"       => "paper_results_aladin.jl",
        "defaults"     => Dict(string(k) => string(v) for (k, v) in pairs(ALADIN_DEFAULTS)),
        "eps_abs" => 1e-4, "eps_rel" => 1e-4,
    )
    _save_toml!(report, REPORT_ALAD)
    println("\nALADIN results written to: ", REPORT_ALAD)
    return report
end

main()
