# ============================================================
# tune_case14_aladin_v2.jl
# ============================================================
# A focused ALADIN tuning sweep for case14_ieee that addresses the
# shortcomings of tune_case14.jl's one-shot ρ sweep. ALADIN has 4
# load-bearing knobs (ρ, μ, rho_growth, mu_growth) plus an on/off
# adaptive switch; we explore a small hand-picked grid covering the
# three strategies that empirically work on case14:
#
#   gentle_adaptive — rho=mu=1e3, growth=1.05/1.1, adaptive=true
#   aggressive_adap — rho=mu=1e4, growth=1.3/2.0, adaptive=true (ALADIN default)
#   fixed_mid       — rho=mu=1e2..1e4, adaptive=false
#
# For each (variant, base, sectioning, coupling) we try each strategy
# and keep the first that converges with gap < 1 %.
#
# Output: `paper/cases/_tune_case14_aladin_v2.toml` and updates to
# `ALADIN_OVERRIDES` via `_apply_tune_overrides.jl`.
# ============================================================

using Lirion
using JuMP: optimizer_with_attributes
using Printf
using Dates
using TOML

const _HAVE_GUROBI   = try; eval(:(using Gurobi));   true; catch; false; end
const _HAVE_IPOPT    = try; eval(:(using Ipopt));    true; catch; false; end
const _HAVE_CLARABEL = try; eval(:(using Clarabel)); true; catch; false; end

const TIME_LIMIT_S = 120.0
const GAP_TARGET   = 1.0
const EPS_ABS      = 1e-4
const EPS_REL      = 1e-4
const MAX_ITER     = 20_000

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

_base_available(base::Symbol) = begin
    base === :dc && return true
    base === :ac && return _HAVE_IPOPT
    base === :socwr && return _HAVE_GUROBI || _HAVE_CLARABEL
    false
end

const PAPER_DIR  = @__DIR__
const CASES_DIR  = joinpath(PAPER_DIR, "cases")
# `LIRION_TUNE_CASE` env var chooses which PGLib-OPF case to tune.
const CASE       = get(ENV, "LIRION_TUNE_CASE", "case14_ieee")
const OUT        = joinpath(CASES_DIR, "_tune_$(CASE)_aladin_v2.toml")

const VARIANTS = ["std", "api", "sad"]
const ALGO_VARIANTS = [
    (:tie_line,  :state),  (:tie_line,  :power),  (:tie_line,  :hybrid),
    (:dummy_bus, :power),  (:dummy_bus, :hybrid), (:split_bus, :hybrid),
]

function _base_specs()
    specs = [(:dc, dc), (:socwr, socwr)]
    _HAVE_IPOPT && insert!(specs, 2, (:ac, ac))
    return specs
end

# Hand-picked strategies (name, kwargs-to-aladin). Order matters — first
# strategy that converges with gap<1% wins for that combo.
const STRATEGIES = [
    ("gentle_adaptive", (rho = 1e3, mu = 1e3,
        rho_growth = 1.05, mu_growth = 1.1, adaptive = true)),
    ("gentle_hi_rho",  (rho = 1e4, mu = 1e3,
        rho_growth = 1.05, mu_growth = 1.1, adaptive = true)),
    ("aggressive_adap", (rho = 1e3, mu = 1e4,
        rho_growth = 1.3, mu_growth = 2.0, adaptive = true)),
    ("fixed_mid",       (rho = 1e2, mu = 1e2, adaptive = false)),
    ("fixed_1e3",       (rho = 1e3, mu = 1e3, adaptive = false)),
    ("fixed_1e4",       (rho = 1e4, mu = 1e3, adaptive = false)),
    ("fixed_1e2_mu1",   (rho = 1e2, mu = 1e1, adaptive = false)),
]

_partition_file(variant::String, sectioning::Symbol) = begin
    dir = joinpath(CASES_DIR, CASE, variant)
    sectioning === :split_bus ?
        joinpath(dir, "$(CASE)_$(variant)_splitbus.m") :
        joinpath(dir, "$(CASE)_$(variant)_tie.m")
end

_original_file(variant::String) =
    joinpath(CASES_DIR, CASE, variant, "$(CASE)_$(variant).m")

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

const _REF_CACHE = Dict{Tuple{String, Symbol}, Float64}()

function _central_reference(variant::String, base::Symbol, base_ctor)
    key = (variant, base)
    haskey(_REF_CACHE, key) && return _REF_CACHE[key]
    _base_available(base) || (return _REF_CACHE[key] = NaN)
    src = _original_file(variant)
    isfile(src) || (return _REF_CACHE[key] = NaN)
    opt = _solver_for_base(base)
    try
        r = opt === nothing ?
            Lirion.solve(src; model = base_ctor(), algorithm = centralized(),
                outputs = false, plots = false) :
            Lirion.solve(src; model = base_ctor(), algorithm = centralized(),
                optimizer = opt, outputs = false, plots = false)
        if _classify(r.status) == "converged" && isfinite(r.objective) && r.objective > 0
            return _REF_CACHE[key] = float(r.objective)
        end
    catch
    end
    return _REF_CACHE[key] = NaN
end

function _probe(variant, base, base_ctor, sect, coup, strategy_kw, ref)
    path = _partition_file(variant, sect)
    t0 = time()
    status, obj, iters = "error", NaN, 0
    try
        alg = aladin(; coupling = coup, sectioning = sect, max_iter = MAX_ITER,
            eps_abs = EPS_ABS, eps_rel = EPS_REL, strategy_kw...)
        opt = _solver_for_base(base)
        r = opt === nothing ?
            Lirion.solve(path; model = base_ctor(), algorithm = alg,
                outputs = false, plots = false) :
            Lirion.solve(path; model = base_ctor(), algorithm = alg,
                optimizer = opt, outputs = false, plots = false)
        status = _classify(r.status)
        obj = isfinite(r.objective) ? r.objective : NaN
        iters = isa(r.iter_done, Integer) ? Int(r.iter_done) : 0
    catch err
        status = "error"
    end
    dt = time() - t0
    gap = (isfinite(obj) && isfinite(ref) && ref > 0) ?
          100.0 * (obj - ref) / ref : NaN
    return (status = status, obj = isfinite(obj) ? obj : 0.0,
            gap_pct = gap, iters = iters, time_s = dt)
end

function main()
    println("╔════════════════════════════════════════════════════════════╗")
    println("║  Lirion paper — case14 ALADIN v2 strategy sweep            ║")
    @printf("║  eps=%.0e max_iter=%d target gap<%.1f%%\n", EPS_ABS, MAX_ITER, GAP_TARGET)
    println("║  Output: ", OUT)
    println("╚════════════════════════════════════════════════════════════╝")
    flush(stdout)

    isdir(CASES_DIR) || mkpath(CASES_DIR)
    report = isfile(OUT) ? TOML.parsefile(OUT) : Dict{String, Any}()

    for variant in VARIANTS
        haskey(report, variant) || (report[variant] = Dict{String, Any}())
        for (base, base_ctor) in _base_specs()
            ref = _central_reference(variant, base, base_ctor)
            isfinite(ref) || continue
            bkey = string(base)
            haskey(report[variant], bkey) || (report[variant][bkey] = Dict{String, Any}())
            for (sect, coup) in ALGO_VARIANTS
                # SOCWR + DummyBus is unsupported in ALADIN (see
                # ALADINConfig guard in src/types.jl) — skip.
                base === :socwr && sect === :dummy_bus && continue
                key = "$(sect)_$(coup)"
                haskey(report[variant][bkey], key) &&
                    get(report[variant][bkey][key], "hit_target", false) &&
                    continue
                cell = Dict{String, Any}("strategies" => Dict{String, Any}())
                winner = nothing
                for (sname, skw) in STRATEGIES
                    print(@sprintf("  [%s/%s/%s/%s / %-18s] ",
                            variant, base, sect, coup, sname))
                    flush(stdout)
                    r = _probe(variant, base, base_ctor, sect, coup, skw, ref)
                    cell["strategies"][sname] = Dict(
                        "rho" => get(skw, :rho, NaN), "mu" => get(skw, :mu, NaN),
                        "rho_growth" => get(skw, :rho_growth, NaN),
                        "mu_growth" => get(skw, :mu_growth, NaN),
                        "adaptive" => get(skw, :adaptive, false),
                        "status" => r.status, "gap_pct" => r.gap_pct,
                        "iters" => r.iters, "time_s" => r.time_s,
                    )
                    @printf("%s iter=%d t=%.2fs gap=%s\n", r.status, r.iters, r.time_s,
                        isfinite(r.gap_pct) ? @sprintf("%+6.3f%%", r.gap_pct) : "NaN")
                    flush(stdout)
                    if r.status == "converged" && isfinite(r.gap_pct) && abs(r.gap_pct) < GAP_TARGET
                        winner = (name = sname, kw = skw, gap = r.gap_pct, iters = r.iters)
                        break
                    end
                end
                if winner !== nothing
                    cell["hit_target"] = true
                    cell["winner"] = Dict(
                        "strategy" => winner.name,
                        "rho" => get(winner.kw, :rho, NaN),
                        "mu" => get(winner.kw, :mu, NaN),
                        "rho_growth" => get(winner.kw, :rho_growth, NaN),
                        "mu_growth" => get(winner.kw, :mu_growth, NaN),
                        "adaptive" => get(winner.kw, :adaptive, false),
                        "gap_pct" => winner.gap, "iters" => winner.iters,
                    )
                else
                    cell["hit_target"] = false
                end
                report[variant][bkey][key] = cell
                _save_toml!(report, OUT)
            end
        end
    end

    report["_meta"] = Dict(
        "generated_at" => string(now()),
        "script" => "tune_case14_aladin_v2.jl",
        "eps_abs" => EPS_ABS, "eps_rel" => EPS_REL,
        "max_iter" => MAX_ITER, "gap_target" => GAP_TARGET,
    )
    _save_toml!(report, OUT)

    # Summary
    ok = 0; total = 0
    for v in VARIANTS, b in ("dc", "ac", "socwr")
        haskey(report, v) && haskey(report[v], b) || continue
        for (sect, coup) in ALGO_VARIANTS
            k = "$(sect)_$(coup)"
            haskey(report[v][b], k) || continue
            total += 1
            get(report[v][b][k], "hit_target", false) && (ok += 1)
        end
    end
    println("\nSummary: $ok/$total combos with gap<$(GAP_TARGET)%")
    println("Written to: $OUT")
end

main()
