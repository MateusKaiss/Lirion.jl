# ============================================================
# tune_case14.jl
# ============================================================
# Systematic ρ (and adaptive flag) sweep for ADMM and ALADIN on the
# case14_ieee benchmark. Produces two TOML reports that can be mined to
# populate `ADMM_OVERRIDES` and `ALADIN_OVERRIDES` in the paper_results
# scripts:
#
#   paper/cases/_tune_case14_admm.toml
#   paper/cases/_tune_case14_aladin.toml
#
# For every combination of (variant, base, sectioning, coupling) × (fixed,
# adaptive) we sweep ρ through a log-spaced ladder and record:
#   - first ρ that converges with gap < GAP_TARGET %
#   - best (smallest |gap|) among all sweeps
#   - per-ρ detail so paper plots can be built later.
#
# Stops early per combo when a ρ reaches gap < GAP_TARGET (pragmatic
# tuning, not full grid search). Set EXHAUSTIVE=true to force a full
# sweep of every ρ for every combo.
#
# Usage:
#     julia> include(joinpath(pkgdir(Lirion), "paper", "tune_case14.jl"))
#
# Environment variables (optional):
#   LIRION_TUNE_ALGO     = "admm" | "aladin" | "both"  (default: both)
#   LIRION_TUNE_VARIANT  = "std" | "api" | "sad" | "all" (default: all)
#   LIRION_TUNE_EXHAUSTIVE = "1" to disable early-stop
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

_HAVE_IPOPT || @warn "Ipopt not loaded — AC tuning runs will be skipped."

const TIME_LIMIT_S = 120.0

# Route solvers so Lirion's tuned silent factories kick in: passing
# `Ipopt.Optimizer` / `Clarabel.Optimizer` as a bare type triggers the
# `_silent_factory(::Type{...})` methods defined in
# `src/distributed_algorithms/optimizer.jl`, which carry the critical
# eps/tolerance/adaptive settings ALADIN's coordinator QP needs. A bare
# `OSQP.Optimizer` has NO silent factory registered, so passing `nothing`
# for DC lets Lirion fall through to its hand-tuned OSQP attributes
# (eps_abs=eps_rel=1e-6, adaptive_rho=true, qdldl) which are required to
# keep the DC ALADIN QP bounded.
_gurobi_opt() = optimizer_with_attributes(
    () -> Gurobi.Optimizer(Gurobi.Env()),
    "OutputFlag" => 0, "TimeLimit" => TIME_LIMIT_S,
    "FeasibilityTol" => 1e-7, "OptimalityTol" => 1e-7,
)

_solver_for_base(base::Symbol) = begin
    base === :dc    && return _HAVE_GUROBI ? _gurobi_opt() : nothing
    base === :ac    && return _HAVE_IPOPT  ? Ipopt.Optimizer : nothing
    base === :socwr && return _HAVE_GUROBI ? _gurobi_opt() : (_HAVE_CLARABEL ? Clarabel.Optimizer : nothing)
    error("unknown base $base")
end

# `nothing` from `_solver_for_base(:dc)` means "use Lirion's tuned OSQP
# default"; for AC/SOCWR it means "no solver available, skip this base".
_base_available(base::Symbol) = begin
    base === :dc    && return true  # Lirion's default OSQP is always routed
    base === :ac    && return _HAVE_IPOPT
    base === :socwr && return _HAVE_GUROBI || _HAVE_CLARABEL
    false
end

# ---- Paths ------------------------------------------------------------
const PAPER_DIR  = @__DIR__
const CASES_DIR  = joinpath(PAPER_DIR, "cases")
const CASE       = "case14_ieee"
const OUT_ADMM   = joinpath(CASES_DIR, "_tune_case14_admm.toml")
const OUT_ALADIN = joinpath(CASES_DIR, "_tune_case14_aladin.toml")

# ---- Tuning knobs -----------------------------------------------------
const GAP_TARGET  = 1.0                       # percent, success criterion
const EXHAUSTIVE  = get(ENV, "LIRION_TUNE_EXHAUSTIVE", "0") == "1"
const MAX_ITER    = 5_000
const EPS_ABS     = 1e-4
const EPS_REL     = 1e-4

# ρ ladder: spans 1e1 … 1e8 so the user sees the entire usable range.
# Ordered per coupling so the most likely winner is tried first.
const RHO_ORDERS = Dict(
    :state  => [1e5, 1e4, 1e6, 1e3, 1e7, 1e2, 1e8, 1e1],
    :power  => [1e3, 1e4, 1e2, 1e5, 1e1, 1e6, 1e7, 1e8],
    :hybrid => [1e3, 1e4, 1e2, 1e5, 1e1, 1e6, 1e7, 1e8],
)

# ALADIN ρ ladder is more conservative: adaptive growth pushes ρ up
# during the run, so starting lower is usually sufficient.
const ALADIN_RHO_ORDERS = Dict(
    :state  => [1e3, 1e4, 1e2, 1e5, 1e1, 1e6],
    :power  => [1e3, 1e4, 1e2, 1e5, 1e1, 1e6],
    :hybrid => [1e3, 1e4, 1e2, 1e5, 1e1, 1e6],
)

const ADAPTIVE_FLAGS = [false, true]

const VARIANTS = let v = get(ENV, "LIRION_TUNE_VARIANT", "all")
    v == "all" ? ["std", "api", "sad"] : [v]
end

const ALGO_SELECT = let a = get(ENV, "LIRION_TUNE_ALGO", "both")
    a == "both" ? [:admm, :aladin] : [Symbol(a)]
end

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

# ---- Cache centralised reference so we compute it once per (variant, base).
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

# ---- Single ADMM probe -----------------------------------------------
function _probe_admm(variant::String, base::Symbol, base_ctor,
                     sectioning::Symbol, coupling::Symbol,
                     rho::Float64, adaptive::Bool, ref_obj::Float64)
    path = _partition_file(variant, sectioning)
    t0 = time()
    status = "error"; obj = NaN; iters = 0
    try
        alg = admm(;
            coupling = coupling, sectioning = sectioning,
            max_iter = MAX_ITER, rho = rho,
            eps_abs = EPS_ABS, eps_rel = EPS_REL,
            adaptive_rho = adaptive,
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
    return (status = status, objective = isfinite(obj) ? obj : 0.0,
            gap_pct = gap, iters = iters, time_s = dt)
end

# ---- Single ALADIN probe ---------------------------------------------
function _probe_aladin(variant::String, base::Symbol, base_ctor,
                       sectioning::Symbol, coupling::Symbol,
                       rho::Float64, adaptive::Bool, ref_obj::Float64)
    path = _partition_file(variant, sectioning)
    t0 = time()
    status = "error"; obj = NaN; iters = 0
    # Fix mu = rho/10 to keep the ALADIN coordinator QP strictly
    # convex across the ρ sweep (matches the test suite's working
    # rho=1e4/mu=1e3 ratio on case14/DC/hybrid).
    mu_val = rho / 10
    try
        alg = aladin(;
            coupling = coupling, sectioning = sectioning,
            max_iter = 300, rho = rho, mu = mu_val,
            eps_abs = EPS_ABS, eps_rel = EPS_REL,
            adaptive = adaptive,
            rho_growth = 1.3, mu_growth = 2.0,
            rho_upper = 1e8, mu_upper = 1e8,
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
    return (status = status, objective = isfinite(obj) ? obj : 0.0,
            gap_pct = gap, iters = iters, time_s = dt)
end

# ---- Combo sweep for a single algorithm ------------------------------
function _sweep_combo!(report::Dict, out_path::AbstractString,
                       algo::Symbol,
                       variant::String,
                       base::Symbol, base_ctor,
                       sectioning::Symbol, coupling::Symbol,
                       ref_obj::Float64)
    haskey(report, variant) || (report[variant] = Dict{String, Any}())
    bkey = string(base)
    haskey(report[variant], bkey) || (report[variant][bkey] = Dict{String, Any}())
    combo_key = "$(sectioning)_$(coupling)"
    haskey(report[variant][bkey], combo_key) ||
        (report[variant][bkey][combo_key] = Dict{String, Any}())
    entry = report[variant][bkey][combo_key]

    rho_list = algo === :admm ? RHO_ORDERS[coupling] : ALADIN_RHO_ORDERS[coupling]
    best_fixed    = nothing
    best_adaptive = nothing

    for adaptive in ADAPTIVE_FLAGS
        section_key = adaptive ? "adaptive" : "fixed"
        haskey(entry, section_key) || (entry[section_key] = Dict{String, Any}())
        sweep_entry = entry[section_key]
        haskey(sweep_entry, "sweep") || (sweep_entry["sweep"] = Dict{String, Any}())

        found = false
        best  = nothing
        for rho in rho_list
            rkey = @sprintf("rho=%.0e", rho)
            if haskey(sweep_entry["sweep"], rkey) &&
               haskey(sweep_entry["sweep"][rkey], "status")
                cached = sweep_entry["sweep"][rkey]
                print(@sprintf("  [%s/%s/%s/%s/%-6s/%s/ρ=%7.0e] cached  ",
                        algo, variant, base, sectioning, coupling, section_key, rho))
                if isfinite(get(cached, "gap_pct", NaN))
                    print(@sprintf("gap=%+6.3f%% iters=%d t=%.2fs %s",
                            cached["gap_pct"], get(cached, "iters", 0),
                            get(cached, "time_s", 0.0), get(cached, "status", "?")))
                end
                println()
                flush(stdout)
                res = (
                    status = String(get(cached, "status", "error")),
                    objective = Float64(get(cached, "objective", 0.0)),
                    gap_pct = Float64(get(cached, "gap_pct", NaN)),
                    iters = Int(get(cached, "iters", 0)),
                    time_s = Float64(get(cached, "time_s", 0.0)),
                )
            else
                print(@sprintf("  [%s/%s/%s/%s/%-6s/%s/ρ=%7.0e] ",
                        algo, variant, base, sectioning, coupling, section_key, rho))
                flush(stdout)
                res = algo === :admm ?
                      _probe_admm(variant, base, base_ctor, sectioning, coupling,
                                  rho, adaptive, ref_obj) :
                      _probe_aladin(variant, base, base_ctor, sectioning, coupling,
                                    rho, adaptive, ref_obj)
                sweep_entry["sweep"][rkey] = Dict(
                    "rho"       => rho,
                    "status"    => res.status,
                    "objective" => res.objective,
                    "gap_pct"   => res.gap_pct,
                    "iters"     => res.iters,
                    "time_s"    => res.time_s,
                )
                _save_toml!(report, out_path)
                println(@sprintf("%s iters=%d t=%.2fs gap=%s",
                        res.status, res.iters, res.time_s,
                        isfinite(res.gap_pct) ? @sprintf("%+6.3f%%", res.gap_pct) : "NaN"))
                flush(stdout)
            end

            # Track best |gap| with converged status.
            if res.status == "converged" && isfinite(res.gap_pct)
                candidate = (rho = rho, gap = res.gap_pct,
                             iters = res.iters, time_s = res.time_s)
                if best === nothing || abs(candidate.gap) < abs(best.gap)
                    best = candidate
                end
                if abs(res.gap_pct) < GAP_TARGET
                    found = true
                    if !EXHAUSTIVE
                        break
                    end
                end
            end
        end

        if best !== nothing
            sweep_entry["best"] = Dict(
                "rho"     => best.rho,
                "gap_pct" => best.gap,
                "iters"   => best.iters,
                "time_s"  => best.time_s,
                "hit_target" => found,
            )
        else
            sweep_entry["best"] = Dict(
                "rho" => NaN, "gap_pct" => NaN, "iters" => 0, "time_s" => 0.0,
                "hit_target" => false,
            )
        end

        if adaptive
            best_adaptive = best
        else
            best_fixed = best
        end
        _save_toml!(report, out_path)
    end

    return (fixed = best_fixed, adaptive = best_adaptive)
end

function _print_header(algo::Symbol, out_path::AbstractString)
    println()
    println("╔════════════════════════════════════════════════════════════╗")
    @printf("║  Lirion paper — case14 tuning sweep (%s)\n", uppercase(string(algo)))
    println("║  Case:      ", CASE)
    println("║  Variants:  ", join(VARIANTS, ", "))
    @printf("║  Solvers:   Gurobi=%s  Ipopt=%s  Clarabel=%s  OSQP=%s\n",
            _HAVE_GUROBI, _HAVE_IPOPT, _HAVE_CLARABEL, _HAVE_OSQP)
    @printf("║  eps_abs/rel = %.0e, max_iter = %d, gap target = %.1f%%\n",
            EPS_ABS, MAX_ITER, GAP_TARGET)
    @printf("║  Exhaustive: %s\n", EXHAUSTIVE)
    println("║  Started:   ", now())
    println("║  Output:    ", out_path)
    println("╚════════════════════════════════════════════════════════════╝")
    flush(stdout)
end

function run_algo_sweep(algo::Symbol, out_path::AbstractString)
    _print_header(algo, out_path)
    isdir(CASES_DIR) || mkpath(CASES_DIR)
    report = isfile(out_path) ? TOML.parsefile(out_path) : Dict{String, Any}()

    winners = Tuple[]
    for variant in VARIANTS
        for (base, base_ctor) in _base_specs()
            ref_obj = _central_reference(variant, base, base_ctor)
            if !isfinite(ref_obj)
                @printf("[%s/%-5s] SKIP (centralised reference unavailable)\n", variant, base)
                flush(stdout)
                continue
            end
            @printf("\n=== %s / %s / %s — ref obj = %.2f ===\n",
                    CASE, variant, base, ref_obj)
            flush(stdout)
            for (sectioning, coupling) in ALGO_VARIANTS
                # SOCWR + DummyBus is unsupported in ALADIN (guard in
                # ALADINConfig). Skip cleanly so the sweep does not
                # invoke it and trip the error.
                if algo === :aladin && base === :socwr && sectioning === :dummy_bus
                    continue
                end
                res = _sweep_combo!(report, out_path, algo, variant, base, base_ctor,
                                    sectioning, coupling, ref_obj)
                push!(winners, (variant, base, sectioning, coupling, res))
            end
        end
    end

    report["_meta"] = Dict(
        "generated_at" => string(now()),
        "script"       => "tune_case14.jl",
        "algo"         => string(algo),
        "case"         => CASE,
        "eps_abs"      => EPS_ABS,
        "eps_rel"      => EPS_REL,
        "max_iter"     => MAX_ITER,
        "gap_target"   => GAP_TARGET,
        "exhaustive"   => EXHAUSTIVE,
    )
    _save_toml!(report, out_path)

    # Final summary
    println("\n" * "="^92)
    @printf("SUMMARY — %s — case14 (fixed ρ / adaptive ρ)\n", uppercase(string(algo)))
    println("="^92)
    @printf("%-3s %-5s %-9s %-6s  %14s   %14s\n",
            "var", "base", "sect", "coup",
            "fixed: ρ gap it", "adapt: ρ gap it")
    for (variant, base, sectioning, coupling, res) in winners
        f = res.fixed    === nothing ? "  —    —     —" :
            @sprintf("%.0e %+6.2f%% %4d", res.fixed.rho, res.fixed.gap, res.fixed.iters)
        a = res.adaptive === nothing ? "  —    —     —" :
            @sprintf("%.0e %+6.2f%% %4d", res.adaptive.rho, res.adaptive.gap, res.adaptive.iters)
        @printf("%-3s %-5s %-9s %-6s  %14s   %14s\n",
                variant, base, sectioning, coupling, f, a)
    end
    flush(stdout)
    println("\nWritten to: ", out_path)
    return report
end

function main()
    for algo in ALGO_SELECT
        out_path = algo === :admm ? OUT_ADMM : OUT_ALADIN
        run_algo_sweep(algo, out_path)
    end
end

main()
