# ============================================================
# paper/scripts/experiments/lib/tuning_common.jl
# ============================================================
# Declared-grid tuning engine (D0 smoke + later the full tuning).
# - reads grids/criteria from TOML;
# - runs every parameter attempt with ROBUST exception capture
#   (TaskFailedException / CompositeException / ITERATION_LIMIT /
#    NO_SOLUTION / any solve error) -> a failed attempt row, never a
#   driver crash; the search continues;
# - selects the best attempt per cell by a declared lexicographic rule;
# - emits rows in the shared distributed schema (+ tuning columns).
# Uses the official distributed solver config (C0): DC=OSQP, SOCWR=
# Clarabel(hardened), AC=Ipopt+MA57 via LIRION_HSLLIB.
# ============================================================

include(joinpath(@__DIR__, "distributed_common.jl"))

import MathOptInterface as MOI
using JuMP: optimizer_with_attributes

const CONFIG_DIR = joinpath(PAPER_DIR, "configs")

# ---- Robust exception unwrapping -------------------------------------
"""Recursively unwrap CompositeException / TaskFailedException to the
innermost exception (the real solver/subproblem error)."""
function unwrap_root(e)
    e isa CompositeException && !isempty(e.exceptions) && return unwrap_root(e.exceptions[1])
    if e isa TaskFailedException
        under = try; e.task.result; catch; e; end
        under === e && return e
        return unwrap_root(under)
    end
    return e
end

function failure_reason_of(e)
    root = unwrap_root(e)
    outer = string(nameof(typeof(e)))
    msg = first(split(sprint(showerror, root), "\n"))
    # collapse the multi-line _require_optimal! message to a compact line
    msg = replace(strip(msg), r"\s+" => " ")
    return e === root ? msg : string(outer, " -> ", msg)
end

function classify_attempt(e)
    root = unwrap_root(e)
    low = lowercase(sprint(showerror, root))
    if root isa MethodError || occursin("not implemented", low) ||
       occursin("not supported", low) || occursin("unsupported", low)
        return "not_implemented"
    elseif root isa ArgumentError &&
           (occursin("proximal_mode", low) || occursin("pflow", low) ||
            occursin("coupling=", low))
        return "skipped_by_design"
    else
        return "unstable_or_fails"   # incl. ITERATION_LIMIT / NO_SOLUTION / subproblem
    end
end

# ---- Official distributed solver per formulation (C0) ----------------
function ac_ma57_optimizer()
    have_hsl() || return nothing        # fall back to Lirion default (MUMPS)
    base = Lirion._SILENT_FACTORY_AC
    keep = Pair[]
    for p in base.params
        nm = p.first isa MOI.RawOptimizerAttribute ? p.first.name : ""
        (nm == "linear_solver" || nm == "hsllib") && continue
        push!(keep, p)
    end
    push!(keep, MOI.RawOptimizerAttribute("linear_solver") => "ma57")
    push!(keep, MOI.RawOptimizerAttribute("hsllib") => hsllib_path())
    return optimizer_with_attributes(base.optimizer_constructor, keep...)
end

# DC -> OSQP (nothing = Lirion default); SOCWR -> Clarabel (nothing);
# AC  -> Ipopt + MA57 (custom), or nothing (MUMPS) if LIRION_HSLLIB unset.
official_dist_optimizer(f) = f == "ac" ? ac_ma57_optimizer() : nothing
official_solver_label(f) =
    f == "ac" ? (have_hsl() ? "Ipopt+MA57 (dist, official)" : "Ipopt+MUMPS (dist, fallback)") :
    dist_solver_label(f)

# ---- Grid / criteria loading -----------------------------------------
load_grid(name) = TOML.parsefile(joinpath(CONFIG_DIR, name))
load_criteria(name) = TOML.parsefile(joinpath(CONFIG_DIR, name))

"""Expand a grid TOML into a vector of parameter NamedTuples per algorithm."""
function param_sets(algo, grid)
    g = grid["global"]; G = grid["grid"]
    if algo == "admm"
        return [(; rho = rho, mu = NaN, rho_growth = NaN,
                 max_iter = Int(g["max_iter"]), eps_abs = g["eps_abs"],
                 eps_rel = g["eps_rel"], adaptive_rho = Bool(g["adaptive_rho"]))
                for rho in Float64.(G["rho"])]
    else # aladin
        return [(; rho = rho, mu = mu, rho_growth = Float64(g["rho_growth"]),
                 max_iter = Int(g["max_iter"]), eps_abs = NaN, eps_rel = NaN,
                 adaptive_rho = false)
                for rho in Float64.(G["rho"]) for mu in Float64.(G["mu"])]
    end
end

param_label(algo, p) = algo == "admm" ?
    @sprintf("rho=%.0e", p.rho) :
    @sprintf("rho=%.0e,mu=%.0e,growth=%.2f", p.rho, p.mu, p.rho_growth)

function build_tuned_algorithm(algo, sect, coup, p)
    if algo == "admm"
        return Lirion.admm(; coupling = Symbol(coup), sectioning = Symbol(sect),
            max_iter = p.max_iter, rho = p.rho,
            eps_abs = p.eps_abs, eps_rel = p.eps_rel, adaptive_rho = p.adaptive_rho)
    else
        return Lirion.aladin(; coupling = Symbol(coup), sectioning = Symbol(sect),
            max_iter = p.max_iter, rho = p.rho, mu = p.mu, rho_growth = p.rho_growth)
    end
end

# ---- One attempt (robust) --------------------------------------------
const TUNING_COLUMNS = [
    :case, :variant, :k, :formulation, :algorithm, :sectioning, :coupling,
    :rho, :mu, :rho_growth, :parameter_set, :attempt_id, :status, :objective,
    :centralized_objective, :gap_pct, :iterations, :primal_residual,
    :dual_residual, :runtime_s, :solver, :convergence_flag, :support_status,
    :failure_reason, :is_selected, :timestamp,
]

function run_attempt(cell, p, attempt_id; t1)
    (case, k, f, algo, sect, coup) = cell
    path = partitioned_case_path(case, "std", k)
    cen = get(t1, (case, f), NaN)
    opt = official_dist_optimizer(f)

    status = ""; obj = NaN; iters = 0; pr = NaN; dr = NaN
    rt = NaN; conv = false; support = ""; fail = ""
    t0 = time()
    try
        alg = build_tuned_algorithm(algo, sect, coup, p)
        r = Lirion.solve(path; model = model_ctor(f)(), algorithm = alg,
                         optimizer = opt, outputs = false, plots = false)
        rt = time() - t0
        obj = Float64(r.objective); status = string(r.status)
        iters = r.iter_done === nothing ? 0 : r.iter_done
        pr, dr = extract_residuals(r)
        conv = (r.status == :converged)
        support = isfinite(obj) ? "supported" : "unstable_or_fails"
        isfinite(obj) || (fail = "non-finite objective")
    catch e
        rt = time() - t0
        support = classify_attempt(e); fail = failure_reason_of(e); status = "error"
    end
    gap = (isnan(obj) || isnan(cen) || cen == 0) ? NaN : 100 * (obj - cen) / abs(cen)

    return Dict{Symbol,Any}(
        :case => case, :variant => "std", :k => k, :formulation => f,
        :algorithm => algo, :sectioning => sect, :coupling => coup,
        :rho => p.rho, :mu => p.mu, :rho_growth => p.rho_growth,
        :parameter_set => param_label(algo, p), :attempt_id => attempt_id,
        :status => status, :objective => obj, :centralized_objective => cen,
        :gap_pct => gap, :iterations => iters, :primal_residual => pr,
        :dual_residual => dr, :runtime_s => rt, :solver => official_solver_label(f),
        :convergence_flag => conv, :support_status => support,
        :failure_reason => fail, :is_selected => false,
        :timestamp => string(Dates.now()),
    )
end

# ---- Lexicographic selection -----------------------------------------
function select_best!(rows, gap_target)
    isempty(rows) && return rows
    _res(r) = isnan(r[:primal_residual]) ? Inf : r[:primal_residual]
    _gapok(r) = (!isnan(r[:gap_pct]) && abs(r[:gap_pct]) < gap_target) ? 0 : 1
    key(r) = (r[:convergence_flag] ? 0 : 1, _gapok(r), _res(r),
              r[:iterations] == 0 ? typemax(Int) : r[:iterations], r[:runtime_s])
    best = argmin(map(key, rows))
    rows[best][:is_selected] = true
    return rows
end

# NamedTuple view for the shared CSV/MD writers (which use getproperty).
as_nt(d::Dict{Symbol,Any}) = (; (c => d[c] for c in TUNING_COLUMNS)...)
