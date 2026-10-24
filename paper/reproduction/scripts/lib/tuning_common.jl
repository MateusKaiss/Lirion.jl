# Tuning helpers: official solver configs per formulation, robust attempt
# runner, grid loading, centralized oracle solvers.

include(joinpath(@__DIR__, "distributed_common.jl"))

import MathOptInterface as MOI
using JuMP: optimizer_with_attributes
# central_opt builds Ipopt/Clarabel optimizers, so import them here
using Ipopt, Clarabel

const CONFIG_DIR = joinpath(PAPER_DIR, "configs")

# ---- git guard: clean tree at the expected commit ----------------------
const _GIT_GUARD_PATHS = ("src", "Project.toml", "Manifest.toml",
                          "paper/scripts/experiments", "paper/cases", "paper/configs")
function git_guard!()
    root = try
        readchomp(`git -C $(PAPER_DIR) rev-parse --show-toplevel`)
    catch e
        error("D4 git guard: $(PAPER_DIR) is not inside a git checkout: $(sprint(showerror, e))")
    end
    # the loaded Lirion must be the one this checkout holds
    loaded = normpath(pkgdir(Lirion))
    normpath(root) == loaded || error(
        "pkgdir guard: the loaded Lirion package lives in\n    $(loaded)\n" *
        "but the guarded checkout is\n    $(normpath(root))\n" *
        "(wrong --project / environment). Activate the project's " *
        "environment so pkgdir(Lirion) == the guarded checkout.")
    sha_full = readchomp(`git -C $(root) rev-parse HEAD`)
    st = readchomp(`git -C $(root) status --porcelain --untracked-files=no -- $(_GIT_GUARD_PATHS)`)
    isempty(st) || error(
        "D4 git guard: the working tree is DIRTY in paths that affect the run " *
        "(src/, Project.toml, Manifest.toml or the driver scripts):\n" * st *
        "\nCommit the changes (the reproduction must run committed, frozen code) " *
        "and set LIRION_EXPECT_SHA to that commit.")
    expect = get(ENV, "LIRION_EXPECT_SHA", "")
    isempty(expect) && error(
        "D4 git guard: LIRION_EXPECT_SHA is not set. Set it to the frozen reproduction " *
        "commit (any unique prefix of the full SHA); HEAD is $(sha_full[1:12]).")
    startswith(sha_full, lowercase(strip(expect))) || error(
        "D4 git guard: checkout SHA $(sha_full[1:12]) ≠ expected frozen commit " *
        "LIRION_EXPECT_SHA=$(expect). Refusing to run the wrong code.")
    return sha_full[1:12]
end

# git guard + no trajectory-altering opt-ins in the env
function reproduction_guard!()
    sha = git_guard!()
    haskey(ENV, "LIRION_GAUGEFIX") && error(
        "reproduction guard: LIRION_GAUGEFIX=$(ENV["LIRION_GAUGEFIX"]) is set. It changes " *
        "the ALADIN trajectory and is not part of the reproduction configuration — " *
        "unset it (the runtime manifest records LIRION_GAUGEFIX=unset for official runs).")
    haskey(ENV, "LIRION_ZONE_WARMSTART") && error(
        "reproduction guard: LIRION_ZONE_WARMSTART=$(ENV["LIRION_ZONE_WARMSTART"]) is set. It " *
        "enables the opt-in ADMM zone warm-start and changes the ADMM trajectory; it is " *
        "not part of the reproduction configuration (sealed runs are cold-start) — unset it.")
    return sha
end

# pad a stall-cut trajectory with its last value up to max_iter (no-op otherwise)
function pad_frozen_tail(v::Vector{Float64}, run_manifest, max_iter::Int)
    ex = Int(get(run_manifest, "EARLY_STALL_EXIT", 0))
    (ex > 0 && !isempty(v) && length(v) < max_iter) || return v, false
    return vcat(v, fill(v[end], max_iter - length(v))), true
end

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

# ---- Official distributed solver per formulation ---------------------
function ac_ma57_optimizer()
    if !have_hsl()
        _strict_solvers() && error(
            "STRICT SOLVERS: MA57 unavailable (LIRION_HSLLIB unset/invalid) — " *
            "official AC runs must not silently fall back to MUMPS.")
        return nothing                  # fall back to Lirion default (MUMPS)
    end
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

# DC OSQP fallback: eps=1e-5 + polish (the 1e-6 default hits max_iter on large zones)
dc_osqp_optimizer() = optimizer_with_attributes(Lirion.OSQP.Optimizer,
    "eps_abs" => 1e-5, "eps_rel" => 1e-5, "max_iter" => 100_000, "verbose" => false,
    "adaptive_rho" => true, "linsys_solver" => "qdldl", "polish" => true, "scaling" => 15)

# DC zones: Gurobi for both ADMM and ALADIN (same zone solver keeps the timing
# comparison fair); one shared Env per process; OSQP fallback without a license
const _GRB_MOD = try; @eval import Gurobi; Gurobi; catch; nothing; end
const _GRB_ENV = Ref{Any}(nothing)

# LIRION_STRICT_SOLVERS=1: a solver fallback is a hard error
_strict_solvers() = get(ENV, "LIRION_STRICT_SOLVERS", "0") == "1"

function dc_gurobi_optimizer()
    if _GRB_MOD === nothing
        _strict_solvers() && error(
            "STRICT SOLVERS: Gurobi unavailable (import failed) — official DC " *
            "runs must not silently fall back to OSQP. Fix the Gurobi install " *
            "or unset LIRION_STRICT_SOLVERS for a declared non-official run.")
        return dc_osqp_optimizer()
    end
    try
        _GRB_ENV[] === nothing && (_GRB_ENV[] = _GRB_MOD.Env())
        # TimeLimit: a stuck zone solve returns TIME_LIMIT instead of hanging the run
        return optimizer_with_attributes(() -> _GRB_MOD.Optimizer(_GRB_ENV[]),
                                         "OutputFlag" => 0, "TimeLimit" => 60.0, "Threads" => 2)
    catch e
        _strict_solvers() && error(
            "STRICT SOLVERS: Gurobi Env failed ($(first(split(sprint(showerror, e), '\n')))) — " *
            "official DC runs must not silently fall back to OSQP.")
        return dc_osqp_optimizer()
    end
end

# DC coordinator QP on OSQP (unused now; Clarabel below is the universal one)
aladin_coord_optimizer() = optimizer_with_attributes(Lirion.OSQP.Optimizer,
    "eps_abs" => 1e-4, "eps_rel" => 1e-4, "max_iter" => 100_000, "verbose" => false,
    "adaptive_rho" => true, "polish" => true)

# coordinator QP: Clarabel (interior-point) solves the ill-conditioned large
# QPs where OSQP/Gurobi fail
aladin_coord_ac_optimizer() = optimizer_with_attributes(Lirion.Clarabel.Optimizer,
    "verbose" => false, "tol_feas" => 1e-7, "tol_gap_abs" => 1e-7,
    "tol_gap_rel" => 1e-7, "max_iter" => 2000,
    # time_limit: a degenerate coordinator QP must hand control back to julia
    "time_limit" => 60.0)

# same Clarabel coordinator for every formulation
aladin_coord_for(f, case = "") = aladin_coord_ac_optimizer()

# LIRION_ZONE_TOL_TIGHT=1: tightened zone solvers
_zone_tol_tight() = get(ENV, "LIRION_ZONE_TOL_TIGHT", "0") == "1"

# LIRION_SOCWR_ZONE_TOL: tolerance of the tight SOCWR zone solver (default 1e-9)
socwr_tight_optimizer() = optimizer_with_attributes(Clarabel.Optimizer,
    "verbose" => false,
    "tol_feas" => parse(Float64, get(ENV, "LIRION_SOCWR_ZONE_TOL", "1e-9")),
    "tol_gap_abs" => parse(Float64, get(ENV, "LIRION_SOCWR_ZONE_TOL", "1e-9")),
    "tol_gap_rel" => parse(Float64, get(ENV, "LIRION_SOCWR_ZONE_TOL", "1e-9")),
    "max_iter" => 50000,
    "direct_solve_method" => :qdldl, "equilibrate_enable" => true,
    "time_limit" => 120.0)

# tight AC zones: Ipopt+MA57 at 1e-9
function ac_ma57_tight_optimizer()
    base = ac_ma57_optimizer()
    # Pair{...,Any}: um array literal promoveria 0 (Int, exigido pelo Ipopt
    # para acceptable_iter) a 0.0 junto com os 1e-9 — "valid option, wrong type"
    extra = Pair{MOI.RawOptimizerAttribute, Any}[
        MOI.RawOptimizerAttribute("tol") => 1e-9,
        MOI.RawOptimizerAttribute("dual_inf_tol") => 1e-9,
        MOI.RawOptimizerAttribute("constr_viol_tol") => 1e-9,
        MOI.RawOptimizerAttribute("compl_inf_tol") => 1e-9,
        # 10 (not 0): a zone that cannot hit 1e-9 in the early transient falls
        # back on an acceptable point instead of crashing the run
        MOI.RawOptimizerAttribute("acceptable_iter") => 10,
    ]
    return optimizer_with_attributes(base.optimizer_constructor,
        vcat(base.params, extra)...)
end

function dc_gurobi_tight_optimizer()
    base = dc_gurobi_optimizer()   # STRICT path already errors if Gurobi absent
    return optimizer_with_attributes(base.optimizer_constructor,
        vcat(base.params, [JuMP.MOI.RawOptimizerAttribute("OptimalityTol") => 1e-9,
                           JuMP.MOI.RawOptimizerAttribute("FeasibilityTol") => 1e-9])...)
end

# DC -> Gurobi (both algos); SOCWR -> Clarabel default (nothing); AC -> Ipopt+MA57
function official_dist_optimizer(f, algo = "aladin")
    f == "ac" && return _zone_tol_tight() ? ac_ma57_tight_optimizer() :
                        ac_ma57_optimizer()
    f == "dc" && return _zone_tol_tight() ? dc_gurobi_tight_optimizer() :
                        dc_gurobi_optimizer()
    f == "socwr" && _zone_tol_tight() && return socwr_tight_optimizer()
    # socwrnl: nothing resolves to _socwrnl_tight_factory (Ipopt+MA57, 1e-9);
    # a concrete factory here would error under STRICT
    f == "socwrnl" && return nothing
    return nothing   # socwr -> Lirion default (Clarabel)
end
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

# ---- Centralized oracle solvers (single source of truth for every gap) --

const _CENTRAL_GRB = try; @eval import Gurobi; Gurobi; catch; nothing; end
_central_grbenv = nothing
function central_gurobi_opt()
    if _CENTRAL_GRB === nothing
        _strict_solvers() && error(
            "STRICT SOLVERS: Gurobi unavailable — the official DC oracle must " *
            "not silently fall back to Ipopt.")
        return optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)
    end
    try
        global _central_grbenv
        _central_grbenv === nothing && (_central_grbenv = _CENTRAL_GRB.Env())
        return optimizer_with_attributes(() -> _CENTRAL_GRB.Optimizer(_central_grbenv), "OutputFlag"=>0)
    catch e
        _strict_solvers() && error(
            "STRICT SOLVERS: Gurobi Env failed — official DC oracle must not " *
            "silently fall back to Ipopt.")
        return optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)
    end
end

# the oracle must be tighter than the gaps it referees: 1e-8
central_clarabel_hard() = optimizer_with_attributes(Clarabel.Optimizer, "verbose"=>false,
    "tol_feas"=>1e-8, "tol_gap_abs"=>1e-8, "tol_gap_rel"=>1e-8, "max_iter"=>20000,
    "direct_solve_method"=>:qdldl, "equilibrate_enable"=>true)

central_ipopt_ac() = have_hsl() ?
    optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "linear_solver"=>"ma57", "hsllib"=>hsllib_path()) :
    optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)

"""
    central_opt(f)

Optimizer for the CENTRALIZED oracle of formulation `f` ("dc" | "socwr" | "ac").
Single source of truth -- do not re-define this per script.
"""
# socwrnl: nothing -> bare Ipopt, same as the oracle the tuning drivers used
central_opt(f) = f == "dc"      ? central_gurobi_opt() :
                 f == "socwr"   ? central_clarabel_hard() :
                 f == "socwrnl" ? nothing :
                                  central_ipopt_ac()
