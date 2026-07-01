# ============================================================
# paper/scripts/experiments/run_centralized.jl
# ============================================================
# Stage A — Lirion centralized baseline (T1). One row per
# (case × variant × formulation); centralized is independent of k, so it
# is NOT repeated per number of zones. This baseline is later the oracle
# for the distributed (ADMM/ALADIN) gap.
#
# Modes (env LIRION_CENTRALIZED_MODE):
#   "smoke" (default) — case57_ieee/std + case118_ieee/std, DC/SOCWR/AC
#   "full"            — all distinct clean-set cases (std), DC/SOCWR/AC
#
# Outputs:
#   results/centralized/centralized_baseline_raw.csv   (scratch, git-ignored)
#   artifacts/tables/T1_centralized_baseline.csv       (curated, versioned)
#   artifacts/tables/T1_centralized_baseline.md        (curated, versioned)
#
# No ADMM/ALADIN/tuning. PowerModels validation (T0) is a separate stage.
# Reproducible via: include(".../experiments/run_centralized.jl")
# ============================================================

include(joinpath(@__DIR__, "lib", "experiment_common.jl"))

using Ipopt, Clarabel
using LinearAlgebra
using JuMP: optimizer_with_attributes

# Gurobi is the official DC solver. It loads only when the local install
# is reachable (GUROBI_HOME + GUROBI_JL_USE_GUROBI_JLL=false); if not, the
# driver warns and the DC result is flagged non-official (no silent
# downgrade). No Gurobi binary/license is stored in the repo.
const _GUROBI = try
    @eval import Gurobi
    Gurobi
catch
    nothing
end
_grb_env = nothing

# --- Official reproducible performance config -------------------------
# Deterministic timing for the paper: single-threaded BLAS (more threads
# did not help these sparse OPFs and adds variance). Julia threads are
# set via JULIA_NUM_THREADS=1 in the environment (Ipopt is not
# Julia-threaded).
BLAS.set_num_threads(1)

# --- Solver factories -------------------------------------------------
# AC: Ipopt + MA57 when LIRION_HSLLIB is set (official); otherwise the
# driver refuses to silently downgrade in "full" mode (see assert below).
function _ipopt_ac()
    if have_hsl()
        return optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0,
            "linear_solver" => "ma57", "hsllib" => hsllib_path())
    else
        return optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
    end
end

_ipopt_default() = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)

# SOCWR: hardened Clarabel (mirrors the test suite's DEFAULT_SOCWR_SOLVER)
# so the baseline oracle does not sit at slow_progress / almost_optimal.
_clarabel_hardened() = optimizer_with_attributes(
    Clarabel.Optimizer,
    "verbose" => false,
    "tol_feas" => 1e-4, "tol_gap_abs" => 1e-4, "tol_gap_rel" => 1e-4,
    "max_iter" => 1200,
    "dynamic_regularization_enable" => true,
    "dynamic_regularization_eps" => 1e-10,
    "dynamic_regularization_delta" => 1e-4,
    "static_regularization_enable" => true,
    "equilibrate_enable" => true,
    "equilibrate_min_scaling" => 1e-4, "equilibrate_max_scaling" => 1e4,
    "direct_solve_method" => :qdldl,
    "iterative_refinement_reltol" => 1e-9, "iterative_refinement_abstol" => 1e-9,
    "max_step_fraction" => 0.7, "linesearch_backtrack_step" => 0.5,
    "tol_ktratio" => 1e-5, "min_terminate_step_length" => 1e-6,
    "min_switch_step_length" => 5e-2,
    "reduced_tol_gap_abs" => 5e-4, "reduced_tol_gap_rel" => 5e-4,
    "reduced_tol_feas" => 5e-4, "reduced_tol_infeas_abs" => 5e-11,
    "reduced_tol_infeas_rel" => 5e-4, "reduced_tol_ktratio" => 1e-4,
)

# DC: Gurobi (official). Returns nothing if the Gurobi install/license is
# not reachable — never throws — so the caller can fall back + flag it.
function _gurobi_dc()
    _GUROBI === nothing && return nothing
    try
        global _grb_env
        _grb_env === nothing && (_grb_env = _GUROBI.Env())
        return optimizer_with_attributes(() -> _GUROBI.Optimizer(_grb_env),
            "OutputFlag" => 0)
    catch
        return nothing
    end
end

have_gurobi() = _gurobi_dc() !== nothing

_optimizer(f) =
    f == "socwr" ? _clarabel_hardened() :
    f == "ac"    ? _ipopt_ac() :
    something(_gurobi_dc(), _ipopt_default())   # dc: Gurobi official, Ipopt = non-official fallback

# Case lists (variant = std for the baseline; api/sad added later).
const SMOKE_CASES = [("case57_ieee", "std"), ("case118_ieee", "std")]
const FULL_CASES = [
    ("case57_ieee", "std"), ("case73_ieee_rts", "std"), ("case118_ieee", "std"),
    ("case179_goc", "std"), ("case300_ieee", "std"), ("case500_goc", "std"),
    ("case588_sdet", "std"), ("case1354_pegase", "std"),
]

# One discarded warm-up solve per main formulation, so the first
# JIT-contaminated solve never lands in the reported timings. Uses the
# smallest case; result is thrown away.
function warmup_formulations(; case = "case14_ieee", variant = "std")
    for f in FORMULATIONS
        try
            run_centralized_one(case, variant, f; optimizer = _optimizer(f))
        catch
        end
    end
    @info "Warm-up done (discarded)" case formulations = FORMULATIONS
    return nothing
end

"""Probe whether Ipopt actually loads MA57 (returns true/false + note)."""
function probe_ma57()
    have_hsl() || return (false, "LIRION_HSLLIB unset")
    r = run_centralized_one("case14_ieee", "std", "ac"; optimizer = _ipopt_ac())
    r.status == "error" && return (false, r.notes)
    return (true, "")
end

"""Probe whether Gurobi (official DC solver) is reachable + solves DC."""
function probe_gurobi()
    _GUROBI === nothing &&
        return (false, "Gurobi.jl not loadable — set GUROBI_HOME and GUROBI_JL_USE_GUROBI_JLL=false")
    opt = _gurobi_dc()
    opt === nothing && return (false, "Gurobi env/license unavailable")
    r = run_centralized_one("case14_ieee", "std", "dc"; optimizer = opt)
    r.status == "error" && return (false, first(split(r.notes, "|")))
    return (true, "")
end

const OFFICIAL_NOTE =
    "timing excludes discarded warm-up solve; JULIA_NUM_THREADS=1; BLAS threads=1; " *
    "PowerModels validation (T0) is a separate stage"

function run_centralized_baseline()
    mode = get(ENV, "LIRION_CENTRALIZED_MODE", "smoke")
    cases = mode == "full" ? FULL_CASES : SMOKE_CASES

    ma57_ok, ma57_note = probe_ma57()
    gurobi_ok, gurobi_note = probe_gurobi()
    @info("Centralized baseline config",
        mode, n_cases = length(cases), formulations = FORMULATIONS,
        JULIA_NUM_THREADS = Threads.nthreads(),
        BLAS_threads = BLAS.get_num_threads(),
        LIRION_HSLLIB = hsllib_path(),
        ma57_loaded = ma57_ok, gurobi_loaded = gurobi_ok)
    if !gurobi_ok
        @warn("DC will NOT use Gurobi — falling back to Ipopt. This is NOT the " *
              "official DC configuration. Set GUROBI_HOME and " *
              "GUROBI_JL_USE_GUROBI_JLL=false to enable Gurobi.", reason = gurobi_note)
    end
    if !ma57_ok
        @warn("AC will NOT use MA57 — falling back to Ipopt/MUMPS. This is " *
              "NOT the official AC configuration. Set LIRION_HSLLIB to a valid " *
              "libhsl.so to enable MA57.", reason = ma57_note)
    end

    warmup_formulations()   # discarded; not part of T1 nor reported timings

    rows = NamedTuple[]
    for (name, variant) in cases
        for f in FORMULATIONS
            row = run_centralized_one(name, variant, f; optimizer = _optimizer(f))
            row = merge(row, (; notes = _with_note(row.notes, OFFICIAL_NOTE)))
            push!(rows, row)
            @info @sprintf("  %-16s %-3s %-6s -> %-14s obj=%s solve=%.3fs%s",
                name, variant, f, row.status,
                isnan(row.objective) ? "NaN" : @sprintf("%.4f", row.objective),
                isnan(row.solve_time_s) ? NaN : row.solve_time_s,
                isempty(row.notes) ? "" : "  [$(row.notes)]")
        end
    end

    raw = write_csv(joinpath(RES_CENTRAL, "centralized_baseline_raw.csv"), rows)
    t1c = write_csv(joinpath(TABLES_DIR, "T1_centralized_baseline.csv"), rows)
    t1m = write_markdown(joinpath(TABLES_DIR, "T1_centralized_baseline.md"), rows;
                         title = "T1 — Lirion centralized baseline ($(mode))")
    @info "Wrote" raw = relpath(raw, PAPER_DIR) csv = relpath(t1c, PAPER_DIR) md = relpath(t1m, PAPER_DIR)
    return rows
end

run_centralized_baseline()
