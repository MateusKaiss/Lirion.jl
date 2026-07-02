# ============================================================
# paper/scripts/experiments/run_runtime_alignment.jl
# ============================================================
# Stage C0 — runtime/solver alignment (DIAGNOSTIC ONLY, no tuning).
# Decides the official runtime config before the main experiments:
#   - BLAS threads sweep (JULIA_NUM_THREADS held at 1);
#   - AC distributed: Ipopt/MUMPS (default) vs Ipopt+MA57 (LIRION_HSLLIB);
#   - DC distributed: OSQP (Lirion default) vs Gurobi probe;
#   - SOCWR distributed: confirm hardened Clarabel (Lirion default).
#
# Small representative set only. Objectives must not change with BLAS
# threads (sanity). Outputs:
#   results/runtime_alignment/runtime_alignment_raw.csv  (scratch)
#   artifacts/tables/C0_runtime_alignment.csv / .md      (curated)
#   artifacts/runtime_config.toml                        (curated)
# ============================================================

include(joinpath(@__DIR__, "lib", "distributed_common.jl"))

using LinearAlgebra
import MathOptInterface as MOI
using JuMP: optimizer_with_attributes
import Dates

const RES_RTA = joinpath(RESULTS_DIR, "runtime_alignment")

# ---- Solver factories (mirror the official T1 config) ----------------
_gurobi_env = nothing
const _GUROBI = try
    @eval import Gurobi
    Gurobi
catch
    nothing
end
function gurobi_opt()
    _GUROBI === nothing && return nothing
    try
        global _gurobi_env
        _gurobi_env === nothing && (_gurobi_env = _GUROBI.Env())
        return optimizer_with_attributes(() -> _GUROBI.Optimizer(_gurobi_env), "OutputFlag" => 0)
    catch
        return nothing
    end
end

ipopt_ac_central() = have_hsl() ?
    optimizer_with_attributes(Lirion.Ipopt.Optimizer, "print_level" => 0,
        "linear_solver" => "ma57", "hsllib" => hsllib_path()) :
    optimizer_with_attributes(Lirion.Ipopt.Optimizer, "print_level" => 0)

# AC distributed optimizer: clone Lirion's tuned _SILENT_FACTORY_AC and
# swap only the linear solver (keeps all convergence tuning intact).
function ac_dist_optimizer(; ma57::Bool)
    base = Lirion._SILENT_FACTORY_AC
    keep = Pair[]
    for p in base.params
        nm = p.first isa MOI.RawOptimizerAttribute ? p.first.name : ""
        (nm == "linear_solver" || nm == "hsllib") && continue
        push!(keep, p)
    end
    if ma57 && have_hsl()
        push!(keep, MOI.RawOptimizerAttribute("linear_solver") => "ma57")
        push!(keep, MOI.RawOptimizerAttribute("hsllib") => hsllib_path())
    else
        push!(keep, MOI.RawOptimizerAttribute("linear_solver") => "mumps")
    end
    return optimizer_with_attributes(base.optimizer_constructor, keep...)
end

# ---- Row helpers -----------------------------------------------------
const C0_COLUMNS = [
    :mode, :case, :formulation, :algorithm, :sectioning, :coupling,
    :solver_variant, :blas_threads, :julia_threads, :status, :objective,
    :obj_match_ref, :iterations, :runtime_s, :warnings, :notes,
]

function _central_row(case, f, solver_variant, opt, blas, refobj)
    ctor = model_ctor(f); status = "error"; obj = NaN; rt = NaN; note = ""
    t0 = time()
    try
        r = Lirion.solve(raw_case_path(case, "std"); model = ctor(),
            algorithm = Lirion.centralized(), optimizer = opt,
            outputs = false, plots = false)
        rt = time() - t0; obj = Float64(r.objective); status = string(r.status)
    catch e
        rt = time() - t0; note = first(split(sprint(showerror, e), "\n"))
    end
    m = (isnan(obj) || isnan(refobj)) ? "" :
        string(isapprox(obj, refobj; rtol = 1e-6))
    return (mode = "centralized", case = case, formulation = f, algorithm = "centralized",
        sectioning = "", coupling = "", solver_variant = solver_variant,
        blas_threads = blas, julia_threads = Threads.nthreads(), status = status,
        objective = obj, obj_match_ref = m, iterations = 0, runtime_s = rt,
        warnings = "", notes = note)
end

function _dist_row(case, k, f, solver_variant, opt, blas; sect = "tie_line", coup = "hybrid",
                   algo = "admm", refobj = NaN)
    ctor = model_ctor(f); status = "error"; obj = NaN; rt = NaN; iters = 0; note = ""
    t0 = time()
    try
        alg = build_algorithm(algo, sect, coup; admm_max_iter = 300, aladin_max_iter = 40)
        r = Lirion.solve(partitioned_case_path(case, "std", k); model = ctor(),
            algorithm = alg, optimizer = opt, outputs = false, plots = false)
        rt = time() - t0; obj = Float64(r.objective); status = string(r.status)
        iters = r.iter_done === nothing ? 0 : r.iter_done
    catch e
        rt = time() - t0; note = first(split(sprint(showerror, e), "\n"))
    end
    m = (isnan(obj) || isnan(refobj)) ? "" :
        string(isapprox(obj, refobj; rtol = 1e-6))
    return (mode = "distributed", case = case, formulation = f, algorithm = algo,
        sectioning = sect, coupling = coup, solver_variant = solver_variant,
        blas_threads = blas, julia_threads = Threads.nthreads(), status = status,
        objective = obj, obj_match_ref = m, iterations = iters, runtime_s = rt,
        warnings = "", notes = note)
end

function run_runtime_alignment()
    mkpath(RES_RTA)
    BLAS.set_num_threads(1)
    blas_list = [1, 2, 4, 8, 14]
    gopt = gurobi_opt()
    @info "C0 runtime alignment config" JULIA_NUM_THREADS = Threads.nthreads() gurobi = (gopt !== nothing) ma57 = have_hsl() blas_list

    rows = NamedTuple[]

    # warm-up (discard): one central + one distributed per formulation
    for f in FORMULATIONS
        try; _central_row("case500_goc", f, "warmup", f == "dc" ? gopt : (f == "socwr" ? nothing : ipopt_ac_central()), 1, NaN); catch; end
        try; _dist_row("case57_ieee", 2, f, "warmup", nothing, 1); catch; end
    end
    @info "Warm-up done (discarded)"

    # ---- Centralized BLAS sweep (case1354, official solvers) ----
    central_opt(f) = f == "dc" ? gopt : (f == "socwr" ? nothing : ipopt_ac_central())
    solver_name(f) = f == "dc" ? "Gurobi" : (f == "socwr" ? "Clarabel(hardened)" : (have_hsl() ? "Ipopt+MA57" : "Ipopt+MUMPS"))
    refc = Dict{String,Float64}()
    for f in FORMULATIONS
        for b in blas_list
            BLAS.set_num_threads(b)
            row = _central_row("case1354_pegase", f, solver_name(f), central_opt(f), b,
                               get(refc, f, NaN))
            haskey(refc, f) || (refc[f] = row.objective)
            push!(rows, row)
            @info @sprintf("  [C] case1354 %-6s %-16s BLAS=%-2d -> %-14s obj=%s rt=%.3fs match=%s",
                f, row.solver_variant, b, row.status,
                isnan(row.objective) ? "NaN" : @sprintf("%.2f", row.objective), row.runtime_s, row.obj_match_ref)
        end
    end

    # ---- Distributed BLAS sweep (case118 k2, ADMM tie_line hybrid, Lirion defaults) ----
    refd = Dict{String,Float64}()
    for f in FORMULATIONS
        for b in blas_list
            BLAS.set_num_threads(b)
            row = _dist_row("case118_ieee", 2, f, dist_solver_label(f), nothing, b;
                            refobj = get(refd, f, NaN))
            haskey(refd, f) || (refd[f] = row.objective)
            push!(rows, row)
            @info @sprintf("  [D] case118 %-6s BLAS=%-2d -> %-16s it=%-4d rt=%.3fs match=%s",
                f, b, row.status, row.iterations, row.runtime_s, row.obj_match_ref)
        end
    end

    # ---- AC distributed: MUMPS vs MA57 (BLAS=1) ----
    BLAS.set_num_threads(1)
    for case in ["case57_ieee", "case118_ieee"]
        ref = _dist_row(case, 2, "ac", "Ipopt+MUMPS(dist)", ac_dist_optimizer(ma57 = false), 1)
        push!(rows, ref)
        push!(rows, _dist_row(case, 2, "ac", "Ipopt+MA57(dist)", ac_dist_optimizer(ma57 = true), 1;
                              refobj = ref.objective))
    end

    # ---- DC distributed: OSQP (default) vs Gurobi probe (BLAS=1) ----
    ref = _dist_row("case57_ieee", 2, "dc", "OSQP(default)", nothing, 1)
    push!(rows, ref)
    push!(rows, _dist_row("case57_ieee", 2, "dc", "Gurobi(probe)", gopt, 1; refobj = ref.objective))

    # ---- write ----
    raw = write_csv(joinpath(RES_RTA, "runtime_alignment_raw.csv"), rows; columns = C0_COLUMNS)
    csv = write_csv(joinpath(TABLES_DIR, "C0_runtime_alignment.csv"), rows; columns = C0_COLUMNS)
    md  = write_markdown(joinpath(TABLES_DIR, "C0_runtime_alignment.md"), rows;
                         columns = C0_COLUMNS, title = "C0 — runtime/solver alignment")
    write_toml(joinpath(ARTIFACTS_DIR, "runtime_config.toml"), Dict(
        "_meta" => Dict("generated_at" => string(Dates.now()), "stage" => "C0",
            "policy" => "Solver config is decided per ROLE. Centralized = objective " *
                "oracle; distributed = the study object (algorithm × sectioning × " *
                "coupling × k). The two roles may use different solvers when that " *
                "reflects their role; both are fixed rules per category (no " *
                "case-by-case cherry-picking). Within the distributed results, ADMM " *
                "and ALADIN share the SAME solver per formulation (guaranteed: " *
                "resolve_optimizer dispatches on base_model, not on the algorithm)."),
        "runtime_config" => Dict("JULIA_NUM_THREADS" => 1, "BLAS_threads" => 1),
        "centralized_solver_config" => Dict(
            "dc" => "Gurobi", "socwr" => "Clarabel (hardened)",
            "ac" => "Ipopt + MA57 via LIRION_HSLLIB"),
        "distributed_solver_config" => Dict(
            "dc" => "OSQP (Lirion default; fast warm-started inner QP)",
            "socwr" => "Clarabel (hardened, Lirion default)",
            "ac" => "Ipopt + MA57 via LIRION_HSLLIB (falls back to MUMPS if unset)"),
    ))
    @info "Wrote" raw = relpath(raw, PAPER_DIR) csv = relpath(csv, PAPER_DIR) md = relpath(md, PAPER_DIR)
    return rows
end

run_runtime_alignment()
