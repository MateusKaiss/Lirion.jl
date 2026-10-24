# Distributed-run helpers: partitioned-case paths, oracle loading, one run +
# error classification, result schema.

include(joinpath(@__DIR__, "experiment_common.jl"))

using DataFrames: nrow
using TOML
import Dates

const PARTITIONED_DIR = joinpath(PAPER_DIR, "cases", "partitioned")
const RES_DIST = joinpath(RESULTS_DIR, "distributed")   # scratch (git-ignored)

function write_toml(path::AbstractString, data::AbstractDict)
    mkpath(dirname(path))
    open(path, "w") do io
        TOML.print(io, data; sorted = true)
    end
    return path
end

"""Path to the partitioned case.m for (case, variant, k)."""
partitioned_case_path(case, variant, k) =
    joinpath(PARTITIONED_DIR, case, variant, "k$(k)", "case.m")

# zone solver per formulation when optimizer = nothing
dist_solver_label(f) =
    f == "dc"      ? "OSQP (Lirion distributed default)" :
    f == "socwr"   ? "Clarabel (hardened, Lirion default)" :
    f == "socwrnl" ? "Ipopt+MA57 @ dual_inf_tol=1e-9 (Lirion _socwrnl_tight_factory)" :
    "Ipopt/MUMPS (Lirion distributed default)"   # ac

"""Load the T1 centralized objective per (case, formulation) for variant std.
The T1 CSV has no commas inside fields, so a plain split is safe."""
function load_t1_objectives()
    d = Dict{Tuple{String,String},Float64}()
    path = joinpath(TABLES_DIR, "T1_centralized_baseline.csv")
    isfile(path) || return d
    lines = readlines(path)
    for ln in lines[2:end]
        f = split(ln, ",")
        length(f) >= 8 || continue
        obj = tryparse(Float64, f[8])
        obj === nothing && continue
        d[(String(f[1]), String(f[3]))] = obj   # (case, formulation) => objective
    end
    return d
end

# one generic rho for the smoke runs
const SMOKE_RHO = 1000.0

# ---- Algorithm builder -----------------------------------------------
"""Build an ADMM/ALADIN algorithm; may throw on invalid combos (caught
upstream and classified)."""
function build_algorithm(algorithm, sectioning, coupling;
                         proximal_mode = nothing,
                         admm_max_iter = 300, aladin_max_iter = 40)
    sect = Symbol(sectioning)
    coup = Symbol(coupling)
    if algorithm == "admm"
        return Lirion.admm(; coupling = coup, sectioning = sect,
                           max_iter = admm_max_iter, rho = SMOKE_RHO)
    elseif algorithm == "aladin"
        kw = (; coupling = coup, sectioning = sect, max_iter = aladin_max_iter,
              rho = SMOKE_RHO)
        proximal_mode !== nothing &&
            (kw = merge(kw, (; proximal_mode = Symbol(proximal_mode))))
        return Lirion.aladin(; kw...)
    end
    error("unknown algorithm: $algorithm")
end

# ---- Residuals -------------------------------------------------------
function extract_residuals(r)
    pr = NaN; dr = NaN
    d = r.diagnostics
    if hasproperty(d, :r_primal) && !isempty(d.r_primal); pr = last(d.r_primal); end
    if hasproperty(d, :r_dual) && !isempty(d.r_dual); dr = last(d.r_dual); end
    if isnan(pr) && r.consensus !== nothing && :residual in propertynames(r.consensus) &&
       nrow(r.consensus) > 0
        pr = maximum(abs.(r.consensus.residual))
    end
    return pr, dr
end

# ---- Error classification -------------------------------------------
function classify_error(e)
    msg = first(split(sprint(showerror, e), "\n"))
    low = lowercase(msg)
    if e isa MethodError || occursin("no method matching", low) ||
       occursin("not implemented", low) || occursin("not supported", low) ||
       occursin("unsupported", low)
        return ("not_implemented", msg)         # capability gap (incl. by-design exclusions)
    elseif e isa ArgumentError &&
           (occursin("proximal_mode", low) || occursin("pflow", low) ||
            occursin("coupling=", low))
        return ("skipped_by_design", msg)       # API forbids this combination
    else
        return ("unstable_or_fails", msg)       # ran but failed numerically
    end
end

# ---- One distributed run --------------------------------------------
const DIST_COLUMNS = [
    :case, :variant, :k, :formulation, :algorithm, :sectioning, :coupling,
    :status, :objective, :centralized_objective, :gap_pct, :iterations,
    :primal_residual, :dual_residual, :runtime_s, :solver, :parameter_set,
    :convergence_flag, :failure_reason, :support_status, :notes,
    :case_path, :timestamp,
]

function run_distributed_one(case, variant, k, formulation, algorithm, sectioning,
                             coupling; proximal_mode = nothing, t1 = Dict(),
                             admm_max_iter = 300, aladin_max_iter = 40)
    path = partitioned_case_path(case, variant, k)
    ctor = model_ctor(formulation)
    cen = get(t1, (case, formulation), NaN)

    status = ""; obj = NaN; iters = 0; pr = NaN; dr = NaN
    runtime = NaN; conv = false; support = ""; fail = ""; note = ""
    coup_label = proximal_mode === nothing ? coupling : "$(coupling)+$(proximal_mode)"

    t0 = time()
    try
        algo = build_algorithm(algorithm, sectioning, coupling;
                               proximal_mode = proximal_mode,
                               admm_max_iter = admm_max_iter,
                               aladin_max_iter = aladin_max_iter)
        r = Lirion.solve(path; model = ctor(), algorithm = algo,
                         optimizer = nothing, outputs = false, plots = false)
        runtime = time() - t0
        obj = Float64(r.objective)
        status = string(r.status)
        iters = r.iter_done === nothing ? 0 : r.iter_done
        pr, dr = extract_residuals(r)
        conv = (r.status == :converged)
        if isfinite(obj)
            support = "supported"
            conv || (note = "did not reach tolerance within smoke max_iter")
        else
            support = "unstable_or_fails"; fail = "non-finite objective"
        end
    catch e
        runtime = time() - t0
        support, fail = classify_error(e)
        status = "error"
    end

    gap = (isnan(obj) || isnan(cen) || cen == 0) ? NaN :
          100 * (obj - cen) / abs(cen)

    return (
        case = case, variant = variant, k = k, formulation = formulation,
        algorithm = algorithm, sectioning = sectioning, coupling = coup_label,
        status = status, objective = obj, centralized_objective = cen,
        gap_pct = gap, iterations = iters, primal_residual = pr,
        dual_residual = dr, runtime_s = runtime,
        solver = dist_solver_label(formulation), parameter_set = "smoke_default (rho=1000)",
        convergence_flag = conv, failure_reason = fail,
        support_status = support, notes = note,
        case_path = relpath(path, PAPER_DIR), timestamp = string(Dates.now()),
    )
end
