# Shared helpers: paths, model/solver mapping, one centralized run,
# CSV/Markdown writers.

using Lirion
using DataFrames: nrow
using Printf
using Dates
import JuMP

# ---- Paths ------------------------------------------------------------
const EXP_LIB_DIR   = @__DIR__
const PAPER_DIR     = normpath(joinpath(EXP_LIB_DIR, "..", "..", ".."))
const RAW_DIR       = joinpath(PAPER_DIR, "cases", "raw", "pglib")
const RESULTS_DIR   = joinpath(PAPER_DIR, "results")
const RES_CENTRAL   = joinpath(RESULTS_DIR, "centralized")   # scratch (git-ignored)
const ARTIFACTS_DIR = joinpath(PAPER_DIR, "artifacts")
const TABLES_DIR    = joinpath(ARTIFACTS_DIR, "tables")       # curated (versioned)

"""Path to a copied raw PGLib case for the centralized baseline."""
function raw_case_path(name::AbstractString, variant::AbstractString)
    if variant == "std"
        return joinpath(RAW_DIR, "std", "pglib_opf_$(name).m")
    else
        return joinpath(RAW_DIR, variant, "pglib_opf_$(name)__$(variant).m")
    end
end

# ---- Formulations, model constructors, solvers ------------------------
const FORMULATIONS = ["dc", "socwr", "ac"]   # main-paper order

function model_ctor(f::AbstractString)
    f == "dc"      && return Lirion.dc
    f == "ac"      && return Lirion.ac
    f == "socwr"   && return Lirion.socwr
    f == "socwrnl" && return Lirion.socwrnl
    error("unknown formulation: $f")
end

# MA57 is opt-in via LIRION_HSLLIB; no binary in the repo
hsllib_path() = get(ENV, "LIRION_HSLLIB", "")
have_hsl() = !isempty(hsllib_path()) && isfile(hsllib_path())

"""
    ensure_hsl!(; require=false)

Guarantee AC runs on Ipopt+MA57. If LIRION_HSLLIB is unset/invalid, auto-detect a
licensed libhsl.so (x86_64-linux-gnu preferred) under the usual local install
dirs and set ENV["LIRION_HSLLIB"]. With `require=true`, ERROR if MA57 can't be
found — so the final result/verify runs never silently fall back to MUMPS.
Returns whether MA57 is available. (Path is machine-local; nothing is committed.)
"""
function ensure_hsl!(; require::Bool=false)
    if !have_hsl()
        roots = [joinpath(homedir(), "Projects", "hsl"),
                 joinpath(homedir(), ".julia", "artifacts")]
        cands = String[]
        for r in roots
            isdir(r) || continue
            for (dir, _, files) in walkdir(r)
                for fn in files
                    fn == "libhsl.so" && push!(cands, joinpath(dir, fn))
                end
            end
        end
        # prefer the native x86_64-linux-gnu build
        pick = something(findfirst(p -> occursin("x86_64-linux-gnu", p), cands), 0)
        chosen = pick > 0 ? cands[pick] : (isempty(cands) ? "" : first(cands))
        isempty(chosen) || (ENV["LIRION_HSLLIB"] = chosen)
    end
    ok = have_hsl()
    ok && @info "HSL/MA57 active" libhsl=hsllib_path()
    require && !ok && error("MA57 required but no libhsl.so found — set LIRION_HSLLIB. Refusing to run AC on MUMPS for final results.")
    ok || @warn "MA57 unavailable — AC will use MUMPS (same optimum, slower). Set LIRION_HSLLIB to enable MA57."
    return ok
end

"""Human-readable solver label per formulation (goes into the `solver` column)."""
function solver_label(f::AbstractString)
    f == "socwr"   && return "Clarabel (hardened)"
    # socwrnl oracle = bare Ipopt; the zones use Ipopt+MA57 (see the run manifest)
    f == "socwrnl" && return "Ipopt (SOCWR smooth-NLP centralized oracle)"
    f == "ac"      && return have_hsl() ? "Ipopt + MA57 via LIRION_HSLLIB" : "Ipopt (MUMPS)"
    return "Gurobi (official DC baseline)"   # dc
end

# ---- One centralized run ---------------------------------------------
"""
Solve one (case, variant, formulation) centrally and return a schema row
(NamedTuple). Never throws — solver/build errors are captured in `notes`
with `status = "error"`.
"""
function run_centralized_one(name, variant, formulation; optimizer)
    path = raw_case_path(name, variant)
    ctor = model_ctor(formulation)
    status = "error"; termstatus = ""; obj = NaN
    solvet = NaN; nbus = 0; nbr = 0; ngen = 0; note = ""

    t0 = time()
    try
        r = Lirion.solve(path; model = ctor(), algorithm = Lirion.centralized(),
                          optimizer = optimizer, outputs = false, plots = false)
        obj = Float64(r.objective)
        status = string(r.status)
        solvet = Float64(r.solve_time)
        termstatus = try
            string(JuMP.termination_status(r.model))
        catch
            ""
        end
        nbus = nrow(r.bus); nbr = nrow(r.branch); ngen = nrow(r.gen)
    catch e
        note = first(split(sprint(showerror, e), "\n"))
        status = "error"
    end
    total = time() - t0
    buildt = max(total - (isnan(solvet) ? 0.0 : solvet), 0.0)

    return (
        case = name, variant = variant, formulation = formulation,
        model = uppercase(formulation), solver = solver_label(formulation),
        status = status, termination_status = termstatus,
        objective = obj, solve_time_s = solvet, build_time_s = buildt,
        total_time_s = total, num_buses = nbus, num_branches = nbr,
        num_gens = ngen, case_path = relpath(path, PAPER_DIR),
        pm_objective = "", pm_gap_pct = "", pm_status = "",   # unused
        notes = _with_note(note, formulation_note(formulation)),
    )
end

# Append a note without clobbering an existing (e.g. error) note.
function _with_note(existing, extra)
    isempty(extra) && return existing
    isempty(existing) && return extra
    return string(existing, " | ", extra)
end

# per-formulation solver note (no license info)
function formulation_note(f::AbstractString)
    f == "dc"    && return "DC solved with Gurobi 13.0.2 via Gurobi_jll artifact " *
        "(native LP/QP solver for the convex DC-OPF; 12.0.3->13.0.2 after the " *
        "2026-07-31 license re-issue, benign — DC objective verified)"
    f == "ac"    && return "AC solved with Ipopt + MA57 via LIRION_HSLLIB"
    f == "socwr" && return "SOCWR solved with hardened Clarabel settings"
    f == "socwrnl" && return "SOCWRNL (same W-space SOCWR relaxation as a smooth " *
        "quadratic NLP) solved centrally with Ipopt; distributed zones use " *
        "Ipopt+MA57 at dual_inf_tol=1e-9 (exact Hessian) for Tier-2 measurability"
    return ""
end

# ---- CSV / Markdown writers ------------------------------------------
const T1_COLUMNS = [
    :case, :variant, :formulation, :model, :solver, :status,
    :termination_status, :objective, :solve_time_s, :build_time_s,
    :total_time_s, :num_buses, :num_branches, :num_gens, :case_path,
    :pm_objective, :pm_gap_pct, :pm_status, :notes,
]

_fmt(x::AbstractFloat) = isnan(x) ? "" : @sprintf("%.10g", x)
_fmt(x) = string(x)

function _csv_cell(x)
    s = _fmt(x)
    (occursin(',', s) || occursin('"', s) || occursin('\n', s)) &&
        return '"' * replace(s, '"' => "\"\"") * '"'
    return s
end

function write_csv(path, rows; columns = T1_COLUMNS)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(string.(columns), ","))
        for r in rows
            println(io, join((_csv_cell(getproperty(r, c)) for c in columns), ","))
        end
    end
    return path
end

function write_markdown(path, rows; columns = T1_COLUMNS, title = "")
    mkpath(dirname(path))
    open(path, "w") do io
        isempty(title) || println(io, "# $title\n")
        println(io, "| " * join(string.(columns), " | ") * " |")
        println(io, "|" * repeat(" --- |", length(columns)))
        for r in rows
            println(io, "| " * join((_fmt(getproperty(r, c)) for c in columns), " | ") * " |")
        end
    end
    return path
end
