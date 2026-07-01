# ============================================================
# paper/scripts/experiments/lib/experiment_common.jl
# ============================================================
# Shared helpers for the paper's experiment drivers. Stage A only needs
# the centralized pieces (paths, model/solver mapping, a single-run
# wrapper, CSV/Markdown writers). No ADMM/ALADIN/tuning here.
# ============================================================

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
    f == "dc"    && return Lirion.dc
    f == "ac"    && return Lirion.ac
    f == "socwr" && return Lirion.socwr
    error("unknown formulation: $f")
end

# HSL/MA57 for AC is opt-in via LIRION_HSLLIB (a path to a licensed
# libhsl.so). No binary/license ever lives in the repo.
hsllib_path() = get(ENV, "LIRION_HSLLIB", "")
have_hsl() = !isempty(hsllib_path()) && isfile(hsllib_path())

"""Human-readable solver label per formulation (goes into the `solver` column)."""
function solver_label(f::AbstractString)
    f == "socwr" && return "Clarabel (hardened)"
    f == "ac"    && return have_hsl() ? "Ipopt + MA57 via LIRION_HSLLIB" : "Ipopt (MUMPS)"
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
        pm_objective = "", pm_gap_pct = "", pm_status = "",   # T0 deferred
        notes = _with_note(note, formulation_note(formulation)),
    )
end

# Append a note without clobbering an existing (e.g. error) note.
function _with_note(existing, extra)
    isempty(extra) && return existing
    isempty(existing) && return extra
    return string(existing, " | ", extra)
end

# Per-formulation methodological note (documents the solver choice).
# No LicenseID or license paths — only what is needed to reproduce.
function formulation_note(f::AbstractString)
    f == "dc"    && return "DC solved with Gurobi 12.0.3 via GUROBI_HOME " *
        "(native LP/QP solver for the convex DC-OPF)"
    f == "ac"    && return "AC solved with Ipopt + MA57 via LIRION_HSLLIB"
    f == "socwr" && return "SOCWR solved with hardened Clarabel settings"
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
