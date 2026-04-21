# ============================================================
# paper_results_centralized.jl
# ============================================================
# Centralised reference OPF solves for the Lirion paper benchmark.
# Writes objective values to `paper/cases/_centralized.toml` so that
# `paper_results_admm.jl` and `paper_results_aladin.jl` can compute
# gaps without repeating the centralised solves.
#
# Usage (after `using Lirion`):
#
#     julia> include(joinpath(pkgdir(Lirion), "paper", "paper_results_centralized.jl"))
#
# Auto-detects Gurobi / Ipopt / Clarabel / OSQP. Skips (case, variant,
# base) triples where the centralised formulation is infeasible (common
# for DC on SAD variants of several PGLib cases).
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

_HAVE_IPOPT || @warn "Ipopt not loaded — AC centralised results will be unavailable."

const TIME_LIMIT_S = 300.0

_gurobi_opt() = optimizer_with_attributes(
    () -> Gurobi.Optimizer(Gurobi.Env()),
    "OutputFlag"     => 0,
    "TimeLimit"      => TIME_LIMIT_S,
    "FeasibilityTol" => 1e-7,
    "OptimalityTol"  => 1e-7,
)

_solver_for_base(base::Symbol) = begin
    base === :dc    && return _HAVE_GUROBI ? _gurobi_opt() : nothing
    base === :ac    && return _HAVE_IPOPT  ? Ipopt.Optimizer : nothing
    base === :socwr && return _HAVE_GUROBI ? _gurobi_opt() :
                              (_HAVE_CLARABEL ? Clarabel.Optimizer : nothing)
    error("unknown base $base")
end

const PAPER_DIR  = @__DIR__
const CASES_DIR  = joinpath(PAPER_DIR, "cases")
const REPORT_CEN = joinpath(CASES_DIR, "_centralized.toml")

const CASES = ["case14_ieee", "case30_ieee", "case57_ieee",
               "case73_ieee_rts", "case118_ieee"]
const VARIANTS = ["std", "api", "sad"]

function _base_specs()
    specs = [(:dc, dc), (:socwr, socwr)]
    _HAVE_IPOPT && insert!(specs, 2, (:ac, ac))
    return specs
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

function _central_reference(case::String, variant::String, base::Symbol, base_ctor)
    # For DC, `_solver_for_base` returns `nothing` when Gurobi is
    # absent, which triggers Lirion's built-in tuned OSQP. Only treat
    # `nothing` as a hard failure for bases whose solver package isn't
    # loaded (AC → Ipopt, SOCWR → Clarabel/Gurobi).
    if base === :ac && !_HAVE_IPOPT
        return (status = "no_solver", objective = NaN, time_s = 0.0)
    elseif base === :socwr && !_HAVE_GUROBI && !_HAVE_CLARABEL
        return (status = "no_solver", objective = NaN, time_s = 0.0)
    end
    src = _original_file(case, variant)
    isfile(src) || return (status = "missing_file", objective = NaN, time_s = 0.0)
    opt = _solver_for_base(base)
    t0 = time()
    try
        r = opt === nothing ?
            Lirion.solve(src; model = base_ctor(), algorithm = centralized(),
                outputs = false, plots = false) :
            Lirion.solve(src; model = base_ctor(), algorithm = centralized(),
                optimizer = opt, outputs = false, plots = false)
        return (
            status    = _classify(r.status),
            objective = isfinite(r.objective) ? float(r.objective) : NaN,
            time_s    = time() - t0,
        )
    catch err
        return (status = "error", objective = NaN, time_s = time() - t0)
    end
end

function main()
    println("╔════════════════════════════════════════════════════════════╗")
    println("║  Lirion paper — CENTRALISED reference solves               ║")
    println("║  Cases:    ", join(CASES, ", "))
    println("║  Variants: ", join(VARIANTS, ", "))
    @printf("║  Solvers:  Gurobi=%s  Ipopt=%s  Clarabel=%s  OSQP=%s\n",
            _HAVE_GUROBI, _HAVE_IPOPT, _HAVE_CLARABEL, _HAVE_OSQP)
    println("║  Started:  ", now())
    println("║  Output:   ", REPORT_CEN)
    println("╚════════════════════════════════════════════════════════════╝")
    flush(stdout)

    isdir(CASES_DIR) || mkpath(CASES_DIR)
    report = isfile(REPORT_CEN) ? TOML.parsefile(REPORT_CEN) : Dict{String, Any}()

    for case in CASES, variant in VARIANTS
        haskey(report, case) || (report[case] = Dict{String, Any}())
        haskey(report[case], variant) || (report[case][variant] = Dict{String, Any}())
        for (base, base_ctor) in _base_specs()
            bkey = string(base)
            if haskey(report[case][variant], bkey) &&
               haskey(report[case][variant][bkey], "status")
                rec = report[case][variant][bkey]
                @printf("[%s/%-3s/%-5s] cached  %-12s obj=%s\n",
                        case, variant, base, rec["status"],
                        isfinite(get(rec, "objective", NaN)) ?
                            @sprintf("%12.2f", rec["objective"]) : "NaN")
                flush(stdout)
                continue
            end
            res = _central_reference(case, variant, base, base_ctor)
            report[case][variant][bkey] = Dict(
                "status"    => res.status,
                "objective" => isfinite(res.objective) ? res.objective : 0.0,
                "time_s"    => res.time_s,
            )
            _save_toml!(report, REPORT_CEN)
            @printf("[%s/%-3s/%-5s] %-12s obj=%s t=%.2fs\n",
                    case, variant, base, res.status,
                    isfinite(res.objective) ? @sprintf("%12.2f", res.objective) : "NaN",
                    res.time_s)
            flush(stdout)
        end
    end

    report["_meta"] = Dict(
        "generated_at" => string(now()),
        "script"       => "paper_results_centralized.jl",
    )
    _save_toml!(report, REPORT_CEN)
    println("\nCentralised references written to: ", REPORT_CEN)
    return report
end

main()
