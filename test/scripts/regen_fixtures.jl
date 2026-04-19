# ============================================================
# Regenerate test fixtures.
#
# For each (case_file, model) pair in the test matrix:
#   1. Solve with Lirion.
#   2. Solve with PowerModels (reference).
#   3. If the gap is within tolerance, write Lirion's objective to
#      `test/fixtures/<group>/<case>.toml` under [objectives].<model>.
#
# Run:  julia --project=. test/scripts/regen_fixtures.jl
# ============================================================

using Pkg
Pkg.activate(@__DIR__)
Pkg.develop(; path = dirname(dirname(@__DIR__)))
Pkg.instantiate()

using TOML
using Lirion
using PowerModels
using Ipopt
using Clarabel
using JuMP: optimizer_with_attributes

const CASES_DIR = joinpath(dirname(@__DIR__), "matpower_cases")
const FIXTURES_DIR = joinpath(dirname(@__DIR__), "fixtures", "centralized")

const BASE_CASES = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case5_pjm.m",
    "pglib_opf_case14_ieee.m",
    "pglib_opf_case24_ieee_rts.m",
    "pglib_opf_case30_as.m",
    "pglib_opf_case30_ieee.m",
    "pglib_opf_case39_epri.m",
    "pglib_opf_case57_ieee.m",
    "pglib_opf_case60_c.m",
    "pglib_opf_case73_ieee_rts.m",
    "pglib_opf_case118_ieee.m",
]

const STORAGE_CASES = ["case14_bus_4_bat.m", "case57_bus_6_bat.m"]

const IPOPT_QUIET = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
const CLARABEL_QUIET = optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false)

# Model matrix per case group. SDPWRM excluded past case30_as because
# Clarabel hits NUMERICAL_ERROR on larger dense PSD problems.
const SMALL_CASES = Set([
    "pglib_opf_case3_lmbd.m", "pglib_opf_case5_pjm.m",
    "pglib_opf_case14_ieee.m", "pglib_opf_case24_ieee_rts.m",
    "pglib_opf_case30_as.m",
])

function _lirion_obj(file::AbstractString, model_ctor, optimizer)
    r = Lirion.solve(
        file; model = model_ctor(), algorithm = Lirion.centralized(),
        optimizer = optimizer, outputs = false, plots = false,
    )
    return Float64(r.objective), r.status
end

function _pm_obj(file::AbstractString, tag::Symbol, optimizer, has_storage::Bool)
    pm = PowerModels.parse_file(file)
    tag === :IVR &&
        return Float64(PowerModels.solve_opf_iv(pm, IVRPowerModel, optimizer)["objective"])

    if has_storage
        pm_model =
            tag === :DC ? DCPPowerModel :
            tag === :AC ? ACPPowerModel :
            tag === :ACR ? ACRPowerModel :
            tag === :LPAC ? LPACCPowerModel :
            tag === :SOCWR ? SOCWRPowerModel :
            tag === :SOCBF ? SOCBFPowerModel :
            error("no storage PM model for $tag")
        return Float64(PowerModels._solve_opf_strg(pm, pm_model, optimizer)["objective"])
    end

    pm_model =
        tag === :DC ? DCPPowerModel :
        tag === :AC ? ACPPowerModel :
        tag === :ACR ? ACRPowerModel :
        tag === :LPAC ? LPACCPowerModel :
        tag === :SOCWR ? SOCWRConicPowerModel :
        tag === :QC ? QCRMPowerModel :
        tag === :SDPWRM ? SDPWRMPowerModel :
        tag === :SOCBF ? SOCBFConicPowerModel :
        tag === :BFA ? BFAPowerModel :
        tag === :NFA ? NFAPowerModel :
        error("no PM model for $tag")

    if tag === :AC
        return Float64(PowerModels.solve_ac_opf(pm, optimizer)["objective"])
    elseif tag in (:SOCBF, :BFA)
        return Float64(PowerModels.solve_opf_bf(pm, pm_model, optimizer)["objective"])
    else
        return Float64(PowerModels.solve_opf(pm, pm_model, optimizer)["objective"])
    end
end

function _models_for(case::AbstractString)
    has_storage = occursin("_bat", case)
    if has_storage
        return [
            (:DC, Lirion.dc, IPOPT_QUIET),
            (:AC, Lirion.ac, IPOPT_QUIET),
            (:ACR, Lirion.acr, IPOPT_QUIET),
        ]
    end
    base = [
        (:DC, Lirion.dc, IPOPT_QUIET),
        (:AC, Lirion.ac, IPOPT_QUIET),
        (:ACR, Lirion.acr, IPOPT_QUIET),
        (:IVR, Lirion.ivr, IPOPT_QUIET),
        (:LPAC, Lirion.lpac, IPOPT_QUIET),
        (:SOCWR, Lirion.socwr, CLARABEL_QUIET),
        (:QC, Lirion.qc, IPOPT_QUIET),
        (:SOCBF, Lirion.socbf, CLARABEL_QUIET),
        (:BFA, Lirion.bfa, CLARABEL_QUIET),
    ]
    if case in SMALL_CASES
        push!(base, (:SDPWRM, Lirion.sdpwrm, CLARABEL_QUIET))
    end
    return base
end

function _pm_optimizer_for(tag::Symbol)
    return tag in (:SOCWR, :SOCBF, :BFA, :SDPWRM) ? CLARABEL_QUIET : IPOPT_QUIET
end

function regen_case(case::AbstractString; pm_tol::Real = 0.05)
    file = joinpath(CASES_DIR, case)
    isfile(file) || (println("  SKIP (not found): $case"); return nothing)

    has_storage = occursin("_bat", case)
    group_dir = FIXTURES_DIR
    mkpath(group_dir)

    fx_path = joinpath(group_dir, replace(case, ".m" => ".toml"))
    fixture = Dict{String, Any}(
        "metadata" => Dict{String, Any}(
            "case" => case, "has_storage" => has_storage,
            "regenerated_at" => string(Dates.today()),
        ),
        "objectives" => Dict{String, Any}(),
        "pm_reference" => Dict{String, Any}(),
        "status" => Dict{String, Any}(),
    )

    for (tag, ctor, opt) in _models_for(case)
        print("  $case | $tag … ")
        lirion_obj, lirion_status = try
            _lirion_obj(file, ctor, opt)
        catch e
            println("LIRION FAILED: $(typeof(e))")
            continue
        end
        pm_obj = try
            _pm_obj(file, tag, _pm_optimizer_for(tag), has_storage)
        catch e
            println("PM FAILED: $(typeof(e)) (keeping Lirion obj anyway)")
            NaN
        end

        if !isnan(pm_obj)
            gap = 100.0 * abs(lirion_obj - pm_obj) / max(abs(pm_obj), 1e-6)
            if gap > 100 * pm_tol
                println(
                    "SKIP (gap $(round(gap, digits=3))% > $(round(100*pm_tol, digits=2))%)",
                )
                continue
            end
            println(
                "OK  lirion=$(round(lirion_obj,digits=4)) pm=$(round(pm_obj,digits=4)) gap=$(round(gap,digits=4))%",
            )
            fixture["pm_reference"][lowercase(String(tag))] = pm_obj
        else
            println("OK  lirion=$(round(lirion_obj,digits=4)) (no PM ref)")
        end
        fixture["objectives"][lowercase(String(tag))] = lirion_obj
        fixture["status"][lowercase(String(tag))] = string(lirion_status)
    end

    open(fx_path, "w") do io
        return TOML.print(io, fixture; sorted = true)
    end
    println("  → wrote $fx_path")
    return fx_path
end

function main()
    @info "Regenerating centralized fixtures …"
    for case in BASE_CASES
        regen_case(case; pm_tol = 0.05)
    end
    for case in STORAGE_CASES
        regen_case(case; pm_tol = 0.05)
    end
    @info "Done. Fixtures in $(FIXTURES_DIR)"
    return nothing
end

using Dates
main()
