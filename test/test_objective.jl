# ==============================================
# Lirion vs PowerModels — objective comparison (API solve)
# ==============================================
# ------------------------------------------------------------
# Short list of cases for rapid testing (comment/uncomment)
# ------------------------------------------------------------
const reduced_list = [
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
    # "pglib_opf_case89_pegase.m",
    # "pglib_opf_case162_ieee_dtc.m",
    # "pglib_opf_case179_goc.m",
    # "pglib_opf_case197_snem.m",
    # "pglib_opf_case200_activ.m",
    # "pglib_opf_case240_pserc.m",
    # "pglib_opf_case300_ieee.m",
    # "pglib_opf_case500_goc.m",
    # "pglib_opf_case588_sdet.m",
    # "pglib_opf_case793_goc.m",
    # "pglib_opf_case1354_pegase.m",
    # "pglib_opf_case1803_snem.m",
    # "pglib_opf_case1888_rte.m",
    # "pglib_opf_case1951_rte.m",
    # "pglib_opf_case2000_goc.m",
    # "pglib_opf_case2312_goc.m",
    # "pglib_opf_case2383wp_k.m",
    # "pglib_opf_case2736sp_k.m",
    # "pglib_opf_case2737sop_k.m",
    # "pglib_opf_case2742_goc.m",
    # "pglib_opf_case2746wop_k.m",
    # "pglib_opf_case2746wp_k.m",
    # "pglib_opf_case2848_rte.m",
    # "pglib_opf_case2853_sdet.m",
    # "pglib_opf_case2868_rte.m",
    # "pglib_opf_case2869_pegase.m",
    # "pglib_opf_case3012wp_k.m",
    # "pglib_opf_case3022_goc.m",
    # "pglib_opf_case3120sp_k.m",
    # "pglib_opf_case3375wp_k.m",
    # "pglib_opf_case3970_goc.m",
    # "pglib_opf_case4020_goc.m",
    # "pglib_opf_case4001_goc.m",
    # "pglib_opf_case4619_goc.m",
    # "pglib_opf_case4661_sdet.m",
    # "pglib_opf_case4837_goc.m",
    # "pglib_opf_case4917_goc.m",
    # "pglib_opf_case5658_epigrids.m",
    # "pglib_opf_case6468_rte.m",
    # "pglib_opf_case6470_rte.m",
    # "pglib_opf_case6495_rte.m",
    # "pglib_opf_case6515_rte.m",
    # "pglib_opf_case7336_epigrids.m",
    # "pglib_opf_case8387_pegase.m",
    # "pglib_opf_case9241_pegase.m",
    # "pglib_opf_case9591_goc.m",
    # "pglib_opf_case10000_goc.m",
    # "pglib_opf_case10192_epigrids.m",
    # "pglib_opf_case10480_goc.m",
    # "pglib_opf_case13659_pegase.m",
    # "pglib_opf_case14902_goc.m",
    # "pglib_opf_case20758_epigrids.m",
    # "pglib_opf_case24464_goc.m",
    # "pglib_opf_case30000_goc.m",
    # "pglib_opf_case78484_epigrids.m",
]

# Reduced list for SDP (dense SDPWRM). Dense SDP scales O(n²) variables and
# scales poorly past ~60 buses; additionally, Clarabel often hits
# NUMERICAL_ERROR on certain pathological cases (e.g. case30_ieee) where the
# reported objective becomes meaningless — those cases are skipped here.
const reduced_testlist_sdp = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case5_pjm.m",
    "pglib_opf_case14_ieee.m",
    "pglib_opf_case24_ieee_rts.m",
    "pglib_opf_case30_as.m",
]

# Reduced list for AC only
const reduced_testlist_ac = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case5_pjm.m",
    "pglib_opf_case14_ieee.m",
    "pglib_opf_case24_ieee_rts.m",
    "pglib_opf_case30_as.m",
    "pglib_opf_case39_epri.m",
    "pglib_opf_case57_ieee.m",
    "pglib_opf_case60_c.m",
    "pglib_opf_case73_ieee_rts.m",
    "pglib_opf_case118_ieee.m",
    # "pglib_opf_case162_ieee_dtc.m",
    # "pglib_opf_case179_goc.m",
    # "pglib_opf_case89_pegase.m",
    # "pglib_opf_case197_snem.m",
    # "pglib_opf_case200_activ.m",
]

# Cases with batteries (matpower_cases folder)
const battery_cases = [
    # "case5_bus_2_bat.m",
    "case14_bus_4_bat.m",
    "case57_bus_6_bat.m",
    # "case118_bus_10_bat.m",
    # "case200active_bus_20_bat.m",
]

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
function get_available_optimizer(preferred, fallback)
    try
        preferred()
    catch
        @warn "Preferred optimizer not available. Falling back to $(fallback)"
        return fallback
    end
    return preferred
end

const GLOBAL_FAILED_CASES = Dict{String, Dict{String, Vector{String}}}()

function track_failure(label::String, origin::String, case_name::String)
    if !haskey(GLOBAL_FAILED_CASES, label)
        GLOBAL_FAILED_CASES[label] = Dict("PowerModels" => String[], "Lirion" => String[])
    end
    return push!(GLOBAL_FAILED_CASES[label][origin], case_name)
end

@inline function _objective_lirion(out::Lirion.LirionResult)
    return Float64(out.objective)
end

@inline function _model_from_type(model_type::Type{<:Lirion.OPFModel})
    if model_type == Lirion.DC
        return Lirion.dc()
    elseif model_type == Lirion.AC
        return Lirion.ac()
    elseif model_type == Lirion.ACR
        return Lirion.acr()
    elseif model_type == Lirion.IVR
        return Lirion.ivr()
    elseif model_type == Lirion.LPAC
        return Lirion.lpac()
    elseif model_type == Lirion.SOCWR
        return Lirion.socwr()
    elseif model_type == Lirion.SDPWRM
        return Lirion.sdpwrm()
    elseif model_type == Lirion.SOCBF
        return Lirion.socbf()
    elseif model_type == Lirion.BFA
        return Lirion.bfa()
    elseif model_type == Lirion.NFA
        return Lirion.nfa()
    else
        error("Unsupported Lirion model type: $model_type")
    end
end

@inline function _pm_objective(file::String, model_type::Type{<:Lirion.OPFModel}, optimizer)
    pm_data = PowerModels.parse_file(file)

    has_storage = haskey(pm_data, "storage") && !isempty(pm_data["storage"])

    # NOTE:
    # - For DC:   DCPPowerModel
    # - For AC:   ACPPowerModel (via solve_ac_opf)
    # - For ACR:  ACRPowerModel
    # - For IVR:  IVRPowerModel (via solve_opf_iv)
    # - For LPAC: LPACCPowerModel
    # - For SOCWR: SOCWRConicPowerModel (no storage) / SOCWRPowerModel (storage helper)

    # IVR uses a different build method (build_opf_iv), handle separately.
    if model_type == Lirion.IVR
        return Float64(
            PowerModels.solve_opf_iv(pm_data, IVRPowerModel, optimizer)["objective"],
        )
    end

    if has_storage
        if model_type == Lirion.DC
            return Float64(
                PowerModels._solve_opf_strg(pm_data, DCPPowerModel, optimizer)["objective"],
            )
        elseif model_type == Lirion.AC
            return Float64(
                PowerModels._solve_opf_strg(pm_data, ACPPowerModel, optimizer)["objective"],
            )
        elseif model_type == Lirion.ACR
            return Float64(
                PowerModels._solve_opf_strg(pm_data, ACRPowerModel, optimizer)["objective"],
            )
        elseif model_type == Lirion.LPAC
            return Float64(
                PowerModels._solve_opf_strg(pm_data, LPACCPowerModel, optimizer)["objective"],
            )
        elseif model_type == Lirion.SOCWR
            return Float64(
                PowerModels._solve_opf_strg(pm_data, SOCWRPowerModel, optimizer)["objective"],
            )
        elseif model_type == Lirion.SOCBF
            return Float64(
                PowerModels._solve_opf_strg(pm_data, SOCBFPowerModel, optimizer)["objective"],
            )
        else
            error("Model type not recognized for storage case.")
        end
    else
        if model_type == Lirion.DC
            return Float64(
                PowerModels.solve_opf(pm_data, DCPPowerModel, optimizer)["objective"],
            )
        elseif model_type == Lirion.AC
            return Float64(PowerModels.solve_ac_opf(pm_data, optimizer)["objective"])
        elseif model_type == Lirion.ACR
            return Float64(
                PowerModels.solve_opf(pm_data, ACRPowerModel, optimizer)["objective"],
            )
        elseif model_type == Lirion.LPAC
            return Float64(
                PowerModels.solve_opf(pm_data, LPACCPowerModel, optimizer)["objective"],
            )
        elseif model_type == Lirion.SOCWR
            return Float64(
                PowerModels.solve_opf(pm_data, SOCWRConicPowerModel, optimizer)["objective"],
            )
        elseif model_type == Lirion.SDPWRM
            return Float64(
                PowerModels.solve_opf(pm_data, SDPWRMPowerModel, optimizer)["objective"],
            )
        elseif model_type == Lirion.SOCBF
            return Float64(
                PowerModels.solve_opf_bf(pm_data, SOCBFConicPowerModel, optimizer)["objective"],
            )
        elseif model_type == Lirion.NFA
            return Float64(
                PowerModels.solve_opf(pm_data, NFAPowerModel, optimizer)["objective"],
            )
        elseif model_type == Lirion.BFA
            return Float64(
                PowerModels.solve_opf_bf(pm_data, BFAPowerModel, optimizer)["objective"],
            )
        else
            error("Model type not recognized.")
        end
    end
end

function compare_with_powermodels(
    file::String,
    model_type::Type{<:Lirion.OPFModel};
    label::String,
    optimizer = Ipopt.Optimizer,
    tol = 0.05,
)
    case_name = basename(file)

    # --- Lirion (API solve) ---
    my_obj = try
        out = Lirion.solve(
            file;
            model = _model_from_type(model_type),
            algorithm = Lirion.centralized(),
            optimizer = optimizer,
            outputs = false,
            plots = false,
        )
        _objective_lirion(out)
    catch e
        @warn "Lirion failed on $case_name: $e"
        track_failure(label, "Lirion", case_name)
        return
    end

    # --- PowerModels ---
    pm_obj = try
        _pm_objective(file, model_type, optimizer)
    catch e
        @warn "PowerModels failed on $case_name: $e"
        track_failure(label, "PowerModels", case_name)
        return
    end

    if isnan(pm_obj)
        @warn "PowerModels returned NaN for $case_name"
        track_failure(label, "PowerModels", case_name)
        return
    end

    @info "File $case_name | $(model_type) | Lirion = $(round(my_obj, digits=3)) | PowerModels = $(round(pm_obj, digits=3))"
    @test isapprox(my_obj, pm_obj; rtol = tol)
end

function run_model_tests(
    model_type::Type{<:Lirion.OPFModel},
    label::String;
    optimizer = Ipopt.Optimizer,
    case_list,
    tol = 0.05,
)
    if isempty(case_list)
        @info "Empty case list for $label. No tests will be run."
        return true
    end

    test_dir = joinpath(@__DIR__, "matpower_cases")
    selected_files = [
        joinpath(test_dir, fname)
        for fname in case_list
        if isfile(joinpath(test_dir, fname))
    ]

    n = length(selected_files)
    for (i, file) in enumerate(selected_files)
        case_name = basename(file)
        @info "[$label] Rodando caso $i/$n: $case_name"
        @testset "$(label) - $(case_name)" begin
            compare_with_powermodels(
                file,
                model_type;
                label = label,
                optimizer = optimizer,
                tol = tol,
            )
        end
    end

    return true
end

function run_test_objective()
    @testset "Lirion vs PowerModels Unified Comparison" begin
        dc_optimizer = Ipopt.Optimizer
        run_model_tests(
            Lirion.DC,
            "DC";
            optimizer = dc_optimizer,
            case_list = reduced_list,
        )

        # AC: Ipopt
        run_model_tests(
            Lirion.AC,
            "AC";
            optimizer = Ipopt.Optimizer,
            case_list = reduced_testlist_ac,
        )

        run_model_tests(
            Lirion.ACR,
            "ACR";
            optimizer = Ipopt.Optimizer,
            case_list = reduced_testlist_ac,
        )

        # IVR: Ipopt (same case list as AC)
        run_model_tests(
            Lirion.IVR,
            "IVR";
            optimizer = Ipopt.Optimizer,
            case_list = reduced_testlist_ac,
        )

        # LPAC: Ipopt (same case list as AC)
        run_model_tests(
            Lirion.LPAC,
            "LPAC";
            optimizer = Ipopt.Optimizer,
            case_list = reduced_testlist_ac,
        )

        socwr_optimizer = Clarabel.Optimizer
        run_model_tests(
            Lirion.SOCWR,
            "SOCWR";
            optimizer = socwr_optimizer,
            case_list = reduced_list,
        )

        sdp_optimizer = Clarabel.Optimizer
        run_model_tests(
            Lirion.SDPWRM,
            "SDPWRM";
            optimizer = sdp_optimizer,
            case_list = reduced_testlist_sdp,
            tol = 0.01,
        )

        socbf_optimizer = Clarabel.Optimizer
        run_model_tests(
            Lirion.SOCBF,
            "SOCBF";
            optimizer = socbf_optimizer,
            case_list = reduced_list,
        )

        bfa_optimizer = Clarabel.Optimizer
        run_model_tests(
            Lirion.BFA,
            "BFA";
            optimizer = bfa_optimizer,
            case_list = reduced_list,
        )

        # Batteries
        run_model_tests(
            Lirion.DC,
            "DC com baterias";
            optimizer = dc_optimizer,
            case_list = battery_cases,
        )
        run_model_tests(
            Lirion.AC,
            "AC com baterias";
            optimizer = Ipopt.Optimizer,
            case_list = battery_cases,
        )
        run_model_tests(
            Lirion.ACR,
            "ACR com baterias";
            optimizer = Ipopt.Optimizer,
            case_list = battery_cases,
        )
        # LPAC + storage: PowerModels does not support _solve_opf_strg with
        # LPACCPowerModel, so we skip battery tests for LPAC.
        # run_model_tests(
        #     Lirion.LPAC,
        #     "LPAC com baterias";
        #     optimizer = Ipopt.Optimizer,
        #     case_list = battery_cases,
        # )
        # run_model_tests(
        #     Lirion.SOCWR,
        #     "SOCWR com baterias";
        #     optimizer = socwr_optimizer,
        #     case_list = battery_cases,
        # )
    end
end
