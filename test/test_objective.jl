# Short list of cases for rapid testing (comment or uncomment as you wish)
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
    "pglib_opf_case89_pegase.m",
    "pglib_opf_case118_ieee.m",
    "pglib_opf_case162_ieee_dtc.m",
    "pglib_opf_case179_goc.m",
    "pglib_opf_case197_snem.m",
    "pglib_opf_case200_activ.m",
    "pglib_opf_case240_pserc.m",
    "pglib_opf_case300_ieee.m",
    # "pglib_opf_case500_goc.m",
    # "pglib_opf_case588_sdet.m",
    "pglib_opf_case793_goc.m",
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
    "pglib_opf_case89_pegase.m",
    "pglib_opf_case118_ieee.m",
    "pglib_opf_case162_ieee_dtc.m",
    "pglib_opf_case179_goc.m",
    "pglib_opf_case197_snem.m",
    "pglib_opf_case200_activ.m",
]

# Cases with batteries
const battery_cases = [
    # "case5_bus_2_bat.m",
    # "case14_bus_4_bat.m",
    # "case57_bus_6_bat.m",
    # "case118_bus_10_bat.m",
    # "case200active_bus_20_bat.m",
]

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
        GLOBAL_FAILED_CASES[label] =
            Dict("PowerModels" => String[], "PowerOptimization" => String[])
    end
    return push!(GLOBAL_FAILED_CASES[label][origin], case_name)
end

function compare_with_powermodels(
    file::String,
    model_type::Type{<:PowerOptimization.OPFModel};
    label,
    optimizer = Ipopt.Optimizer,
    tol = 0.05,
)
    case_name = basename(file)

    # --- PowerOptimization ---
    my_obj = try
        args = [
            file,
            "--model", lowercase(string(model_type.name.name)),
            "--optimizer", string(optimizer),
        ]
        result_dict, _ = PowerOptimization.main(args)
        result_dict["objective"]
    catch e
        @warn "PowerOptimization failed on $case_name: $e"
        track_failure(label, "PowerOptimization", case_name)
        return
    end

    # --- PowerModels ---
    pm_obj = try
        pm_data = PowerModels.parse_file(file)
        if haskey(pm_data, "storage") && !isempty(pm_data["storage"])
            if model_type == PowerOptimization.DC
                PowerModels._solve_opf_strg(pm_data, DCPPowerModel, optimizer)["objective"]
            elseif model_type == PowerOptimization.AC
                PowerModels._solve_opf_strg(pm_data, ACPPowerModel, optimizer)["objective"]
            elseif model_type == PowerOptimization.SOCP
                PowerModels._solve_opf_strg(pm_data, SOCWRPowerModel, optimizer)["objective"]
            else
                error("Model type not recognized.")
            end
        else
            if model_type == PowerOptimization.DC
                PowerModels.solve_opf(pm_data, DCPPowerModel, optimizer)["objective"]
            elseif model_type == PowerOptimization.AC
                PowerModels.solve_ac_opf(pm_data, optimizer)["objective"]
            elseif model_type == PowerOptimization.SOCP
                PowerModels.solve_opf(pm_data, SOCWRConicPowerModel, optimizer)["objective"]
            else
                error("Model type not recognized.")
            end
        end
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

    @info "File $case_name | $(model_type) | PowerOptimization = $(round(my_obj, digits=3)) | PowerModels = $(round(pm_obj, digits=3))"
    @test isapprox(my_obj, pm_obj; rtol = tol)
end

function run_model_tests(
    model_type::Type{<:PowerOptimization.OPFModel},
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
        @info "[$label] Rodando caso $i de $n: $case_name"
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

@testset "PowerOptimization vs PowerModels Unified Comparison" begin
    dc_optimizer = get_available_optimizer(Gurobi.Optimizer, Ipopt.Optimizer)
    run_model_tests(
        PowerOptimization.DC,
        "DC";
        optimizer = dc_optimizer,
        case_list = reduced_list,
    )

    run_model_tests(
        PowerOptimization.AC,
        "AC";
        optimizer = Ipopt.Optimizer,
        case_list = reduced_testlist_ac,
    )

    socp_optimizer = get_available_optimizer(Gurobi.Optimizer, SCS.Optimizer)
    run_model_tests(
        PowerOptimization.SOCP,
        "SOCP";
        optimizer = socp_optimizer,
        case_list = reduced_list,
    )

    run_model_tests(
        PowerOptimization.DC,
        "DC com baterias";
        optimizer = dc_optimizer,
        case_list = battery_cases,
    )
    run_model_tests(
        PowerOptimization.AC,
        "AC com baterias";
        optimizer = Ipopt.Optimizer,
        case_list = battery_cases,
    )
    run_model_tests(
        PowerOptimization.SOCP,
        "SOCP com baterias";
        optimizer = socp_optimizer,
        case_list = battery_cases,
    )
end

println("\n===== GLOBAL SUMMARY OF FAILURES BY MODEL =====")
for (label, failures) in GLOBAL_FAILED_CASES
    println("\n[$label]")
    for (source, cases) in failures
        println("  $source error in $(length(cases)) case(s):")
        for c in cases
            println("    - $c")
        end
    end
end
println("==============================================\n")
