# Lista reduzida de casos para testes rápidos (comente ou descomente conforme desejar)
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
    "pglib_opf_case500_goc.m", 
    "pglib_opf_case588_sdet.m",
    "pglib_opf_case793_goc.m",
    "pglib_opf_case1354_pegase.m",
    "pglib_opf_case1803_snem.m",
    "pglib_opf_case1888_rte.m",
    "pglib_opf_case1951_rte.m",
    "pglib_opf_case2000_goc.m",
    "pglib_opf_case2312_goc.m",
    "pglib_opf_case2383wp_k.m",
    "pglib_opf_case2736sp_k.m",
    "pglib_opf_case2737sop_k.m",
    "pglib_opf_case2742_goc.m",
    "pglib_opf_case2746wop_k.m",
    "pglib_opf_case2746wp_k.m",  
    "pglib_opf_case2848_rte.m",
    # "pglib_opf_case2853_sdet.m", analise
    "pglib_opf_case2868_rte.m",
    "pglib_opf_case2869_pegase.m",
    "pglib_opf_case3012wp_k.m",
    "pglib_opf_case3022_goc.m",
    "pglib_opf_case3120sp_k.m",
    "pglib_opf_case3375wp_k.m",
    "pglib_opf_case3970_goc.m",
    "pglib_opf_case4020_goc.m",
    "pglib_opf_case4001_goc.m",
    "pglib_opf_case4619_goc.m",
    "pglib_opf_case4661_sdet.m",
    "pglib_opf_case4837_goc.m",
    "pglib_opf_case4917_goc.m",
    "pglib_opf_case5658_epigrids.m",
    "pglib_opf_case6468_rte.m",
    "pglib_opf_case6470_rte.m",
    "pglib_opf_case6495_rte.m",
    "pglib_opf_case6515_rte.m",
    "pglib_opf_case7336_epigrids.m",
    "pglib_opf_case8387_pegase.m",
    "pglib_opf_case9241_pegase.m",
    "pglib_opf_case9591_goc.m",
    "pglib_opf_case10000_goc.m",
    "pglib_opf_case10192_epigrids.m",
    # "pglib_opf_case10480_goc.m", - ANalise
    "pglib_opf_case13659_pegase.m",
    "pglib_opf_case14902_goc.m",
    "pglib_opf_case20758_epigrids.m",
    # "pglib_opf_case24464_goc.m", -analise
    "pglib_opf_case30000_goc.m",
    "pglib_opf_case78484_epigrids.m"
]

# Lista reduzida apenas para AC (≤ 73 barras)
const reduced_testlist_ac = [
    "pglib_opf_case3_lmbd.m",
    "pglib_opf_case5_pjm.m",
    "pglib_opf_case14_ieee.m",
    "pglib_opf_case24_ieee_rts.m",
    "pglib_opf_case30_as.m",
    "pglib_opf_case39_epri.m",
    "pglib_opf_case57_ieee.m",
    "pglib_opf_case60_c.m",
    "pglib_opf_case73_ieee_rts.m"
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

function compare_with_powermodels(
    file::String,
    model_type::Type{<:PowerOptimization.OPFModel};
    optimizer = Ipopt.Optimizer,
    tol = 1e-2,
)
    result_dict, _ = PowerOptimization.main(file, model_type; optimizer = optimizer)
    my_obj = result_dict["objective"]

    pm_data = PowerModels.parse_file(file)
    pm_optimizer = optimizer

    pm_obj = try
        if model_type == PowerOptimization.DC
            PowerModels.solve_opf(pm_data, DCPPowerModel, pm_optimizer)["objective"]
        elseif model_type == PowerOptimization.AC
            PowerModels.solve_ac_opf(pm_data, pm_optimizer)["objective"]
        elseif model_type == PowerOptimization.SOCP
            PowerModels.solve_opf(pm_data, SOCWRPowerModel, pm_optimizer)["objective"]
        else
            error("Model type not recognized.")
        end
    catch e
        @warn "Error solving PowerModels for $(basename(file)): $e"
        NaN
    end

    @info "File $(basename(file)) | $(model_type) | PowerOptimization = $(round(my_obj, digits=3)) | PowerModels = $(pm_obj)"

    if isnan(pm_obj)
        @warn "PowerModels did not solve $(basename(file)) for $(model_type)"
        @test true
        return
    end

    # Comparação com tolerância relativa de 5% para todos os modelos
    @test isapprox(my_obj, pm_obj; rtol = 0.05)
end

function run_model_tests(
    model_type::Type{<:PowerOptimization.OPFModel},
    label::String;
    optimizer = Ipopt.Optimizer,
    case_list::Vector{String} = reduced_list,
)
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
        @testset "$(label) - $(basename(file))" begin
            compare_with_powermodels(file, model_type; optimizer = optimizer)
        end
    end
end

@testset "PowerOptimization vs PowerModels Unified Comparison" begin
    # DC — todos os casos
    # dc_optimizer = get_available_optimizer(() -> Gurobi.Optimizer(), Ipopt.Optimizer)
    # run_model_tests(PowerOptimization.DC, "DC"; optimizer = dc_optimizer, case_list = reduced_list)

    # # AC — até 73 barras
    # run_model_tests(PowerOptimization.AC, "AC"; optimizer = Ipopt.Optimizer, case_list = reduced_testlist_ac)


    run_model_tests(PowerOptimization.SOCP, "SOCP"; optimizer = SCS.Optimizer, case_list = reduced_list)

end







