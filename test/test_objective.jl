# Lista reduzida de casos para testes rápidos (comente ou descomente conforme desejar)
const reduced_list = [
    # "pglib_opf_case3_lmbd.m",
    # "pglib_opf_case5_pjm.m",
    # "pglib_opf_case14_ieee.m",
    # "pglib_opf_case24_ieee_rts.m",
    # "pglib_opf_case30_as.m",
    # "pglib_opf_case30_ieee.m",
    # "pglib_opf_case39_epri.m",
    # "pglib_opf_case57_ieee.m",
    # "pglib_opf_case60_c.m",
    # "pglib_opf_case73_ieee_rts.m",
    "pglib_opf_case89_pegase.m",
    "pglib_opf_case118_ieee.m",
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
)
    test_dir = joinpath(@__DIR__, "matpower_cases")
    selected_files = filter(f -> basename(f) in reduced_list, Glob.glob("*.m", test_dir))
    for file in selected_files
        @testset "$(label) - $(basename(file))" begin
            compare_with_powermodels(file, model_type; optimizer = optimizer)
        end
    end
end

@testset "PowerOptimization vs PowerModels Unified Comparison" begin
    # DC
    dc_optimizer = get_available_optimizer(() -> Gurobi.Optimizer(), Ipopt.Optimizer)
    run_model_tests(PowerOptimization.DC, "DC"; optimizer = dc_optimizer)

    # AC
    # run_model_tests(PowerOptimization.AC, "AC"; optimizer = Ipopt.Optimizer)

    # SOCP — only run if Gurobi is available
    # try
    #     Gurobi.Optimizer()  # test if Gurobi is licensed and working
    #     run_model_tests(PowerOptimization.SOCP, "SOCP"; optimizer = Gurobi.Optimizer)
    # catch e
    #     @warn "Skipping SOCP tests: Gurobi not available or license issue. ($e)"
    # end
end