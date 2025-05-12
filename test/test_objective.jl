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
    # "pglib_opf_case73_ieee_rts.m",
    # "pglib_opf_case89_pegase.m",
    # "pglib_opf_case118_ieee.m",
]

function compare_with_powermodels(
    file::String,
    model_type::Type{<:PowerOptimization.OPFModel};
    optimizer = Ipopt.Optimizer,
    tol = 1e-2,
)
    result_dict, _ = PowerOptimization.main(file, model_type; optimizer=optimizer)
    @test haskey(result_dict, "objective")
    @test isfinite(result_dict["objective"])
    my_obj = result_dict["objective"]

    pm_data = PowerModels.parse_file(file)
    pm_obj = if model_type == PowerOptimization.DC
        PowerModels.solve_dc_opf(pm_data, Ipopt.Optimizer)["objective"]
    else
        PowerModels.solve_ac_opf(pm_data, Ipopt.Optimizer)["objective"]
    end

    @test isapprox(my_obj, pm_obj; atol=tol)
    @info "Arquivo $(basename(file)) | $(model_type) | PowerOptimization = $(round(my_obj, digits=3)) | PowerModels = $(round(pm_obj, digits=3))"
end

function run_model_tests(
    model_type::Type{<:PowerOptimization.OPFModel},
    label::String;
    optimizer = Ipopt.Optimizer,
)
    @testset "Comparação PowerModels x PowerOptimization ($label)" begin
        test_dir = joinpath(@__DIR__, "matpower_cases")
        all_files = Glob.glob("*.m", test_dir)
        selected_files =
            isempty(reduced_list) ? all_files :
            filter(f -> basename(f) in reduced_list, all_files)

        for file in selected_files
            @testset "$(basename(file))" begin
                compare_with_powermodels(file, model_type; optimizer=optimizer)
            end
        end
    end
end

run_model_tests(PowerOptimization.DC, "DC"; optimizer=Gurobi.Optimizer)
run_model_tests(PowerOptimization.AC, "AC"; optimizer=Ipopt.Optimizer)
# run_model_tests(PowerOptimization.SOCP, "SOCP"; optimizer=Gurobi.Optimizer)