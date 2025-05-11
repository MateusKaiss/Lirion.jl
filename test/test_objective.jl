function compare_with_powermodels(
    file::String,
    model_type::Type{<:PowerOptimization.OPFModel};
    tol = 1e-2,
)
    my_result = PowerOptimization.main(file, model_type)
    @test haskey(my_result, "objective")
    @test isfinite(my_result["objective"])
    my_obj = my_result["objective"]

    pm_data = PowerModels.parse_file(file)
    pm_result = PowerModels.solve_ac_opf(pm_data, Ipopt.Optimizer)
    pm_obj = pm_result["objective"]

    @test isapprox(my_obj, pm_obj; atol = tol)
    @info "Arquivo $(basename(file)) | PowerOptimization = $(round(my_obj, digits=3)) | PowerModels = $(round(pm_obj, digits=3))"
end

const reduced_list = [
    "pglib_opf_case3_lmbd.m",
    # "pglib_opf_case5_pjm.m",
    # "pglib_opf_case14_ieee.m",
    # "pglib_opf_case24_ieee_rts.m",
    # "pglib_opf_case30_as.m",
    # "pglib_opf_case30_ieee.m",
    # "pglib_opf_case39_epri.m",
    # "pglib_opf_case57_ieee.m",
    # "pglib_opf_case60_c.m",
    # "pglib_opf_case73_ieee_rts.m",
    # "pglib_opf_case89_pegase.m",
    # "pglib_opf_case118_ieee.m",
]

@testset "Comparison with PowerModels (.m cases)" begin
    test_dir = joinpath(@__DIR__, "matpower_cases")
    all_files = Glob.glob("*.m", test_dir)
    selected_files =
        isempty(reduced_list) ? all_files :
        filter(f -> basename(f) in reduced_list, all_files)

    for file in selected_files
        @testset "$(basename(file))" begin
            compare_with_powermodels(file, PowerOptimization.AC)
        end
    end
end
