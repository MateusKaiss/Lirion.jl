@testset "Validação de objetivos para case3_timeseries.m" begin
    case_path = joinpath(@__DIR__, "matpower_cases", "case3_timeseries.m")

    @testset "DC" begin
        args = [
            case_path,
            "--model", "dc",
            "--optimizer", string(Gurobi.Optimizer),
        ]
        result_dict, _ = PowerOptimization.main(args)
        obj = result_dict["objective"]
        @info "DC objective = $obj"
        @test isapprox(obj, 17087.68770350661; atol=1e-1)
    end

    @testset "AC" begin
        args = [
            case_path,
            "--model", "ac",
            "--optimizer", string(Ipopt.Optimizer),
        ]
        result_dict, _ = PowerOptimization.main(args)
        obj = result_dict["objective"]
        @info "AC objective = $obj"
        @test isapprox(obj, 17437.9297122154; atol=1e-1)
    end

    @testset "SOCP" begin
        args = [
            case_path,
            "--model", "socp",
            "--optimizer", string(SCS.Optimizer),
        ]
        result_dict, _ = PowerOptimization.main(args)
        obj = result_dict["objective"]
        @info "SOCP objective = $obj"
        @test isapprox(obj, 17205.621331837865; atol=1e-1)
    end
end