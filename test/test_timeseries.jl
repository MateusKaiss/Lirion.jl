@testset "Validação de objetivos para case3_timeseries.m" begin
    case_path = joinpath(@__DIR__, "matpower_cases", "case3_timeseries.m")

    @testset "DC" begin
        args = [
            case_path,
            "--model", "dc",
            "--optimizer", string(SCS.Optimizer),
        ]
        result_dict, _ = PowerOptimization.main(args)
        obj = result_dict["objective"]
        @info "DC objective = $obj"
        @test isapprox(obj, 17087.68770350661; atol = 1e-1)
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
        @test isapprox(obj, 17437.9297122154; atol = 1e-1)
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
        @test isapprox(obj, 17205.621331837865; atol = 1)
    end
end

@testset "Test preprocessing extension for timeseries (Dict-style)" begin
    # Caminho do caso
    case_path = joinpath(@__DIR__, "matpower_cases", "case3_timeseries_mismatch.m")

    # Carregar e pré-processar
    data = PowerOptimization.load_input_data(case_path)
    PowerOptimization.preprocess_case_data!(data)

    expected_times = [1, 2, 3]

    @testset "bus_timeseries" begin
        ts = data["bus_timeseries"]
        @test sort(unique(ts["time"])) == expected_times
        for bus in unique(ts["bus"])
            idxs = findall(x -> x == bus, ts["bus"])
            @test length(idxs) == 3
            @test ts["pd"][idxs[3]] == ts["pd"][idxs[2]]
            @test ts["qd"][idxs[3]] == ts["qd"][idxs[2]]
        end
    end

    @testset "prosumer_timeseries" begin
        ts = data["prosumer_timeseries"]
        @test sort(unique(ts["time"])) == expected_times
        for bus in unique(ts["bus"])
            idxs = findall(x -> x == bus, ts["bus"])
            @test length(idxs) == 3
            for field in ["generation", "demand", "q_generation", "q_demand"]
                @test ts[field][idxs[3]] == ts[field][idxs[2]]
            end
        end
    end

    @testset "consumer_timeseries" begin
        ts = data["consumer_timeseries"]
        @test sort(unique(ts["time"])) == expected_times
        for bus in unique(ts["bus"])
            idxs = findall(x -> x == bus, ts["bus"])
            @test length(idxs) == 3
            for field in ["demand", "q_demand"]
                @test ts[field][idxs[3]] == ts[field][idxs[2]]
            end
        end
    end

    @testset "gencost_timeseries" begin
        ts = data["gencost_timeseries"]
        @test sort(unique(ts["time"])) == expected_times
        for gen in unique(ts["gen_id"])
            idxs = findall(x -> x == gen, ts["gen_id"])
            for field in ["cost_coefficients"]
                @test ts[field][idxs[3]] == ts[field][idxs[2]]
            end
        end
    end
end
