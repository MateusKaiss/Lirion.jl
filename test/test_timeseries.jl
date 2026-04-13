function run_test_timeseries()
    @testset "Objective validation for case3_timeseries.m" begin
        case_path = joinpath(@__DIR__, "matpower_cases", "case3_timeseries.m")

        # --- DC ---
        @testset "DC" begin
            out = Lirion.solve(
                case_path;
                model = Lirion.dc(),
                algorithm = Lirion.centralized(),
                optimizer = Clarabel.Optimizer,
                outputs = false,
                plots = false,
            )
            obj = Float64(out.result["objective"])
            @info "DC objective = $obj"
            @test isapprox(obj, 17087.68770350661; atol = 1e-1)
        end

        # --- AC ---
        @testset "AC" begin
            out = Lirion.solve(
                case_path;
                model = Lirion.ac(),
                algorithm = Lirion.centralized(),
                optimizer = Ipopt.Optimizer,
                outputs = false,
                plots = false,
            )
            obj = Float64(out.result["objective"])
            @info "AC objective = $obj"
            @test isapprox(obj, 17437.9297122154; atol = 1e-1)
        end

        # --- SOCWR ---
        @testset "SOCWR" begin
            out = Lirion.solve(
                case_path;
                model = Lirion.socwr(),
                algorithm = Lirion.centralized(),
                optimizer = Clarabel.Optimizer,
                outputs = false,
                plots = false,
            )
            obj = Float64(out.result["objective"])
            @info "SOCWR objective = $obj"
            @test isapprox(obj, 17205.621331837865; atol = 1.0)
        end
    end

    # ---------------------------------------------------------------------
    # CaseRaw-style preprocessing test (updated from old Dict-style)
    # ---------------------------------------------------------------------
    @testset "Test preprocessing extension for timeseries (CaseRaw-style)" begin
        case_path = joinpath(@__DIR__, "matpower_cases", "case3_timeseries_mismatch.m")

        raw = Lirion.load_input_data(case_path)   # now returns CaseRaw
        Lirion.preprocess_case_data!(raw)

        expected_times = [1, 2, 3]

        @testset "bus_timeseries" begin
            ts = raw.bus_timeseries
            @test ts !== nothing
            @test sort(unique(ts.time)) == expected_times

            for bus in unique(ts.bus)
                idxs = findall(==(bus), ts.bus)
                @test length(idxs) == 3
                @test ts.pd[idxs[3]] == ts.pd[idxs[2]]
                @test ts.qd[idxs[3]] == ts.qd[idxs[2]]
            end
        end

        @testset "prosumer_timeseries" begin
            ts = raw.prosumer_timeseries
            @test ts !== nothing
            @test sort(unique(ts.time)) == expected_times

            for bus in unique(ts.bus)
                idxs = findall(==(bus), ts.bus)
                @test length(idxs) == 3

                # CaseRaw fields (match your ProsumerTimeseriesRaw)
                @test ts.p_pv[idxs[3]] == ts.p_pv[idxs[2]]
                @test ts.p_load[idxs[3]] == ts.p_load[idxs[2]]
                @test ts.p_curtail[idxs[3]] == ts.p_curtail[idxs[2]]
            end
        end

        @testset "consumer_timeseries" begin
            ts = raw.consumer_timeseries
            @test ts !== nothing
            @test sort(unique(ts.time)) == expected_times

            for bus in unique(ts.bus)
                idxs = findall(==(bus), ts.bus)
                @test length(idxs) == 3
                @test ts.demand[idxs[3]] == ts.demand[idxs[2]]
                @test ts.q_demand[idxs[3]] == ts.q_demand[idxs[2]]
            end
        end

        @testset "gencost_timeseries" begin
            ts = raw.gencost_timeseries
            @test ts !== nothing
            @test sort(unique(ts.time)) == expected_times

            for gen in unique(ts.gen_id)
                idxs = findall(==(gen), ts.gen_id)
                @test length(idxs) == 3
                @test ts.cost_coefficients[idxs[3]] == ts.cost_coefficients[idxs[2]]
            end
        end
    end
end
