# ============================================================
#  MatPower
# ============================================================

function _parse_case(data_string::String, ::MatpowerM)::CaseRaw
    res = parse_matlab_string(data_string; meta = true)
    matlab_data = res.data
    func_name = res.function_name

    case = CaseRaw()
    case.name = func_name === nothing ? "no_name_found" : String(func_name)
    case.source_type = "matpower"
    case.source_version =
        haskey(matlab_data, "mpc.version") ? String(matlab_data["mpc.version"]) : "0.0.0+"
    case.baseMVA =
        haskey(matlab_data, "mpc.baseMVA") ? Float64(matlab_data["mpc.baseMVA"]) : 1.0

    if haskey(matlab_data, "mpc.bus")
        case.bus = _bus_from_mat(matlab_data["mpc.bus"]::Vector{Vector{String}})
    else
        @error "No bus table found in matpower file. Missing mpc.bus."
    end

    if haskey(matlab_data, "mpc.gen")
        case.gen = _gen_from_mat(matlab_data["mpc.gen"]::Vector{Vector{String}})
    else
        @error "No gen table found in matpower file. Missing mpc.gen."
    end

    if haskey(matlab_data, "mpc.gencost")
        case.gencost =
            _parse_gencost_table_raw(matlab_data["mpc.gencost"]::Vector{Vector{String}})
    else
        @error "No gencost table found in matpower file. Missing mpc.gencost."
    end

    if haskey(matlab_data, "mpc.branch")
        case.branch = _branch_from_mat(matlab_data["mpc.branch"]::Vector{Vector{String}})
    else
        @error "No branch table found in matpower file. Missing mpc.branch."
    end

    if haskey(matlab_data, "mpc.dcline")
        case.dcline = _dcline_from_mat(matlab_data["mpc.dcline"]::Vector{Vector{String}})
    end

    if haskey(matlab_data, "mpc.storage")
        case.storage = _storage_from_mat(matlab_data["mpc.storage"]::Vector{Vector{String}})
    end

    if haskey(matlab_data, "mpc.switch")
        case.switch = _switch_from_mat(matlab_data["mpc.switch"]::Vector{Vector{String}})
    end

    if haskey(matlab_data, "mpc.bus_name")
        case.bus_name =
            _busname_from_mat(matlab_data["mpc.bus_name"]::Vector{Vector{String}})
        if case.bus !== nothing && length(case.bus_name.name) != length(case.bus.bus_i)
            @error "Incorrect Matpower file: number of bus names != number of buses."
        end
    end

    if haskey(matlab_data, "mpc.prosumer")
        case.prosumer =
            _prosumer_from_mat(matlab_data["mpc.prosumer"]::Vector{Vector{String}})
    end

    if haskey(matlab_data, "mpc.consumer")
        case.consumer =
            _consumer_from_mat(matlab_data["mpc.consumer"]::Vector{Vector{String}})
    end

    if haskey(matlab_data, "mpc.bus_timeseries")
        case.bus_timeseries =
            _bus_ts_from_mat(matlab_data["mpc.bus_timeseries"]::Vector{Vector{String}})
    end

    if haskey(matlab_data, "mpc.prosumer_timeseries")
        case.prosumer_timeseries = _prosumer_ts_from_mat(
            matlab_data["mpc.prosumer_timeseries"]::Vector{Vector{String}},
        )
    end

    if haskey(matlab_data, "mpc.consumer_timeseries")
        case.consumer_timeseries = _consumer_ts_from_mat(
            matlab_data["mpc.consumer_timeseries"]::Vector{Vector{String}},
        )
    end

    if haskey(matlab_data, "mpc.gencost_timeseries")
        case.gencost_timeseries = _parse_gencost_timeseries_table_raw(
            matlab_data["mpc.gencost_timeseries"]::Vector{Vector{String}},
        )
    end

    return case
end

# ============================================================
#  CSV
# ============================================================

function _parse_case(dir::String, ::CsvFolder)::CaseRaw
    d = abspath(String(dir))

    case = CaseRaw(;
        name = splitdir(d)[end],
        source_type = "csv",
        source_version = "2",
        baseMVA = 100.0,
        T = 1,
    )

    base_path = joinpath(d, "baseMVA.csv")
    if isfile(base_path)
        df = _read_csv(base_path)
        if nrow(df) > 0 && _hascol(df, "baseMVA")
            x = df[1, _sym("baseMVA")]
            case.baseMVA = (x === missing) ? 100.0 : Float64(x)
        else
            @warn "baseMVA.csv empty or missing baseMVA. Using 100.0"
        end
    else
        @warn "baseMVA.csv not found. Using 100.0"
    end

    # -------------------------
    # core tables
    # -------------------------

    let p = joinpath(d, "bus.csv")
        if isfile(p)
            tbl = _parse_csv_table_dict(p, _mp_bus_columns_csv)
            case.bus = _bus_from_tbl(tbl)
        end
    end

    let p = joinpath(d, "gen.csv")
        if isfile(p)
            tbl = _parse_csv_table_dict(p, _mp_gen_columns_csv)
            n = length(tbl["gen_bus"])

            if !haskey(tbl, "plant_type")
                tbl["plant_type"] = fill("T", n)
            else
                pt = tbl["plant_type"]
                if !(pt isa Vector{String})
                    tbl["plant_type"] = String.(pt)
                end
                @inbounds for i in 1:n
                    if isempty(tbl["plant_type"][i])
                        tbl["plant_type"][i] = "T"
                    end
                end
            end

            if !haskey(tbl, "curtailment_cost")
                tbl["curtailment_cost"] = fill(0.0, n)
            end

            case.gen = _gen_from_tbl(tbl)
        end
    end

    let p = joinpath(d, "branch.csv")
        if isfile(p)
            tbl = _parse_csv_table_dict(p, _mp_branch_columns_csv)
            # NÃO converte angmin/angmax aqui.
            case.branch = _branch_from_tbl(tbl)
        end
    end

    let p = joinpath(d, "dcline.csv")
        if isfile(p)
            tbl = _parse_csv_table_dict(p, _mp_dcline_columns_csv)
            case.dcline = _dcline_from_tbl(tbl)
        end
    end

    let p = joinpath(d, "storage.csv")
        if isfile(p)
            tbl = _parse_csv_table_dict(p, _mp_storage_columns_csv)
            case.storage = _storage_from_tbl(tbl)
        end
    end

    let p = joinpath(d, "switch.csv")
        if isfile(p)
            tbl = _parse_csv_table_dict(p, _mp_switch_columns_csv)
            case.switch = _switch_from_tbl(tbl)
        end
    end

    let p = joinpath(d, "prosumer.csv")
        if isfile(p)
            tbl = _parse_csv_table_dict(p, _mp_prosumer_columns_csv)
            case.prosumer = _prosumer_from_tbl(tbl)
        end
    end

    let p = joinpath(d, "consumer.csv")
        if isfile(p)
            tbl = _parse_csv_table_dict(p, _mp_consumer_columns_csv)
            case.consumer = _consumer_from_tbl(tbl)
        end
    end

    let p = joinpath(d, "bus_name.csv")
        if isfile(p)
            tbl = _parse_csv_table_dict(p, _mp_bus_name_columns_csv)
            case.bus_name = _busname_from_tbl(tbl)
            if case.bus !== nothing && length(case.bus_name.name) != length(case.bus.bus_i)
                @error "Incorrect CSV case: number of bus names != number of buses."
            end
        end
    end

    let p = joinpath(d, "gencost.csv")
        if isfile(p)
            case.gencost = _parse_gencost_csv(p)
        end
    end

    # -------------------------
    # timeseries folder
    # -------------------------
    tdir = joinpath(d, "time")
    if isdir(tdir)
        let p = joinpath(tdir, "bus_timeseries.csv")
            if isfile(p)
                tbl = _parse_csv_table_dict(p, _mp_bus_timeseries_columns_csv)
                case.bus_timeseries = _bus_ts_from_tbl(tbl)
                if !isempty(case.bus_timeseries.time)
                    case.T = maximum(case.bus_timeseries.time)
                end
            end
        end

        let p = joinpath(tdir, "consumer_timeseries.csv")
            if isfile(p)
                tbl = _parse_csv_table_dict(p, _mp_consumer_timeseries_columns_csv)
                case.consumer_timeseries = _consumer_ts_from_tbl(tbl)
                if !isempty(case.consumer_timeseries.time)
                    case.T = max(case.T, maximum(case.consumer_timeseries.time))
                end
            end
        end

        let p = joinpath(tdir, "prosumer_timeseries.csv")
            if isfile(p)
                df = _read_csv(p)
                tbl = Dict{String, Any}()

                tbl["time"] = _col_or_default(df, "time", Int, 1)
                tbl["bus"] = _col_or_default(df, "bus", Int, 1)

                if _hascol(df, "p_pv") || _hascol(df, "p_load") || _hascol(df, "p_curtail")
                    tbl["p_pv"] = _col_or_default(df, "p_pv", Float64, 0.0)
                    tbl["p_load"] = _col_or_default(df, "p_load", Float64, 0.0)
                    tbl["p_curtail"] = _col_or_default(df, "p_curtail", Float64, 0.0)
                else
                    tbl["generation"] = _col_or_default(df, "generation", Float64, 0.0)
                    tbl["demand"] = _col_or_default(df, "demand", Float64, 0.0)
                end

                case.prosumer_timeseries = _prosumer_ts_from_tbl_csv(tbl)
                if !isempty(case.prosumer_timeseries.time)
                    case.T = max(case.T, maximum(case.prosumer_timeseries.time))
                end
            end
        end

        let p = joinpath(tdir, "gencost_timeseries.csv")
            if isfile(p)
                case.gencost_timeseries = _parse_gencost_timeseries_csv(p)
                if !isempty(case.gencost_timeseries.time)
                    case.T = max(case.T, maximum(case.gencost_timeseries.time))
                end
            end
        end
    end

    return case
end
