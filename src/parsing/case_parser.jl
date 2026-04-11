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
#  CSV — implementation lives in ext/CSVExt.jl
# ============================================================

# Extension hook: CSVExt defines _csv_parse_case(dir::String)::CaseRaw
function _csv_parse_case end

function _parse_case(dir::String, ::CsvFolder)::CaseRaw
    if !applicable(_csv_parse_case, dir)
        error("""
        CSV folder parsing requires CSV.jl and DataFrames.jl.
        Add them to your project and load them first:

            using CSV, DataFrames
        """)
    end
    return _csv_parse_case(dir)
end
