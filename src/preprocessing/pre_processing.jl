function load_input_data(path::String)
    result = if endswith(lowercase(path), ".m")
        open(path) do io
            return parse_matpower_string_silent(read(io, String))
        end
    elseif isdir(path)
        parse_case_from_csv_silent(path)
    else
        error("Unrecognized input format. Use .m or folder with .csv files.")
    end
    global_logger(ConsoleLogger(stderr, Logging.Info))
    return result
end

function preprocess_case_data!(data::Dict{String, Any})
    T_all = get_model_times(data)
    T = isempty(T_all) ? 1 : maximum(T_all)
    data["T"] = T

    check_time_alignment(data)

    ensure_bus_timeseries!(data, T)
    ensure_gencost_timeseries!(data, T)
    ensure_prosumer_timeseries!(data, T)

    extend_bus_timeseries!(data, T)
    extend_gencost_timeseries!(data, T)
    extend_prosumer_timeseries!(data, T)

    add_case_index_maps!(data)
    convert_to_pu!(data)
    convert_timeseries_to_pu!(data)
    return nothing
end

function parse_case_from_csv_silent(dir::String)
    with_logger(NullLogger()) do
        return parse_case_from_csv(dir)
    end
end

function parse_matpower_string_silent(data_string::String)
    with_logger(NullLogger()) do
        return _parse_matpower_string(data_string)
    end
end

function map_items_by_bus(bus_vec::Vector{T}) where {T}
    mapping = Dict{Int, Vector{Int}}()
    for (i, b) in enumerate(bus_vec)
        push!(get!(mapping, b, Int[]), i)
    end
    return mapping
end

function add_case_index_maps!(case::Dict{String, Any})
    if haskey(case, "bus")
        bus_ids = case["bus"]["bus_i"]
        case["bus_idx_map"] = Dict(b => i for (i, b) in enumerate(bus_ids))
    end

    if haskey(case, "gen")
        gen_bus = case["gen"]["gen_bus"]
        case["gen_idx_map"] = Dict(i => b for (i, b) in enumerate(gen_bus))
        case["bus_to_gens"] = map_items_by_bus(gen_bus)
    end

    if haskey(case, "branch")
        f_bus = case["branch"]["f_bus"]
        t_bus = case["branch"]["t_bus"]
        case["bus_to_branches_from"] = map_items_by_bus(f_bus)
        case["bus_to_branches_to"] = map_items_by_bus(t_bus)
        f_bus = case["branch"]["f_bus"]
        t_bus = case["branch"]["t_bus"]
        bus_to_neighbors = Dict{Int, Vector{Int}}()
        for (f, t) in zip(f_bus, t_bus)
            push!(get!(bus_to_neighbors, f, Int[]), t)
            push!(get!(bus_to_neighbors, t, Int[]), f)
        end
        case["bus_to_neighbors"] = bus_to_neighbors
    end

    if haskey(case, "prosumer")
        case["bus_to_prosumers"] = map_items_by_bus(case["prosumer"]["bus"])
    end

    if haskey(case, "storage")
        case["bus_to_batteries"] = map_items_by_bus(case["storage"]["storage_bus"])
    end

    case["bus_params"] = get_bus_parameters(case["bus"])
    case["branch_params"] = get_branch_parameters(case["branch"], case)
    case["generator_params"] = get_generator_parameters(case["gen"])
    case["power_balance_params"] = get_power_balance_parameters(case)
    case["reactive_balance_params"] = get_reactive_power_parameters(case)
    case["objective_params"] = merge(
        get_objective_parameters(
            case,
            case["gen"],
            get(case, "storage", Dict()),
            case["gencost"],
            get(case, "prosumer", Dict()),
            get(case, "gencost_timeseries", nothing),
        ),
        (; sbase = case["baseMVA"]),
    )

    return case
end

function convert_to_pu!(case::Dict{String, Any})
    Sbase = case["baseMVA"]

    for field in ["pd", "qd", "gs", "bs"]
        case["bus"][field] ./= Sbase
    end

    for field in ["pg", "qg", "pmax", "pmin", "qmax", "qmin"]
        case["gen"][field] ./= Sbase
    end

    for field in ["rate_a", "rate_b", "rate_c"]
        case["branch"][field] ./= Sbase
    end

    if haskey(case, "storage")
        for field in
            ["ps", "qs", "energy", "energy_rating", "charge_rating", "discharge_rating",
            "thermal_rating", "qmin", "qmax", "q_loss", "p_loss"]
            case["storage"][field] ./= Sbase
        end
    end

    return case
end

function convert_timeseries_to_pu!(case::Dict{String, Any})
    Sbase = case["baseMVA"]

    if haskey(case, "bus_timeseries")
        for field in ["pd", "qd", "gs", "bs"]
            if haskey(case["bus_timeseries"], field)
                case["bus_timeseries"][field] ./= Sbase
            end
        end
    end

    if haskey(case, "prosumer_timeseries")
        for field in ["p_load", "p_pv", "p_curtail"]
            if haskey(case["prosumer_timeseries"], field)
                case["prosumer_timeseries"][field] ./= Sbase
            end
        end
    end

    if haskey(case, "gencost_timeseries")
        for field in ["pmax", "pmin"]
            if haskey(case["gencost_timeseries"], field)
                case["gencost_timeseries"][field] ./= Sbase
            end
        end
    end

    return case
end
