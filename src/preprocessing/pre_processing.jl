function load_input_data(path::String)
    if endswith(lowercase(path), ".m")
        return open(path) do io
            return parse_matpower_string_silent(read(io, String))
        end
    elseif isdir(path)
        return parse_case_from_csv_silent(path)
    else
        error("Unrecognized input format. Use .m or folder with .csv files.")
    end
end

function preprocess_case_data!(data::Dict{String, Any})
    add_case_index_maps!(data)
    return convert_to_pu!(data)
end

function parse_case_from_csv_silent(dir::String)
    original_logger = current_logger()
    logger = NullLogger()
    global_logger(logger)
    result = parse_case_from_csv(dir)
    global_logger(original_logger)
    return result
end

function parse_matpower_string_silent(data_string::String)
    original_logger = current_logger()
    global_logger(NullLogger())
    case = _parse_matpower_string(data_string)
    global_logger(original_logger)
    return case
end

function map_items_by_bus(bus_vec::Vector{Any})
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
            case["gen"],
            get(case, "storage", Dict()),
            case["gencost"],
            get(case, "prosumer", Dict()),
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
