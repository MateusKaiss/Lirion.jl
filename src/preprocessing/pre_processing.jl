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

function preprocess_case_data!(data)
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

function add_case_index_maps!(case::Dict{String, Any})
    if haskey(case, "bus")
        bus_ids = case["bus"]["bus_i"]
        case["bus_idx_map"] = Dict(b => i for (i, b) in enumerate(bus_ids))
    end

    if haskey(case, "gen")
        gen_bus = case["gen"]["gen_bus"]
        case["gen_idx_map"] = Dict(i => b for (i, b) in enumerate(gen_bus))

        bus_to_gens = Dict{Int, Vector{Int}}()
        for (i, b) in enumerate(gen_bus)
            push!(get!(bus_to_gens, b, Int[]), i)
        end
        case["bus_to_gens"] = bus_to_gens
    end

    if haskey(case, "branch")
        f_bus = case["branch"]["f_bus"]
        t_bus = case["branch"]["t_bus"]

        bus_to_branches_from = Dict{Int, Vector{Int}}()
        bus_to_branches_to = Dict{Int, Vector{Int}}()

        for (i, f) in enumerate(f_bus)
            push!(get!(bus_to_branches_from, f, Int[]), i)
        end
        for (i, t) in enumerate(t_bus)
            push!(get!(bus_to_branches_to, t, Int[]), i)
        end

        case["bus_to_branches_from"] = bus_to_branches_from
        case["bus_to_branches_to"] = bus_to_branches_to
    end

    if haskey(case, "prosumer")
        prosumer_bus = case["prosumer"]["bus"]
        bus_to_prosumers = Dict{Int, Vector{Int}}()
        for (i, b) in enumerate(prosumer_bus)
            push!(get!(bus_to_prosumers, b, Int[]), i)
        end
        case["bus_to_prosumers"] = bus_to_prosumers
    end

    if haskey(case, "storage")
        storage_bus = case["storage"]["storage_bus"]
        bus_to_batteries = Dict{Int, Vector{Int}}()
        for (i, b) in enumerate(storage_bus)
            push!(get!(bus_to_batteries, b, Int[]), i)
        end
        case["bus_to_batteries"] = bus_to_batteries
    end

    case["bus_params"] = get_bus_parameters(case["bus"])
    case["branch_params"] = get_branch_parameters(case["branch"], case)
    case["generator_params"] = get_generator_parameters(case["gen"])
    case["power_balance_params"] = get_power_balance_parameters(case)
    case["reactive_balance_params"] = get_reactive_power_parameters(case)
    case["objective_params"] = get_objective_parameters(
        case["gen"],
        get(case, "storage", Dict()),
        case["gencost"],
        get(case, "prosumer", Dict()),
    )

    case["objective_params"] = merge(
        case["objective_params"],
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
