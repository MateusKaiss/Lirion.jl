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
        case["bus_to_branches_to"]   = map_items_by_bus(t_bus)

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

    return case
end

function convert_to_pu!(case::Dict{String, Any})
    Sbase = case["baseMVA"]

    # Barras
    for field in ["pd", "qd", "gs", "bs"]
        case["bus"][field] ./= Sbase
    end

    # Geradores
    for field in ["pg", "qg", "pmax", "pmin", "qmax", "qmin"]
        case["gen"][field] ./= Sbase
    end

    # Ramos (limites térmicos)
    for field in ["rate_a", "rate_b", "rate_c"]
        case["branch"][field] ./= Sbase
    end

    # Armazenamento
    if haskey(case, "storage")
        for field in [
            "ps",
            "qs",
            "energy",
            "energy_rating",
            "charge_rating",
            "discharge_rating",
            "thermal_rating",
            "qmin",
            "qmax",
            "q_loss",
            "p_loss",
        ]
            case["storage"][field] ./= Sbase
        end
    end

    return case
end

function convert_timeseries_to_pu!(case::Dict{String, Any})
    Sbase = case["baseMVA"]

    # Timeseries de barra
    if haskey(case, "bus_timeseries")
        for field in ["pd", "qd", "gs", "bs"]
            if haskey(case["bus_timeseries"], field)
                case["bus_timeseries"][field] ./= Sbase
            end
        end
    end

    # Timeseries de prosumer
    if haskey(case, "prosumer_timeseries")
        for field in ["p_load", "p_pv", "p_curtail"]
            if haskey(case["prosumer_timeseries"], field)
                case["prosumer_timeseries"][field] ./= Sbase
            end
        end
    end

    # Timeseries de consumidor
    if haskey(case, "consumer_timeseries")
        for field in ["demand", "q_demand"]
            if haskey(case["consumer_timeseries"], field)
                case["consumer_timeseries"][field] ./= Sbase
            end
        end
    end

    return case
end

function preprocess_case_data!(
    data::Dict{String, Any},
    ::Type{T} = OPFModel,
) where {T <: OPFModel}
    # TODO Define horizonte temporal (exceto para modelos distribuídos
    # onde cada partição já vem com T bem definido).
    if T != DistributedOPFModel
        all_time_steps = get_model_times(data)
        time_horizon = isempty(all_time_steps) ? 1 : maximum(all_time_steps)
        data["T"] = time_horizon
    end

    time_horizon = data["T"]

    # Consistência básica de tempo
    check_time_alignment(data)

    # Garante tabelas timeseries mínimas
    ensure_bus_timeseries!(data, time_horizon)
    ensure_gencost_timeseries!(data, time_horizon)
    ensure_prosumer_timeseries!(data, time_horizon)
    ensure_consumer_timeseries!(data, time_horizon)

    # Estende se faltar tempo
    extend_bus_timeseries!(data, time_horizon)
    extend_gencost_timeseries!(data, time_horizon)
    extend_prosumer_timeseries!(data, time_horizon)
    extend_consumer_timeseries!(data, time_horizon)

    # Ordena por (time, id)
    sort_timeseries_by_time_and_id!(data)

    # Mapas de índice / adjacência
    add_case_index_maps!(data)

    # Converte tudo para pu
    convert_to_pu!(data)
    convert_timeseries_to_pu!(data)

    return nothing
end
getfieldsafe(data::Dict, field::String, default = missing) = get(data, field, default)

function get_bus_index_map(bus_ids::AbstractVector)
    return Dict(b => i for (i, b) in enumerate(bus_ids))
end