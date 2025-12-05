struct PowerFlowData
    bus_ids::Vector{Int}
    bus_idx_map::Dict{Int, Int}

    generator_bus_index::Vector{Int}
    battery_bus_index::Vector{Int}
    prosumer_bus_index::Vector{Int}
    demand_bus_index::Vector{Int}

    branch_bus_from::Vector{Int}
    branch_bus_to::Vector{Int}
    branch_status::Vector{Int}

    shunt_conductance::Vector{Float64}
    shunt_susceptance::Vector{Float64}

    generator_status::Vector{Int}
    is_virtual::Vector{Bool}
    has_timeseries::Bool

    # Agora vetores indexados por índice da barra (1:n_bus)
    bus_to_generators::Vector{Vector{Int}}
    bus_to_batteries::Vector{Vector{Int}}
    bus_to_prosumers::Vector{Vector{Int}}
    bus_to_branches_from::Vector{Vector{Int}}
    bus_to_branches_to::Vector{Vector{Int}}
    bus_to_inter_in_branches::Vector{Vector{Int}}
    bus_to_inter_out_branches::Vector{Vector{Int}}
end

function dict_to_busvec(
    bus_ids::AbstractVector{Int},
    bus_idx_map::Dict{Int, Int},
    mapping::Dict{Int, <:AbstractVector{Int}},
)
    nbus = length(bus_ids)
    result = [Int[] for _ in 1:nbus]

    @inbounds for (bus_id, items) in mapping
        idx = bus_idx_map[bus_id]             # mapeia bus_id -> índice interno
        result[idx] = Int.(items)             # garante Vector{Int}
    end

    return result
end

function get_power_flow_parameters(matpower_data::Dict{String, Any})
    bus_data = matpower_data["bus"]
    branch = matpower_data["branch"]
    gen = matpower_data["gen"]

    # Garante que seja Vector{Int}
    bus_ids = Int.(bus_data["bus_i"])
    bus_idx_map = get_bus_index_map(bus_ids)
    nbus = length(bus_ids)

    # Índices de barramento por tipo de elemento
    generator_bus_index = Int.(gen["gen_bus"])
    battery_bus_index = getfieldsafe(
        getfieldsafe(matpower_data, "storage", Dict()),
        "storage_bus",
        Int[],
    )
    prosumer_bus_index = getfieldsafe(
        getfieldsafe(matpower_data, "prosumer", Dict()),
        "bus",
        Int[],
    )
    demand_bus_index = bus_ids

    branch_bus_from = Int.(branch["f_bus"])
    branch_bus_to = Int.(branch["t_bus"])
    branch_status = Int.(branch["br_status"])

    shunt_conductance = Float64.(bus_data["gs"])
    shunt_susceptance = Float64.(bus_data["bs"])

    generator_status = Int.(gen["gen_status"])
    is_virtual = get(
        bus_data,
        "is_virtual",
        fill(false, nbus),
    )
    has_timeseries = get(matpower_data, "T", 1) > 1

    # -------------------------------------------
    # Construção dos mapas barra -> lista de índices
    # -------------------------------------------

    # 1) bus_to_generators
    raw_bus_to_gens = getfieldsafe(
        matpower_data,
        "bus_to_gens",
        Dict{Int, Vector{Int}}(),
    )

    bus_to_generators =
        if !isempty(raw_bus_to_gens)
            dict_to_busvec(bus_ids, bus_idx_map, raw_bus_to_gens)
        else
            # constrói a partir de generator_bus_index
            v = [Int[] for _ in 1:nbus]
            @inbounds for (g, bus_id) in enumerate(generator_bus_index)
                idx = bus_idx_map[bus_id]
                push!(v[idx], g)
            end
            v
        end

    # 2) bus_to_batteries
    raw_bus_to_bat = getfieldsafe(
        matpower_data,
        "bus_to_batteries",
        Dict{Int, Vector{Int}}(),
    )
    bus_to_batteries =
        if !isempty(raw_bus_to_bat)
            dict_to_busvec(bus_ids, bus_idx_map, raw_bus_to_bat)
        else
            v = [Int[] for _ in 1:nbus]
            @inbounds for (b, bus_id) in enumerate(battery_bus_index)
                idx = bus_idx_map[bus_id]
                push!(v[idx], b)
            end
            v
        end

    # 3) bus_to_prosumers
    raw_bus_to_pros = getfieldsafe(
        matpower_data,
        "bus_to_prosumers",
        Dict{Int, Vector{Int}}(),
    )
    bus_to_prosumers =
        if !isempty(raw_bus_to_pros)
            dict_to_busvec(bus_ids, bus_idx_map, raw_bus_to_pros)
        else
            v = [Int[] for _ in 1:nbus]
            @inbounds for (p, bus_id) in enumerate(prosumer_bus_index)
                idx = bus_idx_map[bus_id]
                push!(v[idx], p)
            end
            v
        end

    # 4) bus_to_branches_from
    raw_bus_to_br_from = getfieldsafe(
        matpower_data,
        "bus_to_branches_from",
        Dict{Int, Vector{Int}}(),
    )
    bus_to_branches_from =
        if !isempty(raw_bus_to_br_from)
            dict_to_busvec(bus_ids, bus_idx_map, raw_bus_to_br_from)
        else
            v = [Int[] for _ in 1:nbus]
            @inbounds for (l, bus_id) in enumerate(branch_bus_from)
                idx = bus_idx_map[bus_id]
                push!(v[idx], l)
            end
            v
        end

    # 5) bus_to_branches_to
    raw_bus_to_br_to = getfieldsafe(
        matpower_data,
        "bus_to_branches_to",
        Dict{Int, Vector{Int}}(),
    )
    bus_to_branches_to =
        if !isempty(raw_bus_to_br_to)
            dict_to_busvec(bus_ids, bus_idx_map, raw_bus_to_br_to)
        else
            v = [Int[] for _ in 1:nbus]
            @inbounds for (l, bus_id) in enumerate(branch_bus_to)
                idx = bus_idx_map[bus_id]
                push!(v[idx], l)
            end
            v
        end
    raw_inter_in = getfieldsafe(
        matpower_data,
        "bus_to_inter_in_branches",
        Dict{Int, Vector{Int}}(),
    )
    bus_to_inter_in_branches =
        isempty(raw_inter_in) ? [Int[] for _ in 1:nbus] :
        dict_to_busvec(bus_ids, bus_idx_map, raw_inter_in)

    raw_inter_out = getfieldsafe(
        matpower_data,
        "bus_to_inter_out_branches",
        Dict{Int, Vector{Int}}(),
    )
    bus_to_inter_out_branches =
        isempty(raw_inter_out) ? [Int[] for _ in 1:nbus] :
        dict_to_busvec(bus_ids, bus_idx_map, raw_inter_out)

    return (
        bus_ids = bus_ids,
        bus_idx_map = bus_idx_map,
        generator_bus_index = generator_bus_index,
        battery_bus_index = battery_bus_index,
        prosumer_bus_index = prosumer_bus_index,
        demand_bus_index = demand_bus_index,
        branch_bus_from = branch_bus_from,
        branch_bus_to = branch_bus_to,
        branch_status = branch_status,
        shunt_conductance = shunt_conductance,
        shunt_susceptance = shunt_susceptance,
        generator_status = generator_status,
        is_virtual = is_virtual,
        has_timeseries = has_timeseries,
        bus_to_generators = bus_to_generators,
        bus_to_batteries = bus_to_batteries,
        bus_to_prosumers = bus_to_prosumers,
        bus_to_branches_from = bus_to_branches_from,
        bus_to_branches_to = bus_to_branches_to,
        bus_to_inter_in_branches = bus_to_inter_in_branches,
        bus_to_inter_out_branches = bus_to_inter_out_branches,
    )
end

function build_power_flow_data(data::Dict{String, Any})::PowerFlowData
    p = get_power_flow_parameters(data)

    return PowerFlowData(
        p.bus_ids,
        p.bus_idx_map,
        p.generator_bus_index,
        p.battery_bus_index,
        p.prosumer_bus_index,
        p.demand_bus_index,
        p.branch_bus_from,
        p.branch_bus_to,
        p.branch_status,
        p.shunt_conductance,
        p.shunt_susceptance,
        p.generator_status,
        p.is_virtual,
        p.has_timeseries,
        p.bus_to_generators,
        p.bus_to_batteries,
        p.bus_to_prosumers,
        p.bus_to_branches_from,
        p.bus_to_branches_to,
        p.bus_to_inter_in_branches,
        p.bus_to_inter_out_branches,
    )
end
