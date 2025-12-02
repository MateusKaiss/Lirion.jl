struct PowerFlowData
    bus_ids::Vector{Int}
    bus_idx_map::Dict{Int,Int}
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
    bus_to_generators::Dict{Int,Vector{Int}}
    bus_to_batteries::Dict{Int,Vector{Int}}
    bus_to_prosumers::Dict{Int,Vector{Int}}
    bus_to_branches_from::Dict{Int,Vector{Int}}
    bus_to_branches_to::Dict{Int,Vector{Int}}
    bus_to_inter_in_branches::Dict{Int,Vector{Int}}
    bus_to_inter_out_branches::Dict{Int,Vector{Int}}
end

function get_power_flow_parameters(matpower_data)
    bus_ids = matpower_data["bus"]["bus_i"]

    return (
        bus_ids         = bus_ids,
        bus_idx_map     = get_bus_index_map(bus_ids),

        generator_bus_index = matpower_data["gen"]["gen_bus"],
        battery_bus_index   = getfieldsafe(
            getfieldsafe(matpower_data, "storage", Dict()),
            "storage_bus",
            Int[],
        ),
        prosumer_bus_index  = getfieldsafe(
            getfieldsafe(matpower_data, "prosumer", Dict()),
            "bus",
            Int[],
        ),
        demand_bus_index    = bus_ids,

        branch_bus_from = matpower_data["branch"]["f_bus"],
        branch_bus_to   = matpower_data["branch"]["t_bus"],
        branch_status   = matpower_data["branch"]["br_status"],

        shunt_conductance = matpower_data["bus"]["gs"],
        shunt_susceptance = matpower_data["bus"]["bs"],

        generator_status = matpower_data["gen"]["gen_status"],
        is_virtual = get(
            matpower_data["bus"],
            "is_virtual",
            fill(false, length(matpower_data["bus"]["bus_i"])),
        ),
        has_timeseries = get(matpower_data, "T", 1) > 1,

        bus_to_generators = getfieldsafe(
            matpower_data,
            "bus_to_gens",
            Dict{Int, Vector{Int}}(),
        ),
        bus_to_batteries = getfieldsafe(
            matpower_data,
            "bus_to_batteries",
            Dict{Int, Vector{Int}}(),
        ),
        bus_to_prosumers = getfieldsafe(
            matpower_data,
            "bus_to_prosumers",
            Dict{Int, Vector{Int}}(),
        ),
        bus_to_branches_from = getfieldsafe(
            matpower_data,
            "bus_to_branches_from",
            Dict{Int, Vector{Int}}(),
        ),
        bus_to_branches_to = getfieldsafe(
            matpower_data,
            "bus_to_branches_to",
            Dict{Int, Vector{Int}}(),
        ),
        bus_to_inter_in_branches = getfieldsafe(
            matpower_data,
            "bus_to_inter_in_branches",
            Dict{Int, Vector{Int}}(),
        ),
        bus_to_inter_out_branches = getfieldsafe(
            matpower_data,
            "bus_to_inter_out_branches",
            Dict{Int, Vector{Int}}(),
        ),
    )
end

function build_power_flow_data(data::Dict{String,Any})::PowerFlowData
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