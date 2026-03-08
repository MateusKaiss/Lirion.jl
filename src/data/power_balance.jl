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

    bus_to_generators::Vector{Vector{Int}}
    bus_to_batteries::Vector{Vector{Int}}
    bus_to_prosumers::Vector{Vector{Int}}
    bus_to_branches_from::Vector{Vector{Int}}
    bus_to_branches_to::Vector{Vector{Int}}
    bus_to_inter_in_branches::Vector{Vector{Int}}
    bus_to_inter_out_branches::Vector{Vector{Int}}
end

function _dict_to_busvec(
    bus_ids::Vector{Int},
    bus_idx_map::Dict{Int, Int},
    mapping::Dict{Int, Vector{Int}},
)
    nbus = length(bus_ids)
    result = [Int[] for _ in 1:nbus]
    @inbounds for (bus_id, items) in mapping
        idx = bus_idx_map[Int(bus_id)]
        result[idx] = Int.(items)
    end
    return result
end

function get_power_flow_parameters(case_raw::CaseRaw)
    @assert case_raw.bus !== nothing && case_raw.branch !== nothing

    bus = case_raw.bus::BusRaw{Float64}
    br = case_raw.branch::BranchRaw{Float64}

    bus_ids = bus.bus_i
    nbus = length(bus_ids)

    # prefer maps from preprocess
    bus_idx_map = case_raw.maps.bus_idx_map

    gen_bus = case_raw.gen === nothing ? Int[] : case_raw.gen.gen_bus
    bat_bus = case_raw.storage === nothing ? Int[] : case_raw.storage.storage_bus
    pros_bus = case_raw.prosumer === nothing ? Int[] : case_raw.prosumer.bus

    branch_bus_from = br.f_bus
    branch_bus_to = br.t_bus
    branch_status = br.br_status

    shunt_gs = bus.gs
    shunt_bs = bus.bs

    gen_status = case_raw.gen === nothing ? Int[] : case_raw.gen.gen_status
    is_virtual = bus.is_virtual
    has_timeseries = case_raw.T > 1

    # maps -> busvec
    bus_to_gens =
        isempty(case_raw.maps.bus_to_gens) ? [Int[] for _ in 1:nbus] :
        _dict_to_busvec(bus_ids, bus_idx_map, case_raw.maps.bus_to_gens)

    bus_to_bats =
        isempty(case_raw.maps.bus_to_batteries) ? [Int[] for _ in 1:nbus] :
        _dict_to_busvec(bus_ids, bus_idx_map, case_raw.maps.bus_to_batteries)

    bus_to_pros =
        isempty(case_raw.maps.bus_to_prosumers) ? [Int[] for _ in 1:nbus] :
        _dict_to_busvec(bus_ids, bus_idx_map, case_raw.maps.bus_to_prosumers)

    bus_to_br_from =
        isempty(case_raw.maps.bus_to_branches_from) ? [Int[] for _ in 1:nbus] :
        _dict_to_busvec(bus_ids, bus_idx_map, case_raw.maps.bus_to_branches_from)

    bus_to_br_to =
        isempty(case_raw.maps.bus_to_branches_to) ? [Int[] for _ in 1:nbus] :
        _dict_to_busvec(bus_ids, bus_idx_map, case_raw.maps.bus_to_branches_to)

    # these may not exist in your preprocess yet
    bus_to_inter_in = [Int[] for _ in 1:nbus]
    bus_to_inter_out = [Int[] for _ in 1:nbus]

    return (
        bus_ids = bus_ids,
        bus_idx_map = bus_idx_map,
        generator_bus_index = gen_bus,
        battery_bus_index = bat_bus,
        prosumer_bus_index = pros_bus,
        demand_bus_index = bus_ids,
        branch_bus_from = branch_bus_from,
        branch_bus_to = branch_bus_to,
        branch_status = branch_status,
        shunt_conductance = shunt_gs,
        shunt_susceptance = shunt_bs,
        generator_status = gen_status,
        is_virtual = is_virtual,
        has_timeseries = has_timeseries,
        bus_to_generators = bus_to_gens,
        bus_to_batteries = bus_to_bats,
        bus_to_prosumers = bus_to_pros,
        bus_to_branches_from = bus_to_br_from,
        bus_to_branches_to = bus_to_br_to,
        bus_to_inter_in_branches = bus_to_inter_in,
        bus_to_inter_out_branches = bus_to_inter_out,
    )
end

function build_power_flow_data(case_raw::CaseRaw)::PowerFlowData
    p = get_power_flow_parameters(case_raw)
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
