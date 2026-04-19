struct PowerFlowData
    bus_ids::Vector{Int}
    bus_idx_map::Vector{Int}

    generator_bus_index::Vector{Int}
    storage_bus_index::Vector{Int}
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
    bus_to_storages::Vector{Vector{Int}}
    bus_to_branches_from::Vector{Vector{Int}}
    bus_to_branches_to::Vector{Vector{Int}}
    bus_to_inter_in_branches::Vector{Vector{Int}}
    bus_to_inter_out_branches::Vector{Vector{Int}}

    non_virtual_bus_ids::Vector{Int}

    has_bus_ts::Bool
    pd_ts::Matrix{Float64}
    qd_ts::Matrix{Float64}
end

# ------------------------------------------------------------
# Timeseries snapshot builder (bus × time dense aggregation)
# ------------------------------------------------------------

@inline function _aggregate_bus_time!(
    out::Matrix{Float64},
    bus_idx_map::Vector{Int},
    row_bus::Vector{Int},
    row_time::Vector{Int},
    row_value::Vector{Float64},
)
    _, T = size(out)
    @inbounds for i in eachindex(row_bus)
        bus_id = row_bus[i]
        t = row_time[i]
        (bus_id < 1 || bus_id > length(bus_idx_map)) && continue
        bus_idx = bus_idx_map[bus_id]
        (bus_idx == 0 || t < 1 || t > T) && continue
        out[bus_idx, t] += row_value[i]
    end
    return nothing
end

function get_power_flow_parameters(case_raw::CaseRaw)
    @assert case_raw.bus !== nothing && case_raw.branch !== nothing

    bus = case_raw.bus::BusRaw{Float64}
    br = case_raw.branch::BranchRaw{Float64}

    bus_ids = bus.bus_i
    nbus = length(bus_ids)

    bus_idx_map = case_raw.maps.bus_idx_map

    gen_bus = case_raw.gen === nothing ? Int[] : case_raw.gen.gen_bus
    stor_bus = case_raw.storage === nothing ? Int[] : case_raw.storage.storage_bus

    branch_bus_from = br.f_bus
    branch_bus_to = br.t_bus
    branch_status = br.br_status

    shunt_gs = bus.gs
    shunt_bs = bus.bs

    gen_status = case_raw.gen === nothing ? Int[] : case_raw.gen.gen_status
    is_virtual = bus.is_virtual
    has_timeseries = case_raw.T > 1

    bus_to_gens =
        isempty(case_raw.maps.bus_to_gens) ? [Int[] for _ in 1:nbus] :
        case_raw.maps.bus_to_gens

    bus_to_stors =
        isempty(case_raw.maps.bus_to_storages) ? [Int[] for _ in 1:nbus] :
        case_raw.maps.bus_to_storages

    bus_to_br_from =
        isempty(case_raw.maps.bus_to_branches_from) ? [Int[] for _ in 1:nbus] :
        case_raw.maps.bus_to_branches_from

    bus_to_br_to =
        isempty(case_raw.maps.bus_to_branches_to) ? [Int[] for _ in 1:nbus] :
        case_raw.maps.bus_to_branches_to

    bus_to_inter_in = [Int[] for _ in 1:nbus]
    bus_to_inter_out = [Int[] for _ in 1:nbus]

    non_virtual_bus_ids = Int[]
    sizehint!(non_virtual_bus_ids, nbus)
    @inbounds for i in 1:nbus
        is_virtual[i] || push!(non_virtual_bus_ids, bus_ids[i])
    end

    T = case_raw.T
    bus_ts_raw = case_raw.bus_timeseries

    has_bus_ts = bus_ts_raw !== nothing && !isempty(bus_ts_raw.time)

    T_dim = max(T, 1)
    pd_ts = zeros(Float64, nbus, T_dim)
    qd_ts = zeros(Float64, nbus, T_dim)

    if has_bus_ts
        _aggregate_bus_time!(
            pd_ts, bus_idx_map, bus_ts_raw.bus, bus_ts_raw.time, bus_ts_raw.pd,
        )
        _aggregate_bus_time!(
            qd_ts, bus_idx_map, bus_ts_raw.bus, bus_ts_raw.time, bus_ts_raw.qd,
        )
    end

    return (
        bus_ids = bus_ids,
        bus_idx_map = bus_idx_map,
        generator_bus_index = gen_bus,
        storage_bus_index = stor_bus,
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
        bus_to_storages = bus_to_stors,
        bus_to_branches_from = bus_to_br_from,
        bus_to_branches_to = bus_to_br_to,
        bus_to_inter_in_branches = bus_to_inter_in,
        bus_to_inter_out_branches = bus_to_inter_out,
        non_virtual_bus_ids = non_virtual_bus_ids,
        has_bus_ts = has_bus_ts,
        pd_ts = pd_ts,
        qd_ts = qd_ts,
    )
end

function build_power_flow_data(case_raw::CaseRaw)::PowerFlowData
    params = get_power_flow_parameters(case_raw)
    return PowerFlowData(
        params.bus_ids,
        params.bus_idx_map,
        params.generator_bus_index,
        params.storage_bus_index,
        params.demand_bus_index,
        params.branch_bus_from,
        params.branch_bus_to,
        params.branch_status,
        params.shunt_conductance,
        params.shunt_susceptance,
        params.generator_status,
        params.is_virtual,
        params.has_timeseries,
        params.bus_to_generators,
        params.bus_to_storages,
        params.bus_to_branches_from,
        params.bus_to_branches_to,
        params.bus_to_inter_in_branches,
        params.bus_to_inter_out_branches,
        params.non_virtual_bus_ids,
        params.has_bus_ts,
        params.pd_ts,
        params.qd_ts,
    )
end
