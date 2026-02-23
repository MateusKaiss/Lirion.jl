# ============================================================
# SystemData
# ============================================================

struct SystemData
    baseMVA::Float64
    T::Int

    bus::BusData{Float64}
    branch::BranchData
    gen::GenData{Float64}
    battery::BatteryData
    prosumer::ProsumerData

    bus_ts::BusTimeseriesData
    prosumer_ts::ProsumerTimeseriesData
    consumer_ts::ConsumerTimeseriesData
    gencost_ts::GencostTimeseriesData

    powerflow::PowerFlowData
    objective::ObjectiveData

    raw::Dict{String, Any}
end

# ==============================================================================
# Distributed system
# ==============================================================================

struct DistributedSystem
    zones::Vector{Int}
    parts::Vector{SystemData}

    zone_branch_idx_in_full::Vector{Vector{Int}}
    full_to_zone_branch_idx::Vector{Vector{Int}}
    shared_buses_zone::Vector{Vector{Int}}

    shared_buses::Vector{Int}
    interzonal_branches::Vector{Int}
end

# ==============================================================================
# Utils
# ==============================================================================

function _raw_snapshot(raw::CaseRaw)
    d = Dict{String, Any}()
    d["name"] = raw.name
    d["baseMVA"] = raw.baseMVA
    d["T"] = raw.T
    d["maps"] = raw.maps
    return d
end

# ============================================================
# Builders
# ============================================================

function build_system_data(raw::CaseRaw)::SystemData
    bus_data = build_bus_data(raw)
    branch_data = build_branch_data(raw)
    gen_data = build_gen_data(raw)
    battery_data = build_battery_data(raw)
    prosumer_data = build_prosumer_data(raw)

    bus_ts_data = build_bus_timeseries_data(raw)
    prosumer_ts_data = build_prosumer_timeseries_data(raw)
    consumer_ts_data = build_consumer_timeseries_data(raw)
    gencost_ts_data = build_gencost_timeseries_data(raw)

    powerflow_data = build_power_flow_data(raw)
    objective_data = build_objective_data(raw)

    return SystemData(
        raw.baseMVA,
        raw.T,
        bus_data,
        branch_data,
        gen_data,
        battery_data,
        prosumer_data,
        bus_ts_data,
        prosumer_ts_data,
        consumer_ts_data,
        gencost_ts_data,
        powerflow_data,
        objective_data,
        _raw_snapshot(raw),
    )
end

function build_distributed_system(raw::CaseRaw, sectioning::Sectioning)::DistributedSystem
    parts_dict = prepare_distributed_partitions(raw, sectioning)
    interzonal_branches = _find_interzonal_branches(raw, sectioning)

    zones = sort!(collect(keys(parts_dict)))
    n_zones = length(zones)

    parts = Vector{SystemData}(undef, n_zones)
    zone_branch_idx_in_full = Vector{Vector{Int}}(undef, n_zones)
    full_to_zone_branch_idx = Vector{Vector{Int}}(undef, n_zones)
    shared_buses_zone = Vector{Vector{Int}}(undef, n_zones)

    @inbounds for i in 1:n_zones
        z = zones[i]
        zone_raw = parts_dict[z]::CaseRaw

        parts[i] = build_system_data(zone_raw)

        zone_branch_idx_in_full[i] = zone_raw.dist.zone_branch_idx_in_full
        full_to_zone_branch_idx[i] = zone_raw.dist.full_to_zone_branch_idx
        shared_buses_zone[i] = zone_raw.dist.zone_shared_bus_ids
    end

    shared_buses = raw.dist.full_shared_bus_ids

    return DistributedSystem(
        zones,
        parts,
        zone_branch_idx_in_full,
        full_to_zone_branch_idx,
        shared_buses_zone,
        shared_buses,
        interzonal_branches,
    )
end
