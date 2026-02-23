# ============================================================
# ADMM zone meta (unchanged shape)
# ============================================================

struct ADMMZoneMeta
    zone_branch_idx_in_full::Vector{Int}
    full_to_zone_branch_idx::Vector{Int}
    shared_buses_zone::Vector{Int}
end

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

function _raw_snapshot(raw::CaseRaw)
    d = Dict{String, Any}()
    d["name"] = raw.name
    d["baseMVA"] = raw.baseMVA
    d["T"] = raw.T
    d["maps"] = raw.maps
    return d
end

# ============================================================
# Builder: CaseRaw -> SystemData
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
