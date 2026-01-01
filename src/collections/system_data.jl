struct ADMMZoneMeta
    branch_local_indices::Vector{Int}
    branch_global_to_local::Dict{Int, Int}
    shared_buses_local::Vector{Int}
end

struct SystemData
    baseMVA::Float64
    T::Int

    bus::BusData
    branch::BranchData
    gen::GenData
    battery::BatteryData
    prosumer::ProsumerData

    bus_ts::BusTimeseriesData
    prosumer_ts::ProsumerTimeseriesData
    consumer_ts::ConsumerTimeseriesData
    gencost_ts::GencostTimeseriesData

    powerflow::PowerFlowData
    objective::ObjectiveData

    admm::Union{Nothing, ADMMZoneMeta}

    raw::Dict{String, Any}
end

function build_system_data(
    data::Dict{String, Any};
    admm_meta::Union{Nothing, ADMMZoneMeta} = nothing,
)::SystemData
    bus_data = build_bus_data(data)
    branch_data = build_branch_data(data)
    gen_data = build_gen_data(data)
    battery_data = build_battery_data(data)
    prosumer_data = build_prosumer_data(data)

    bus_ts_data = build_bus_timeseries_data(data)
    prosumer_ts_data = build_prosumer_timeseries_data(data)
    consumer_ts_data = build_consumer_timeseries_data(data)
    gencost_ts_data = build_gencost_timeseries_data(data)

    powerflow_data = build_power_flow_data(data)
    objective_data = build_objective_data(data)

    return SystemData(
        data["baseMVA"],
        data["T"],
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
        admm_meta,
        data,
    )
end
