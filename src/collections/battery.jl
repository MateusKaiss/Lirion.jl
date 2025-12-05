struct BatteryData
    n::Int
    bus::Vector{Int}
    min_storage::Vector{Float64}
    max_storage::Vector{Float64}
    charge_ramp::Vector{Float64}
    discharge_ramp::Vector{Float64}
    initial_storage::Vector{Float64}
    qmin::Vector{Float64}
    qmax::Vector{Float64}
    thermal_limit::Vector{Float64}
    charge_efficiency::Vector{Float64}
    discharge_efficiency::Vector{Float64}
    r::Vector{Float64}
    x::Vector{Float64}
    p_loss::Vector{Float64}
    q_loss::Vector{Float64}
end

# ---------------------------------------------------------------------
# Battery Parameters
# ---------------------------------------------------------------------
function get_battery_parameters(battery_data)
    return (
        number_of_batteries = length(battery_data["storage_bus"]),
        min_storage = battery_data["energy"],
        max_storage = battery_data["energy_rating"],
        charge_ramp = battery_data["charge_rating"],
        discharge_ramp = battery_data["discharge_rating"],
        initial_storage = battery_data["energy"],
        qmin = battery_data["qmin"],
        qmax = battery_data["qmax"],
        thermal_limit = battery_data["thermal_rating"],
        charge_efficiency = battery_data["charge_efficiency"],
        discharge_efficiency = battery_data["discharge_efficiency"],
        r = battery_data["r"],
        x = battery_data["x"],
        p_loss = battery_data["p_loss"],
        q_loss = battery_data["q_loss"],
        bus = battery_data["storage_bus"],
    )
end

function build_battery_data(data::Dict{String, Any})::BatteryData
    if !haskey(data, "storage")
        return BatteryData(
            0,
            Int[],
            Float64[], Float64[],
            Float64[], Float64[],
            Float64[],
            Float64[], Float64[],
            Float64[],
            Float64[], Float64[],
            Float64[], Float64[],
            Float64[], Float64[],
        )
    end

    params = get_battery_parameters(data["storage"])
    n = params.number_of_batteries

    return BatteryData(
        n,
        params.bus,
        params.min_storage,
        params.max_storage,
        params.charge_ramp,
        params.discharge_ramp,
        params.initial_storage,
        params.qmin,
        params.qmax,
        params.thermal_limit,
        params.charge_efficiency,
        params.discharge_efficiency,
        params.r,
        params.x,
        params.p_loss,
        params.q_loss,
    )
end
