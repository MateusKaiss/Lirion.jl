# ============================================================
# Battery (Storage) data
# ============================================================

struct BatteryData
    n::Int
    bus::Vector{Int}

    min_storage::Vector{Float64}
    max_storage::Vector{Float64}
    energy_start::Vector{Float64}

    charge_ramp::Vector{Float64}
    discharge_ramp::Vector{Float64}

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

function get_battery_parameters(case_raw::CaseRaw)
    if case_raw.storage === nothing
        return (
            n = 0,
            bus = Int[],
            min_storage = Float64[],
            max_storage = Float64[],
            energy_start = Float64[],
            charge_ramp = Float64[],
            discharge_ramp = Float64[],
            qmin = Float64[],
            qmax = Float64[],
            thermal_limit = Float64[],
            charge_efficiency = Float64[],
            discharge_efficiency = Float64[],
            r = Float64[],
            x = Float64[],
            p_loss = Float64[],
            q_loss = Float64[],
        )
    end

    battery_raw = case_raw.storage::StorageRaw{Float64}

    n = length(battery_raw.storage_bus)

    return (
        n = n,
        bus = battery_raw.storage_bus,
        min_storage = battery_raw.energy,
        max_storage = battery_raw.energy_rating,
        energy_start = battery_raw.energy,
        charge_ramp = battery_raw.charge_rating,
        discharge_ramp = battery_raw.discharge_rating,
        qmin = battery_raw.qmin,
        qmax = battery_raw.qmax,
        thermal_limit = battery_raw.thermal_rating,
        charge_efficiency = battery_raw.charge_efficiency,
        discharge_efficiency = battery_raw.discharge_efficiency,
        r = battery_raw.r,
        x = battery_raw.x,
        p_loss = battery_raw.p_loss,
        q_loss = battery_raw.q_loss,
    )
end

function build_battery_data(case_raw::CaseRaw)::BatteryData
    params = get_battery_parameters(case_raw)

    return BatteryData(
        params.n,
        params.bus,
        params.min_storage,
        params.max_storage,
        params.energy_start,
        params.charge_ramp,
        params.discharge_ramp,
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
