# ============================================================
# Battery (Storage) data
# ============================================================

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
# Battery Parameters (CaseRaw)
# ---------------------------------------------------------------------
function get_battery_parameters(raw::CaseRaw)
    if raw.storage === nothing
        return (
            n = 0,
            bus = Int[],
            min_storage = Float64[],
            max_storage = Float64[],
            charge_ramp = Float64[],
            discharge_ramp = Float64[],
            initial_storage = Float64[],
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

    s = raw.storage::StorageRaw{Float64}
    n = length(s.storage_bus)

    return (
        n = n,
        min_storage = s.energy,
        max_storage = s.energy_rating,
        charge_ramp = s.charge_rating,
        discharge_ramp = s.discharge_rating,
        initial_storage = s.energy,
        qmin = s.qmin,
        qmax = s.qmax,
        thermal_limit = s.thermal_rating,
        charge_efficiency = s.charge_efficiency,
        discharge_efficiency = s.discharge_efficiency,
        r = s.r,
        x = s.x,
        p_loss = s.p_loss,
        q_loss = s.q_loss,
        bus = s.storage_bus,
    )
end

function build_battery_data(raw::CaseRaw)::BatteryData
    p = get_battery_parameters(raw)

    return BatteryData(
        p.n,
        p.bus,
        p.min_storage,
        p.max_storage,
        p.charge_ramp,
        p.discharge_ramp,
        p.initial_storage,
        p.qmin,
        p.qmax,
        p.thermal_limit,
        p.charge_efficiency,
        p.discharge_efficiency,
        p.r,
        p.x,
        p.p_loss,
        p.q_loss,
    )
end
