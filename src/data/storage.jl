# ============================================================
# Storage data
#
# Internal representation of the `mpc.storage` table. The name "storage"
# matches both the matpower input convention and the more general
# physical model (batteries, pumped hydro, flywheels, …). The
# parametrisation here assumes a battery-like charge/discharge unit with
# efficiencies and ramp limits, but the struct is not restricted to
# electrochemical storage.
# ============================================================

struct StorageData
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

function get_storage_parameters(case_raw::CaseRaw)
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

    storage_raw = case_raw.storage::StorageRaw{Float64}

    n = length(storage_raw.storage_bus)

    return (
        n = n,
        bus = storage_raw.storage_bus,
        min_storage = zeros(Float64, n),
        max_storage = storage_raw.energy_rating,
        energy_start = storage_raw.energy,
        charge_ramp = storage_raw.charge_rating,
        discharge_ramp = storage_raw.discharge_rating,
        qmin = storage_raw.qmin,
        qmax = storage_raw.qmax,
        thermal_limit = storage_raw.thermal_rating,
        charge_efficiency = storage_raw.charge_efficiency,
        discharge_efficiency = storage_raw.discharge_efficiency,
        r = storage_raw.r,
        x = storage_raw.x,
        p_loss = storage_raw.p_loss,
        q_loss = storage_raw.q_loss,
    )
end

function build_storage_data(case_raw::CaseRaw)::StorageData
    params = get_storage_parameters(case_raw)

    return StorageData(
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
