# =======================================================================
# STORAGE — energy / charge / discharge / (P, Q) variables.
# =======================================================================

function add_storage_energy_variables!(model::JuMP.Model, stor::StorageData)
    T = get(model.ext, :time_periods, 1)
    B = stor.n

    @variable(
        model,
        storage[b in 1:B, t in 1:T],
        lower_bound = 0.0,
        upper_bound = stor.max_storage[b],
        start = stor.energy_start[b],
    )

    model[:storage] = storage
    return nothing
end

function add_charge_discharge_variables!(model::JuMP.Model, stor::StorageData)
    T = get(model.ext, :time_periods, 1)
    B = stor.n

    @variable(
        model,
        charge[b in 1:B, t in 1:T],
        lower_bound = 0.0,
        upper_bound = stor.charge_ramp[b],
        start = 0.0,
    )

    @variable(
        model,
        discharge[b in 1:B, t in 1:T],
        lower_bound = 0.0,
        upper_bound = stor.discharge_ramp[b],
        start = 0.0,
    )

    model[:charge] = charge
    model[:discharge] = discharge

    return nothing
end

function add_storage_power_variables!(model::JuMP.Model, stor::StorageData)
    T = get(model.ext, :time_periods, 1)
    B = stor.n

    @variable(
        model,
        p_storage_injected[b in 1:B, t in 1:T],
        lower_bound = -stor.thermal_limit[b],
        upper_bound = stor.thermal_limit[b],
        start = 0.0,
    )

    @variable(
        model,
        q_storage_injected[b in 1:B, t in 1:T],
        lower_bound = stor.qmin[b],
        upper_bound = stor.qmax[b],
        start = 0.0,
    )

    @variable(
        model,
        qsc[b in 1:B, t in 1:T],
        lower_bound = stor.qmin[b],
        upper_bound = stor.qmax[b],
        start = 0.0,
    )

    model[:p_storage_injected] = p_storage_injected
    model[:q_storage_injected] = q_storage_injected
    model[:qsc] = qsc

    return nothing
end
