# ============================================================
# Storage constraints (shared helpers)
# ============================================================

function _storage_energy_balance!(model::JuMP.Model, stor::StorageData)
    T = get(model.ext, :time_periods, 1)
    B = stor.n

    # t = 1
    @constraint(model, [b in 1:B],
        model[:storage][b, 1] ==
        stor.energy_start[b] +
        model[:charge][b, 1] * stor.charge_efficiency[b] -
        model[:discharge][b, 1] / stor.discharge_efficiency[b]
    )

    # t > 1
    @constraint(model, [b in 1:B, t in 2:T],
        model[:storage][b, t] ==
        model[:storage][b, t-1] +
        model[:charge][b, t] * stor.charge_efficiency[b] -
        model[:discharge][b, t] / stor.discharge_efficiency[b]
    )

    return nothing
end

function _storage_no_simultaneous_charge_discharge!(
    model::JuMP.Model, stor::StorageData,
)
    T = get(model.ext, :time_periods, 1)
    B = stor.n

    @constraint(model, [b in 1:B, t in 1:T],
        model[:charge][b, t] * model[:discharge][b, t] == 0
    )

    return nothing
end

function _storage_thermal_limit!(model::JuMP.Model, stor::StorageData)
    T = get(model.ext, :time_periods, 1)
    B = stor.n

    @constraint(model, [b in 1:B, t in 1:T],
        model[:p_storage_injected][b, t]^2 +
        model[:q_storage_injected][b, t]^2 <= stor.thermal_limit[b]^2
    )

    return nothing
end
