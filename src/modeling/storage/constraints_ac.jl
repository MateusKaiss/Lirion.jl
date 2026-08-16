# ============================================================
# AC storage constraints
# ============================================================

function add_ac_storage_constraints!(model::JuMP.Model, stor::StorageData)
    _storage_no_simultaneous_charge_discharge!(model, stor)
    _storage_thermal_limit!(model, stor)
    _storage_energy_balance!(model, stor)

    T = get(model.ext, :time_periods, 1)
    B = stor.n

    @constraint(model, [b in 1:B, t in 1:T],
        model[:p_storage_injected][b, t] +
        model[:discharge][b, t] - model[:charge][b, t] ==
        stor.p_loss[b] +
        stor.r[b] * (
            model[:p_storage_injected][b, t]^2 +
            model[:q_storage_injected][b, t]^2
        ) / model[:v_magnitude][stor.bus[b], t]^2
    )

    @constraint(model, [b in 1:B, t in 1:T],
        model[:q_storage_injected][b, t] ==
        model[:qsc][b, t] + stor.q_loss[b] +
        stor.x[b] * (
            model[:p_storage_injected][b, t]^2 +
            model[:q_storage_injected][b, t]^2
        ) / model[:v_magnitude][stor.bus[b], t]^2
    )

    return nothing
end
