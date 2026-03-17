# =======================================================================
# DC battery constraints
# =======================================================================

function add_dc_battery_constraints!(model::JuMP.Model, bat::BatteryData)
    _battery_no_simultaneous_charge_discharge!(model, bat)
    _battery_thermal_limit!(model, bat)
    _battery_energy_balance!(model, bat)

    T = get(model.ext, :time_periods, 1)
    B = bat.n

    @constraint(model, [b in 1:B, t in 1:T],
        model[:p_battery_injected][b, t] +
        model[:discharge][b, t] - model[:charge][b, t] ==
        bat.p_loss[b]
    )

    return nothing
end
