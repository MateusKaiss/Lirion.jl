function _battery_energy_balance!(model::JuMP.Model, bat::BatteryData)
    T = get(model.ext, :time_periods, 1)
    B = bat.n

    # t = 1
    @constraint(model, [b in 1:B],
        model[:storage][b, 1] ==
        bat.energy_start[b] +
        model[:charge][b, 1] * bat.charge_efficiency[b] -
        model[:discharge][b, 1] / bat.discharge_efficiency[b]
    )

    # t > 1
    @constraint(model, [b in 1:B, t in 2:T],
        model[:storage][b, t] ==
        model[:storage][b, t-1] +
        model[:charge][b, t] * bat.charge_efficiency[b] -
        model[:discharge][b, t] / bat.discharge_efficiency[b]
    )

    return nothing
end

function _battery_no_simultaneous_charge_discharge!(model::JuMP.Model, bat::BatteryData)
    T = get(model.ext, :time_periods, 1)
    B = bat.n

    @constraint(model, [b in 1:B, t in 1:T],
        model[:charge][b, t] * model[:discharge][b, t] == 0
    )

    return nothing
end

function _battery_thermal_limit!(model::JuMP.Model, bat::BatteryData)
    T = get(model.ext, :time_periods, 1)
    B = bat.n

    @constraint(model, [b in 1:B, t in 1:T],
        model[:p_battery_injected][b, t]^2 +
        model[:q_battery_injected][b, t]^2 <= bat.thermal_limit[b]^2
    )

    return nothing
end
