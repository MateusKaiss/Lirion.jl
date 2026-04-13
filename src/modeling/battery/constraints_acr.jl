# =======================================================================
# ACR battery constraints
# -----------------------------------------------------------------------
# Mirrors `add_ac_battery_constraints!` but uses rectangular voltage:
#   |V_bus|² = v_real[bus]² + v_imag[bus]²
# =======================================================================

function add_acr_battery_constraints!(model::JuMP.Model, bat::BatteryData)
    _battery_no_simultaneous_charge_discharge!(model, bat)
    _battery_thermal_limit!(model, bat)
    _battery_energy_balance!(model, bat)

    T = get(model.ext, :time_periods, 1)
    B = bat.n

    @constraint(model, [b in 1:B, t in 1:T],
        model[:p_battery_injected][b, t] +
        model[:discharge][b, t] - model[:charge][b, t] ==
        bat.p_loss[b] +
        bat.r[b] * (
            model[:p_battery_injected][b, t]^2 +
            model[:q_battery_injected][b, t]^2
        ) / (
            model[:v_real][bat.bus[b], t]^2 +
            model[:v_imag][bat.bus[b], t]^2
        )
    )

    @constraint(model, [b in 1:B, t in 1:T],
        model[:q_battery_injected][b, t] ==
        model[:qsc][b, t] + bat.q_loss[b] +
        bat.x[b] * (
            model[:p_battery_injected][b, t]^2 +
            model[:q_battery_injected][b, t]^2
        ) / (
            model[:v_real][bat.bus[b], t]^2 +
            model[:v_imag][bat.bus[b], t]^2
        )
    )

    return nothing
end
