# =======================================================================
# SOCP battery constraints
# =======================================================================
function create_socp_battery_constraints!(model::JuMP.Model, bat::BatteryData)
    T = get(model.ext, :time_periods, 1)
    B = bat.n

    @variable(model, z_b[b in 1:B, t in 1:T], Bin)

    @variable(
        model,
        ccms[b in 1:B, t in 1:T],
        lower_bound = 0.0,
        start = 0.0,
    )

    vmin_global = get(model.ext, :v_min, 0.9)
    @inbounds for b in 1:B, t in 1:T
        ub = (bat.thermal_limit[b] / vmin_global)^2
        set_upper_bound(ccms[b, t], ub)
    end

    @constraint(model, [b in 1:B, t in 1:T],
        model[:charge][b, t] <= z_b[b, t] * bat.charge_ramp[b]
    )
    @constraint(model, [b in 1:B, t in 1:T],
        model[:discharge][b, t] <= (1 - z_b[b, t]) * bat.discharge_ramp[b]
    )

    _battery_energy_balance!(model, bat)

    @constraint(model, [b in 1:B, t in 1:T],
        model[:p_battery_injected][b, t]^2 +
        model[:q_battery_injected][b, t]^2 <=
        model[:v_magnitude_sq][bat.bus[b], t] * ccms[b, t]
    )

    @constraint(model, [b in 1:B, t in 1:T],
        model[:p_battery_injected][b, t] +
        (model[:discharge][b, t] - model[:charge][b, t]) ==
        bat.p_loss[b] + bat.r[b] * ccms[b, t]
    )

    @constraint(model, [b in 1:B, t in 1:T],
        model[:q_battery_injected][b, t] ==
        model[:qsc][b, t] + bat.q_loss[b] + bat.x[b] * ccms[b, t]
    )

    return nothing
end
