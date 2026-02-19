# =======================================================================
# Battery variables
# =======================================================================
function create_battery_variables!(model::JuMP.Model, bat::BatteryData)
    T = get(model.ext, :time_periods, 1)
    B = bat.n

    @variable(
        model,
        storage[b in 1:B, t in 1:T],
        lower_bound = 0.0,
        upper_bound = bat.max_storage[b],
        start = bat.initial_storage[b],
    )

    @variable(
        model,
        charge[b in 1:B, t in 1:T],
        lower_bound = 0.0,
        upper_bound = bat.charge_ramp[b],
        start = 0.0,
    )

    @variable(
        model,
        discharge[b in 1:B, t in 1:T],
        lower_bound = 0.0,
        upper_bound = bat.discharge_ramp[b],
        start = 0.0,
    )

    @variable(
        model,
        p_battery_injected[b in 1:B, t in 1:T],
        lower_bound = -bat.thermal_limit[b],
        upper_bound = bat.thermal_limit[b],
        start = 0.0,
    )

    @variable(
        model,
        q_battery_injected[b in 1:B, t in 1:T],
        lower_bound = bat.qmin[b],
        upper_bound = bat.qmax[b],
        start = 0.0,
    )

    @variable(
        model,
        qsc[b in 1:B, t in 1:T],
        lower_bound = bat.qmin[b],
        upper_bound = bat.qmax[b],
        start = 0.0,
    )

    return nothing
end
