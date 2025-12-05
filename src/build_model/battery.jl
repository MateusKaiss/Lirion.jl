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

# Balanço de energia (comum a DC/AC/SOCP)
function _battery_energy_balance!(model::JuMP.Model, bat::BatteryData)
    T = get(model.ext, :time_periods, 1)
    B = bat.n

    # t = 1
    @constraint(model, [b in 1:B],
        model[:storage][b, 1] ==
        bat.initial_storage[b] +
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

# =======================================================================
# DC battery constraints
# =======================================================================
function create_dc_battery_constraints!(model::JuMP.Model, bat::BatteryData)
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

# =======================================================================
# AC battery constraints
# =======================================================================
function create_ac_battery_constraints!(model::JuMP.Model, bat::BatteryData)
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
        ) / model[:v_magnitude][bat.bus[b], t]^2
    )

    @constraint(model, [b in 1:B, t in 1:T],
        model[:q_battery_injected][b, t] ==
        model[:qsc][b, t] + bat.q_loss[b] +
        bat.x[b] * (
            model[:p_battery_injected][b, t]^2 +
            model[:q_battery_injected][b, t]^2
        ) / model[:v_magnitude][bat.bus[b], t]^2
    )

    return nothing
end

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

    vmin_global = get(model.ext, :vmin, 0.9)
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

# =======================================================================
# Dispatch by model type (DC / AC / SOCP)
# =======================================================================
function batteries!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{DC},
    ::Type{SubproblemBuild},
)
    bat = sys.battery
    if bat.n == 0
        @info "No battery information in SystemData (DC). Skipping storage modeling."
        return nothing
    end

    create_battery_variables!(model, bat)
    create_dc_battery_constraints!(model, bat)
    return nothing
end

function batteries!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{AC},
    ::Type{SubproblemBuild},
)
    bat = sys.battery
    if bat.n == 0
        @info "No battery information in SystemData (AC). Skipping storage modeling."
        return nothing
    end

    create_battery_variables!(model, bat)
    create_ac_battery_constraints!(model, bat)
    return nothing
end

function batteries!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{SOCP},
    ::Type{SubproblemBuild},
)
    bat = sys.battery
    if bat.n == 0
        @info "No battery information in SystemData (SOCP). Skipping storage modeling."
        return nothing
    end

    create_battery_variables!(model, bat)
    create_socp_battery_constraints!(model, bat)
    return nothing
end
