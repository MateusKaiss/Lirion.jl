function create_battery_variables!(model, battery_data)
    params = get_battery_parameters(battery_data)
    T = get(model.ext, :time_periods, 1)

    @variable(
        model,
        storage[b = 1:params.number_of_batteries, t = 1:T],
        lower_bound = 0,
        upper_bound = params.max_storage[b],
        start = params.initial_storage[b],
    )

    @variable(model,
        charge[b = 1:params.number_of_batteries, t = 1:T],
        lower_bound = 0,
        upper_bound = params.charge_ramp[b],
        start = 1,
    )

    @variable(model,
        discharge[b = 1:params.number_of_batteries, t = 1:T],
        lower_bound = 0,
        upper_bound = params.discharge_ramp[b],
        start = 1
    )

    @variable(model, p_battery_injected[b = 1:params.number_of_batteries, t = 1:T],
        lower_bound = -params.thermal_limit[b],
        upper_bound = params.thermal_limit[b],
        start = 0.0)

    @variable(model,
        q_battery_injected[b = 1:params.number_of_batteries, t = 1:T],
        lower_bound = params.qmin[b],
        upper_bound = params.qmax[b],
        start = 0.0,
    )

    @variable(model,
        qsc[b = 1:params.number_of_batteries, t = 1:T],
        lower_bound = params.qmin[b],
        upper_bound = params.qmax[b],
        start = 0.0,
    )

    return nothing
end

function create_dc_battery_constraints!(model, battery_data)
    params = get_battery_parameters(battery_data)
    T = get(model.ext, :time_periods, 1)

    @constraint(model,
        [b = 1:params.number_of_batteries, t = 1:T],
        model[:charge][b, t] * model[:discharge][b, t] == 0
    )

    @constraint(model,
        [b = 1:params.number_of_batteries, t = 1:T],
        model[:p_battery_injected][b, t]^2 + model[:q_battery_injected][b, t]^2 <=
        params.thermal_limit[b]^2
    )

    @constraint(model,
        [b = 1:params.number_of_batteries],
        model[:storage][b, 1] ==
        params.initial_storage[b] +
        model[:charge][b, 1] * params.charge_efficiency[b] -
        model[:discharge][b, 1] / params.discharge_efficiency[b]
    )

    @constraint(model,
        [b = 1:params.number_of_batteries, t = 2:T],
        model[:storage][b, t] ==
        model[:storage][b, t-1] +
        model[:charge][b, t] * params.charge_efficiency[b] -
        model[:discharge][b, t] / params.discharge_efficiency[b]
    )

    @constraint(model, [b = 1:params.number_of_batteries, t = 1:T],
        model[:p_battery_injected][b, t] + model[:discharge][b, t] - model[:charge][b, t] ==
        params.p_loss[b]
    )

    return nothing
end

function create_ac_battery_constraints!(model, battery_data)
    params = get_battery_parameters(battery_data)
    T = get(model.ext, :time_periods, 1)

    @constraint(model,
        [b = 1:params.number_of_batteries, t = 1:T],
        model[:charge][b, t] * model[:discharge][b, t] == 0
    )

    @constraint(model,
        [b = 1:params.number_of_batteries, t = 1:T],
        model[:p_battery_injected][b, t]^2 + model[:q_battery_injected][b, t]^2 <=
        params.thermal_limit[b]^2
    )

    @constraint(model,
        [b = 1:params.number_of_batteries],
        model[:storage][b, 1] ==
        params.initial_storage[b] +
        model[:charge][b, 1] * params.charge_efficiency[b] -
        model[:discharge][b, 1] / params.discharge_efficiency[b]
    )

    @constraint(model,
        [b = 1:params.number_of_batteries, t = 2:T],
        model[:storage][b, t] ==
        model[:storage][b, t-1] +
        model[:charge][b, t] * params.charge_efficiency[b] -
        model[:discharge][b, t] / params.discharge_efficiency[b]
    )

    @constraint(model,
        [b = 1:params.number_of_batteries, t = 1:T],
        model[:p_battery_injected][b, t] + model[:discharge][b, t] - model[:charge][b, t] ==
        params.p_loss[b] +
        params.r[b] * (
            model[:p_battery_injected][b, t]^2 + model[:q_battery_injected][b, t]^2
        ) / model[:v_magnitude][params.bus[b], t]^2
    )

    @constraint(model,
        [b = 1:params.number_of_batteries, t = 1:T],
        model[:q_battery_injected][b, t] ==
        model[:qsc][b, t] + params.q_loss[b] +
        params.x[b] * (
            model[:p_battery_injected][b, t]^2 + model[:q_battery_injected][b, t]^2
        ) / model[:v_magnitude][params.bus[b], t]^2
    )

    return nothing
end

function create_socp_battery_constraints!(model, battery_data)
    params = get_battery_parameters(battery_data)
    T = get(model.ext, :time_periods, 1)

    @variable(model, z_b[1:params.number_of_batteries, 1:T], Bin)

    @variable(model, ccms[b = 1:params.number_of_batteries, t = 1:T], lower_bound = 0,
        start = 0
    )

    for b in 1:params.number_of_batteries
        for t in 1:T
            vmin = get(model.ext, :vmin, 0.9)
            ub = (params.thermal_limit[b] / vmin)^2
            set_upper_bound(ccms[b, t], ub)
        end
    end

    @constraint(model, [b = 1:params.number_of_batteries, t = 1:T],
        model[:charge][b, t] <= z_b[b, t] * params.charge_ramp[b]
    )
    @constraint(model, [b = 1:params.number_of_batteries, t = 1:T],
        model[:discharge][b, t] <= (1 - z_b[b, t]) * params.discharge_ramp[b]
    )

    @constraint(model, [b = 1:params.number_of_batteries],
        model[:storage][b, 1] ==
        params.initial_storage[b] +
        model[:charge][b, 1] * params.charge_efficiency[b] -
        model[:discharge][b, 1] / params.discharge_efficiency[b]
    )

    @constraint(model, [b = 1:params.number_of_batteries, t = 2:T],
        model[:storage][b, t] ==
        model[:storage][b, t-1] +
        model[:charge][b, t] * params.charge_efficiency[b] -
        model[:discharge][b, t] / params.discharge_efficiency[b]
    )

    @constraint(model, [b = 1:params.number_of_batteries, t = 1:T],
        model[:p_battery_injected][b, t]^2 + model[:q_battery_injected][b, t]^2 <=
        model[:v_magnitude_sq][params.bus[b], t] * ccms[b, t]
    )

    @constraint(model, [b = 1:params.number_of_batteries, t = 1:T],
        model[:p_battery_injected][b, t] +
        (model[:discharge][b, t] - model[:charge][b, t]) ==
        params.p_loss[b] + params.r[b] * ccms[b, t]
    )

    @constraint(model, [b = 1:params.number_of_batteries, t = 1:T],
        model[:q_battery_injected][b, t] ==
        model[:qsc][b, t] + params.q_loss[b] + params.x[b] * ccms[b, t]
    )

    return nothing
end

# ------------------------------------------------------------------------
# Dispatching Battery Variables,Expressions and Constraints by Model Type
# ------------------------------------------------------------------------
function batteries!(model::JuMP.Model, matpower_data, ::Type{DC}, ::Type{SubproblemBuild})
    if haskey(matpower_data, "storage") && !isempty(matpower_data["storage"])
        battery_data = matpower_data["storage"]
        create_battery_variables!(model, battery_data)
        create_dc_battery_constraints!(model, battery_data)
    else
        @warn "No battery information found in file (DC). Skipping storage modeling."
    end
    return nothing
end

function batteries!(model::JuMP.Model, matpower_data, ::Type{AC}, ::Type{SubproblemBuild})
    if haskey(matpower_data, "storage") && !isempty(matpower_data["storage"])
        battery_data = matpower_data["storage"]
        create_battery_variables!(model, battery_data)
        create_ac_battery_constraints!(model, battery_data)
    else
        @warn "No battery information found in file (AC). Skipping storage modeling."
    end
    return nothing
end

function batteries!(model::JuMP.Model, matpower_data, ::Type{SOCP}, ::Type{SubproblemBuild})
    if haskey(matpower_data, "storage") && !isempty(matpower_data["storage"])
        battery_data = matpower_data["storage"]
        create_battery_variables!(model, battery_data)
        create_socp_battery_constraints!(model, battery_data)
    else
        @warn "No battery information found in file (SOCP). Skipping storage modeling."
    end
    return nothing
end
