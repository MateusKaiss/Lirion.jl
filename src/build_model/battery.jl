function get_battery_parameters(battery_data)
    return (
        number_of_batteries = length(battery_data["storage_bus"]),
        min_storage = battery_data["energy"],
        max_storage = battery_data["energy_rating"],
        charge_ramp = battery_data["charge_rating"],
        discharge_ramp = battery_data["discharge_rating"],
        qmin = battery_data["qmin"],
        qmax = battery_data["qmax"],
        thermal_limit = battery_data["thermal_rating"],
        charge_efficiency = battery_data["charge_efficiency"],
        discharge_efficiency = battery_data["discharge_efficiency"],
    )
end

function create_battery_variables!(model, battery_data)
    params = get_battery_parameters(battery_data)
    T = get(model.ext, :time_periods, 1)

    @variable(
        model,
        storage[b = 1:params.number_of_batteries, t = 1:T],
        lower_bound = params.min_storage[b],
        upper_bound = params.max_storage[b]
    )
    @variable(
        model,
        charge[b = 1:params.number_of_batteries, t = 1:T],
        lower_bound = 0,
        upper_bound = params.charge_ramp[b]
    )
    @variable(
        model,
        qcharge[b = 1:params.number_of_batteries, t = 1:T],
        lower_bound = params.qmin[b],
        upper_bound = params.qmax[b]
    )
    @variable(
        model,
        discharge[b = 1:params.number_of_batteries, t = 1:T],
        lower_bound = 0,
        upper_bound = params.discharge_ramp[b]
    )
    @variable(
        model,
        qdischarge[b = 1:params.number_of_batteries, t = 1:T],
        lower_bound = params.qmin[b],
        upper_bound = params.qmax[b]
    )
    @variable(model, p_battery_injected[b = 1:params.number_of_batteries, t = 1:T])
    @variable(model, q_battery_injected[b = 1:params.number_of_batteries, t = 1:T])

    return nothing
end

function create_ac_dc_battery_constraints!(model, battery_data)
    params = get_battery_parameters(battery_data)
    T = get(model.ext, :time_periods, 1)

    @constraint(
        model,
        [b = 1:params.number_of_batteries, t = 1:T],
        model[:charge][b, t]^2 + model[:qcharge][b, t]^2 <= params.thermal_limit[b]^2
    )
    @constraint(
        model,
        [b = 1:params.number_of_batteries, t = 1:T],
        model[:discharge][b, t]^2 + model[:qdischarge][b, t]^2 <= params.thermal_limit[b]^2
    )
    @constraint(
        model,
        [b = 1:params.number_of_batteries, t = 1:T],
        model[:charge][b, t] * model[:discharge][b, t] == 0
    )
    @constraint(
        model,
        [b = 1:params.number_of_batteries, t = 1:T],
        model[:qcharge][b, t] * model[:qdischarge][b, t] == 0
    )
    @constraint(
        model,
        [b = 1:params.number_of_batteries, t = 1:T],
        model[:p_battery_injected][b, t] == model[:discharge][b, t] - model[:charge][b, t]
    )
    @constraint(
        model,
        [b = 1:params.number_of_batteries, t = 1:T],
        model[:q_battery_injected][b, t] ==
        model[:qdischarge][b, t] - model[:qcharge][b, t]
    )
    @constraint(
        model,
        [b = 1:params.number_of_batteries, t = 1:T],
        model[:p_battery_injected][b, t]^2 + model[:q_battery_injected][b, t]^2 <=
        params.thermal_limit[b]^2
    )

    @constraint(
        model,
        [b = 1:params.number_of_batteries],
        model[:storage][b, 1] == params.min_storage[b]
    )

    @constraint(
        model,
        [b = 1:params.number_of_batteries, t = 2:T],
        model[:storage][b, t] ==
        model[:storage][b, t-1] +
        model[:charge][b, t] * params.charge_efficiency[b] -
        model[:discharge][b, t] / params.discharge_efficiency[b]
    )

    return nothing
end

function create_socp_battery_constraints!(model, battery_data)
    params = get_battery_parameters(battery_data)
    T = get(model.ext, :time_periods, 1)

    @variable(model, p_aux[b = 1:params.number_of_batteries, t = 1:T] >= 0)
    @variable(model, q_aux[b = 1:params.number_of_batteries, t = 1:T] >= 0)
    @variable(model, z_b[b = 1:params.number_of_batteries, t = 1:T], Bin)

    @constraint(
        model,
        [b = 1:params.number_of_batteries, t = 1:T],
        [params.thermal_limit[b], p_aux[b, t], q_aux[b, t]] in JuMP.SecondOrderCone()
    )

    @constraint(
        model,
        [b = 1:params.number_of_batteries, t = 1:T],
        p_aux[b, t] >= abs(model[:p_battery_injected][b, t])
    )

    @constraint(
        model,
        [b = 1:params.number_of_batteries, t = 1:T],
        q_aux[b, t] >= abs(model[:q_battery_injected][b, t])
    )

    @constraint(
        model,
        [b = 1:params.number_of_batteries, t = 1:T],
        model[:charge][b, t] <= z_b[b, t] * params.charge_ramp[b]
    )

    @constraint(
        model,
        [b = 1:params.number_of_batteries, t = 1:T],
        model[:discharge][b, t] <= (1 - z_b[b, t]) * params.discharge_ramp[b]
    )

    @constraint(
        model,
        [b = 1:params.number_of_batteries],
        model[:storage][b, 1] == params.min_storage[b]
    )

    @constraint(
        model,
        [b = 1:params.number_of_batteries, t = 2:T],
        model[:storage][b, t] ==
        model[:storage][b, t-1] +
        model[:charge][b, t] * params.charge_efficiency[b] -
        model[:discharge][b, t] / params.discharge_efficiency[b]
    )

    @constraint(
        model,
        [b = 1:params.number_of_batteries, t = 1:T],
        model[:p_battery_injected][b, t] == model[:discharge][b, t] - model[:charge][b, t]
    )

    @constraint(
        model,
        [b = 1:params.number_of_batteries, t = 1:T],
        model[:q_battery_injected][b, t] ==
        model[:qdischarge][b, t] - model[:qcharge][b, t]
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
        create_ac_dc_battery_constraints!(model, battery_data)
    else
        @warn "No battery information found in file (DC). Skipping storage modeling."
    end
    return nothing
end

function batteries!(model::JuMP.Model, matpower_data, ::Type{AC}, ::Type{SubproblemBuild})
    if haskey(matpower_data, "storage") && !isempty(matpower_data["storage"])
        battery_data = matpower_data["storage"]
        create_battery_variables!(model, battery_data)
        create_ac_dc_battery_constraints!(model, battery_data)
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
