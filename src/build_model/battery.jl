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

    @variable(model, storage[b = 1:params.number_of_batteries], lower_bound = params.min_storage[b], upper_bound = params.max_storage[b])
    @variable(model, charge[b = 1:params.number_of_batteries], lower_bound = 0, upper_bound = params.charge_ramp[b])
    @variable(model, qcharge[b = 1:params.number_of_batteries], lower_bound = params.qmin[b], upper_bound = params.qmax[b])
    @variable(model, discharge[b = 1:params.number_of_batteries], lower_bound = 0, upper_bound = params.discharge_ramp[b])
    @variable(model, qdischarge[b = 1:params.number_of_batteries], lower_bound = params.qmin[b], upper_bound = params.qmax[b])
    @variable(model, p_battery_injected[b = 1:params.number_of_batteries])
    @variable(model, q_battery_injected[b = 1:params.number_of_batteries])

    return nothing
end

function create_ac_dc_battery_constraints!(model, battery_data)
    params = get_battery_parameters(battery_data)

    @constraint(model, [b = 1:params.number_of_batteries], model[:charge][b]^2 + model[:qcharge][b]^2 <= params.thermal_limit[b]^2)
    @constraint(model, [b = 1:params.number_of_batteries], model[:discharge][b]^2 + model[:qdischarge][b]^2 <= params.thermal_limit[b]^2)
    @constraint(model, [b = 1:params.number_of_batteries], model[:charge][b] * model[:discharge][b] == 0)
    @constraint(model, [b = 1:params.number_of_batteries], model[:qcharge][b] * model[:qdischarge][b] == 0)
    @constraint(
        model,
        [b = 1:params.number_of_batteries],
        model[:storage][b] ==
        get_initial_storage(inputs, b) +
        model[:charge][b] * params.charge_efficiency[b] -
        model[:discharge][b] / params.discharge_efficiency[b]
    )
    @constraint(model, [b = 1:params.number_of_batteries], model[:p_battery_injected][b] == model[:discharge][b] - model[:charge][b])
    @constraint(model, [b = 1:params.number_of_batteries], model[:q_battery_injected][b] == model[:qdischarge][b] - model[:qcharge][b])
    @constraint(model, [b = 1:params.number_of_batteries], model[:p_battery_injected][b]^2 + model[:q_battery_injected][b]^2 <= params.thermal_limit[b]^2)

    return nothing
end

function create_socp_battery_constraints!(model, battery_data)
    params = get_battery_parameters(battery_data)

    @variable(model, p_aux[b = 1:params.number_of_batteries] >= 0)
    @variable(model, q_aux[b = 1:params.number_of_batteries] >= 0)
    @constraint(model, [b = 1:params.number_of_batteries], [params.thermal_limit[b], p_aux[b], q_aux[b]] in JuMP.SecondOrderCone())
    @constraint(model, [b = 1:params.number_of_batteries], p_aux[b] >= abs(model[:p_battery_injected][b]))
    @constraint(model, [b = 1:params.number_of_batteries], q_aux[b] >= abs(model[:q_battery_injected][b]))
    @variable(model, z_b[b = 1:params.number_of_batteries], Bin)
    @constraint(model, [b = 1:params.number_of_batteries], model[:charge][b] <= z_b[b] * params.charge_ramp[b])
    @constraint(model, [b = 1:params.number_of_batteries], model[:discharge][b] <= (1 - z_b[b]) * params.discharge_ramp[b])
    @constraint(
        model,
        [b = 1:params.number_of_batteries],
        model[:storage][b] ==
        get_initial_storage(inputs, b) +
        model[:charge][b] * params.charge_efficiency[b] -
        model[:discharge][b] / params.discharge_efficiency[b]
    )
    @constraint(model, [b = 1:params.number_of_batteries], model[:p_battery_injected][b] == model[:discharge][b] - model[:charge][b])
    @constraint(model, [b = 1:params.number_of_batteries], model[:q_battery_injected][b] == model[:qdischarge][b] - model[:qcharge][b])

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
        @warn "No battery information found in file. Skipping storage modeling."
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
