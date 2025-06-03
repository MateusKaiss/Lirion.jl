function get_battery_parameters(battery_data)
    return (
        number_of_batteries = length(battery_data["storage_bus"]),
        min_storage = battery_data["energy"],
        max_storage = battery_data["energy_rating"],
        charge_ramp = battery_data["charge_rating"],
        discharge_ramp = battery_data["discharge_rating"],
        initial_storage = battery_data["energy"],
        qmin = battery_data["qmin"],
        qmax = battery_data["qmax"],
        thermal_limit = battery_data["thermal_rating"],
        charge_efficiency = battery_data["charge_efficiency"],
        discharge_efficiency = battery_data["discharge_efficiency"],
        r = battery_data["r"],
        x = battery_data["x"],
        p_loss = battery_data["p_loss"],
        q_loss = battery_data["q_loss"],
        bus = battery_data["storage_bus"],
    )
end

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

    # Exclusividade entre carga e descarga
    @constraint(model,
        [b = 1:params.number_of_batteries, t = 1:T],
        model[:charge][b, t] * model[:discharge][b, t] == 0
    )

    # Limite térmico (cone quadrático)
    @constraint(model,
        [b = 1:params.number_of_batteries, t = 1:T],
        model[:p_battery_injected][b, t]^2 + model[:q_battery_injected][b, t]^2 <=
        params.thermal_limit[b]^2
    )

    # Dinâmica do armazenamento (t = 1)
    @constraint(model,
        [b = 1:params.number_of_batteries],
        model[:storage][b, 1] ==
        params.initial_storage[b] +
        model[:charge][b, 1] * params.charge_efficiency[b] -
        model[:discharge][b, 1] / params.discharge_efficiency[b]
    )

    # Dinâmica do armazenamento (t > 1)
    @constraint(model,
        [b = 1:params.number_of_batteries, t = 2:T],
        model[:storage][b, t] ==
        model[:storage][b, t-1] +
        model[:charge][b, t] * params.charge_efficiency[b] -
        model[:discharge][b, t] / params.discharge_efficiency[b]
    )

    # Potência ativa com perdas (AC model)
    @constraint(model, [b = 1:params.number_of_batteries, t = 1:T],
        model[:p_battery_injected][b, t] + model[:discharge][b, t] - model[:charge][b, t] ==
        params.p_loss[b]
    )

    return nothing
end

function create_ac_battery_constraints!(model, battery_data)
    params = get_battery_parameters(battery_data)
    T = get(model.ext, :time_periods, 1)

    # Exclusividade entre carga e descarga
    @constraint(model,
        [b = 1:params.number_of_batteries, t = 1:T],
        model[:charge][b, t] * model[:discharge][b, t] == 0
    )

    # Limite térmico (cone quadrático)
    @constraint(model,
        [b = 1:params.number_of_batteries, t = 1:T],
        model[:p_battery_injected][b, t]^2 + model[:q_battery_injected][b, t]^2 <=
        params.thermal_limit[b]^2
    )

    # Dinâmica do armazenamento (t = 1)
    @constraint(model,
        [b = 1:params.number_of_batteries],
        model[:storage][b, 1] ==
        params.initial_storage[b] +
        model[:charge][b, 1] * params.charge_efficiency[b] -
        model[:discharge][b, 1] / params.discharge_efficiency[b]
    )

    # Dinâmica do armazenamento (t > 1)
    @constraint(model,
        [b = 1:params.number_of_batteries, t = 2:T],
        model[:storage][b, t] ==
        model[:storage][b, t-1] +
        model[:charge][b, t] * params.charge_efficiency[b] -
        model[:discharge][b, t] / params.discharge_efficiency[b]
    )

    # Potência ativa com perdas (AC model)
    @constraint(model,
        [b = 1:params.number_of_batteries, t = 1:T],
        model[:p_battery_injected][b, t] + model[:discharge][b, t] - model[:charge][b, t] ==
        params.p_loss[b] +
        params.r[b] * (
            model[:p_battery_injected][b, t]^2 + model[:q_battery_injected][b, t]^2
        ) / model[:v_magnitude][params.bus[b], t]^2
    )

    # Potência reativa com perdas (AC model)
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

    # Variável binária para exclusividade
    @variable(model, z_b[1:params.number_of_batteries, 1:T], Bin)

    # Variável de corrente ao quadrado (ccms)
    @variable(model, ccms[b = 1:params.number_of_batteries, t = 1:T], lower_bound = 0,
        start = 0
    )

    # Limite superior para ccms baseado no thermal_limit e vmin
    for b in 1:params.number_of_batteries
        for t in 1:T
            vmin = get(model.ext, :vmin, 0.9)
            ub = (params.thermal_limit[b] / vmin)^2
            set_upper_bound(ccms[b, t], ub)
        end
    end

    # Exclusividade entre carga e descarga
    @constraint(model, [b = 1:params.number_of_batteries, t = 1:T],
        model[:charge][b, t] <= z_b[b, t] * params.charge_ramp[b]
    )
    @constraint(model, [b = 1:params.number_of_batteries, t = 1:T],
        model[:discharge][b, t] <= (1 - z_b[b, t]) * params.discharge_ramp[b]
    )

    # Dinâmica do armazenamento
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

    # # Restrições de perdas e cone rotacionado
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
