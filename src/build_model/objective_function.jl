function get_objective_parameters(
    gen_data,
    battery_data,
    gencost_data,
    prosumer_data = Dict(),
)
    return (
        number_of_generators = get(gen_data, "gen_bus", []) |> length,
        number_of_batteries = get(battery_data, "storage_bus", []) |> length,
        generator_costs = get(gencost_data, "cost_coefficients", []),
        num_polynomial_coefs = get(gencost_data, "ncost", []),
        curtailment_cost = get(gen_data, "curtailment_cost", []),
        plant_type = get(gen_data, "plant_type", []),
        generator_status = get(gen_data, "gen_status", []),
        number_of_prosumers = get(prosumer_data, "bus", []) |> length,
        prosumer_costs = get(prosumer_data, "cost_injection", []),
    )
end

function compute_generation_cost(model, gen_data, gencost_data)
    params = get_objective_parameters(gen_data, Dict(), gencost_data)
    T = get(model.ext, :time_periods, 1)

    @expression(
        model,
        generation_cost,
        sum(
            sum(
                params.generator_costs[g][i] *
                model[:p_generation][g, t]^(params.num_polynomial_coefs[g] - i)
                for i in 1:params.num_polynomial_coefs[g]
            )
            for g in 1:params.number_of_generators, t in 1:T
            if params.generator_status[g] == 1
        )
    )
end

function compute_curtailment_cost(model, gen_data)
    params = get_objective_parameters(gen_data, Dict(), Dict())
    T = get(model.ext, :time_periods, 1)

    @expression(
        model,
        curtailment_cost,
        sum(
            params.curtailment_cost[g] * model[:curtailment][g, t]
            for g in 1:params.number_of_generators, t in 1:T
            if params.plant_type[g] == "R" && params.generator_status[g] == 1
        )
    )
end

function compute_prosumer_cost(model, prosumer_data)
    params = get_objective_parameters(Dict(), Dict(), Dict(), prosumer_data)
    if params.number_of_prosumers == 0
        return nothing
    end
    T = get(model.ext, :time_periods, 1)

    @expression(
        model,
        prosumer_cost,
        sum(
            params.prosumer_costs[p] * model[:p_prosumer_injected][p, t]
            for p in 1:params.number_of_prosumers, t in 1:T
        )
    )
end

function compute_reactive_power_penalty(model, gen_data, gencost_data)
    params = get_objective_parameters(gen_data, Dict(), gencost_data)
    T = get(model.ext, :time_periods, 1)

    @expression(
        model,
        reactive_power_penalty,
        sum(
            1e-9 * sum(
                params.generator_costs[g][i] *
                model[:q_generation][g, t]^(params.num_polynomial_coefs[g] - i)
                for i in 1:length(params.generator_costs[g])
            )
            for g in 1:params.number_of_generators, t in 1:T
            if params.plant_type[g] == "T" && params.generator_status[g] == 1
        )
    )
end

function compute_battery_reactive_injection(model, battery_data)
    params = get_objective_parameters(Dict(), battery_data, Dict())
    if params.number_of_batteries == 0
        return nothing
    end
    T = get(model.ext, :time_periods, 1)

    @expression(
        model,
        battery_reactive_injection,
        -sum(
            1e-13 * model[:q_battery_injected][b, t]
            for b in 1:params.number_of_batteries, t in 1:T
        )
    )
end

function objective!(model::JuMP.Model, matpower_data, ::Type{DC}, ::Type{SubproblemBuild})
    gen_data = matpower_data["gen"]
    gencost_data = matpower_data["gencost"]
    battery_data = haskey(matpower_data, "storage") ? matpower_data["storage"] : Dict()
    prosumer_data = haskey(matpower_data, "prosumer") ? matpower_data["prosumer"] : Dict()

    compute_generation_cost(model, gen_data, gencost_data)
    compute_curtailment_cost(model, gen_data)
    compute_prosumer_cost(model, prosumer_data)

    @objective(
        model,
        Min,
        model[:generation_cost] + model[:curtailment_cost] +
        (haskey(model, :prosumer_cost) ? model[:prosumer_cost] : 0.0)
    )
end

function objective!(model::JuMP.Model, matpower_data, ::Type{AC}, ::Type{SubproblemBuild})
    gen_data = matpower_data["gen"]
    gencost_data = matpower_data["gencost"]
    battery_data = haskey(matpower_data, "storage") ? matpower_data["storage"] : Dict()
    prosumer_data = haskey(matpower_data, "prosumer") ? matpower_data["prosumer"] : Dict()

    compute_generation_cost(model, gen_data, gencost_data)
    compute_curtailment_cost(model, gen_data)
    compute_reactive_power_penalty(model, gen_data, gencost_data)
    compute_prosumer_cost(model, prosumer_data)

    if !isempty(battery_data)
        compute_battery_reactive_injection(model, battery_data)
    end

    @objective(
        model,
        Min,
        model[:generation_cost] + model[:curtailment_cost] +
        model[:reactive_power_penalty] +
        (
            haskey(model, :battery_reactive_injection) ?
            model[:battery_reactive_injection] : 0.0
        ) +
        (haskey(model, :prosumer_cost) ? model[:prosumer_cost] : 0.0)
    )
end

function objective!(model::JuMP.Model, matpower_data, ::Type{SOCP}, ::Type{SubproblemBuild})
    gen_data = matpower_data["gen"]
    gencost_data = matpower_data["gencost"]
    battery_data = haskey(matpower_data, "storage") ? matpower_data["storage"] : Dict()
    prosumer_data = haskey(matpower_data, "prosumer") ? matpower_data["prosumer"] : Dict()

    compute_generation_cost(model, gen_data, gencost_data)
    compute_curtailment_cost(model, gen_data)
    compute_reactive_power_penalty(model, gen_data, gencost_data)
    compute_prosumer_cost(model, prosumer_data)

    if !isempty(battery_data)
        compute_battery_reactive_injection(model, battery_data)
    end

    @objective(
        model,
        Min,
        model[:generation_cost] + model[:curtailment_cost] +
        # model[:reactive_power_penalty] +  # Comentado
        (
            haskey(model, :battery_reactive_injection) ?
            model[:battery_reactive_injection] : 0.0
        ) +
        (haskey(model, :prosumer_cost) ? model[:prosumer_cost] : 0.0)
    )
end