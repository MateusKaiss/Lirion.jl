function get_objective_parameters(gen_data, battery_data)
    return (
        number_of_generators = length(gen_data["gen_bus"]),
        number_of_batteries = length(battery_data["storage_bus"]),
        generator_costs = gen_data["cost_coefficients"],
        num_polynomial_coefs = gen_data["ncost"],
        curtailment_cost = gen_data["shutdown"],
        plant_type = gen_data["apf"],
    )
end

function compute_generation_cost(model, gen_data)
    params = get_objective_parameters(gen_data, Dict())
    @expression(
        model,
        generation_cost,
        sum(
            sum(
                params.generator_costs[g][i] *
                model[:p_generation][g]^(params.num_polynomial_coefs[g] - i)
                for i in 1:params.num_polynomial_coefs[g];
                init = 0,
            )
            for g in 1:params.number_of_generators;
            init = 0,
        )
    )
end

function compute_curtailment_cost(model, gen_data)
    params = get_objective_parameters(gen_data, Dict())
    @expression(
        model,
        curtailment_cost,
        sum(
            params.curtailment_cost[g] * model[:curtailment][g]
            for g in 1:params.number_of_generators if params.plant_type[g] == "R";
            init = 0,
        )
    )
end

function compute_reactive_power_penalty(model, gen_data)
    params = get_objective_parameters(gen_data, Dict())
    @expression(
        model,
        reactive_power_penalty,
        sum(
            10e-10 * sum(
                gen_data["cost_coefficients"][g][i] *
                model[:q_generation][g]^(params.num_polynomial_coefs[g] - i)
                for i in 1:length(gen_data["cost_coefficients"][g]);
                init = 0,
            )
            for g in 1:params.number_of_generators if params.plant_type[g] == "T";
            init = 0,
        )
    )
end

function compute_battery_reactive_injection(model, battery_data)
    params = get_objective_parameters(Dict(), battery_data)
    @expression(
        model,
        battery_reactive_injection,
        -sum(
            10e-14 * model[:q_battery_injected][b]
            for b in 1:params.number_of_batteries;
            init = 0,
        )
    )
end

function objective!(model::JuMP.Model, matpower_data, battery_data, ::Type{DC})
    gen_data = matpower_data["gen"]
    battery_data = matpower_data["storage"]
    compute_generation_cost(model, gen_data)
    compute_curtailment_cost(model, gen_data)

    @objective(model, Min, model[:generation_cost] + model[:curtailment_cost])
end

function objective!(model::JuMP.Model, gen_data, battery_data, ::Type{AC})
    gen_data = matpower_data["gen"]
    battery_data = matpower_data["storage"]
    compute_generation_cost(model, gen_data)
    compute_curtailment_cost(model, gen_data)
    compute_reactive_power_penalty(model, gen_data)
    compute_battery_reactive_injection(model, battery_data)

    @objective(
        model,
        Min,
        model[:generation_cost] + model[:curtailment_cost] +
        model[:reactive_power_penalty] + model[:battery_reactive_injection]
    )
end

function objective!(model::JuMP.Model, gen_data, battery_data, ::Type{SOCP})
    gen_data = matpower_data["gen"]
    battery_data = matpower_data["storage"]
    compute_generation_cost(model, gen_data)
    compute_curtailment_cost(model, gen_data)
    compute_reactive_power_penalty(model, gen_data)
    compute_battery_reactive_injection(model, battery_data)

    @objective(
        model,
        Min,
        model[:generation_cost] + model[:curtailment_cost] +
        model[:reactive_power_penalty] + model[:battery_reactive_injection]
    )
end