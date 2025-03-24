function get_objective_parameters(gen_data, battery_data, gencost_data)
    return (
        number_of_generators = length(gen_data["gen_bus"]),
        number_of_batteries = haskey(battery_data, "storage_bus") ? length(battery_data["storage_bus"]) : 0,
        generator_costs = haskey(gencost_data, "cost_coefficients") ? gencost_data["cost_coefficients"] : [],
        num_polynomial_coefs = haskey(gencost_data, "ncost") ? gencost_data["ncost"] : [],
        curtailment_cost = gen_data["curtailment_cost"],
        plant_type = gen_data["plant_type"],
    )
end

function compute_generation_cost(model, gen_data, gencost_data)
    params = get_objective_parameters(gen_data, Dict(), gencost_data)
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
    params = get_objective_parameters(gen_data, Dict(), Dict())
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

function compute_reactive_power_penalty(model, gen_data, gencost_data)
    params = get_objective_parameters(gen_data, Dict(), gencost_data)
    @expression(
        model,
        reactive_power_penalty,
        sum(
            1e-9 * sum(
                params.generator_costs[g][i] *
                model[:q_generation][g]^(params.num_polynomial_coefs[g] - i)
                for i in 1:length(params.generator_costs[g]);
                init = 0,
            )
            for g in 1:params.number_of_generators if params.plant_type[g] == "T";
            init = 0,
        )
    )
end

function compute_battery_reactive_injection(model, battery_data)
    params = get_objective_parameters(Dict(), battery_data, Dict())
    if params.number_of_batteries == 0
        return nothing
    end

    @expression(
        model,
        battery_reactive_injection,
        -sum(
            1e-13 * model[:q_battery_injected][b]
            for b in 1:params.number_of_batteries;
            init = 0,
        )
    )
end

function objective!(model::JuMP.Model, matpower_data, ::Type{DC}, ::Type{SubproblemBuild})
    @show matpower_data["gencost"]["cost_coefficients"]
    gen_data = matpower_data["gen"]
    gencost_data = matpower_data["gencost"]
    battery_data = haskey(matpower_data, "storage") ? matpower_data["storage"] : Dict()

    compute_generation_cost(model, gen_data, gencost_data)
    compute_curtailment_cost(model, gen_data)

    @objective(model, Min, model[:generation_cost] + model[:curtailment_cost])
end

function objective!(model::JuMP.Model, matpower_data, ::Type{AC}, ::Type{SubproblemBuild})
    gen_data = matpower_data["gen"]
    gencost_data = matpower_data["gencost"]
    battery_data = haskey(matpower_data, "storage") ? matpower_data["storage"] : Dict()

    compute_generation_cost(model, gen_data, gencost_data)
    compute_curtailment_cost(model, gen_data)
    compute_reactive_power_penalty(model, gen_data, gencost_data)

    if !isempty(battery_data)
        compute_battery_reactive_injection(model, battery_data)
    end

    @objective(
        model,
        Min,
        model[:generation_cost] + model[:curtailment_cost] +
        model[:reactive_power_penalty] +
        (haskey(model, :battery_reactive_injection) ? model[:battery_reactive_injection] : 0)
    )
end

function objective!(model::JuMP.Model, matpower_data, ::Type{SOCP}, ::Type{SubproblemBuild})
    gen_data = matpower_data["gen"]
    gencost_data = matpower_data["gencost"]
    battery_data = haskey(matpower_data, "storage") ? matpower_data["storage"] : Dict()

    compute_generation_cost(model, gen_data, gencost_data)
    compute_curtailment_cost(model, gen_data)
    compute_reactive_power_penalty(model, gen_data, gencost_data)

    if !isempty(battery_data)
        compute_battery_reactive_injection(model, battery_data)
    end

    @objective(
        model,
        Min,
        model[:generation_cost] + model[:curtailment_cost] +
        model[:reactive_power_penalty] +
        (haskey(model, :battery_reactive_injection) ? model[:battery_reactive_injection] : 0)
    )
end