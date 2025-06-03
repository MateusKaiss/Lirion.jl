# ---------------------------------------------------------------------
# Getting Objective Parameters (pré-processado no início)
# ---------------------------------------------------------------------
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

# ---------------------------------------------------------------------
# Cost Component Expressions
# ---------------------------------------------------------------------
function compute_generation_cost(model, params)
    T = get(model.ext, :time_periods, 1)
    active_generators = findall(params.generator_status .== 1)

    @expression(model, generation_cost,
        sum(
            sum(
                params.generator_costs[g][i] *
                (
                    params.sbase * model[:p_generation][g, t]
                )^(params.num_polynomial_coefs[g] - i)
                for i in 1:params.num_polynomial_coefs[g]
            )
            for g in active_generators, t in 1:T
        )
    )
end

function compute_curtailment_cost(model, params)
    T = get(model.ext, :time_periods, 1)
    renewable_generators = findall(
        g -> params.plant_type[g] == "R" && params.generator_status[g] == 1,
        1:params.number_of_generators,
    )

    @expression(model, curtailment_cost,
        sum(
            params.curtailment_cost[g] * model[:curtailment][g, t]
            for g in renewable_generators, t in 1:T
        )
    )
end

function compute_prosumer_cost(model, params)
    if params.number_of_prosumers == 0
        return nothing
    end
    T = get(model.ext, :time_periods, 1)
    @expression(model, prosumer_cost,
        sum(
            @inbounds params.prosumer_costs[p] * model[:p_prosumer_injected][p, t]
            for p in 1:params.number_of_prosumers, t in 1:T
        )
    )
end

function compute_reactive_power_penalty(model, params)
    T = get(model.ext, :time_periods, 1)
    thermal_generators = findall(
        g -> params.plant_type[g] == "T" && params.generator_status[g] == 1,
        1:params.number_of_generators,
    )

    @expression(model, reactive_power_penalty,
        sum(
            1e-9 * sum(
                params.generator_costs[g][i] *
                model[:q_generation][g, t]^(params.num_polynomial_coefs[g] - i)
                for i in 1:params.num_polynomial_coefs[g]
            )
            for g in thermal_generators, t in 1:T
        )
    )
end

function compute_battery_reactive_injection(model, params)
    if params.number_of_batteries == 0
        return nothing
    end
    T = get(model.ext, :time_periods, 1)
    @expression(model, battery_reactive_injection,
        -sum(
            1e-13 * model[:q_battery_injected][b, t]
            for b in 1:params.number_of_batteries, t in 1:T
        )
    )
end

# ---------------------------------------------------------------------
# Objective Assembly by Model Type (usando :objective_params)
# ---------------------------------------------------------------------
function objective!(model::JuMP.Model, matpower_data, ::Type{DC}, ::Type{SubproblemBuild})
    params = matpower_data["objective_params"]

    compute_generation_cost(model, params)

    @objective(model, Min,
        model[:generation_cost]
    )
end

function objective!(model::JuMP.Model, matpower_data, ::Type{AC}, ::Type{SubproblemBuild})
    params = matpower_data["objective_params"]

    compute_generation_cost(model, params)

    @objective(model, Min,
        model[:generation_cost]
    )
end

function objective!(model::JuMP.Model, matpower_data, ::Type{SOCP}, ::Type{SubproblemBuild})
    params = matpower_data["objective_params"]
    
    compute_generation_cost(model, params)

    @objective(model, Min,
        model[:generation_cost]
    )
end
