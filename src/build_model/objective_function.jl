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
