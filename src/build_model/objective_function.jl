function compute_generation_cost(model, params)
    T = get(model.ext, :time_periods, 1)
    active_generators = findall(params.generator_status .== 1)

    @expression(model, generation_cost,
        sum(
            begin
                coefs = params.generator_costs[t][g]
                sum(
                    coefs[i] *
                    (params.sbase * model[:p_generation][g, t])^(length(coefs) - i)
                    for i in eachindex(coefs)
                )
            end
            for t in 1:T, g in active_generators
            if haskey(params.generator_costs, t) && haskey(params.generator_costs[t], g)
        )
    )
end

# function compute_curtailment_cost(model, params)
#     T = get(model.ext, :time_periods, 1)
#     renewable_generators = findall(
#         g -> params.plant_type[g] == "R" && params.generator_status[g] == 1,
#         1:params.number_of_generators,
#     )

#     @expression(model, curtailment_cost,
#         sum(
#             params.curtailment_cost[g] * model[:curtailment][g, t]
#             for g in renewable_generators, t in 1:T
#         )
#     )
# end

# ---------------------------------------------------------------------
# Objective Assembly by Model Type (usando :objective_params)
# ---------------------------------------------------------------------
function objective!(model::JuMP.Model, matpower_data, ::Type{DC}, ::Type{SubproblemBuild})
    params = matpower_data["objective_params"]
    compute_generation_cost(model, params)

    T = get(model.ext, :time_periods, 1)
    virtual_cost_coeff = get(matpower_data, "virtual_cost_coeff", 750.0)
    penalty_expr = 0.0

    if haskey(model.ext, :p_virtual)
        p_virtual = model.ext[:p_virtual]
        if isa(virtual_cost_coeff, Number)
            @expression(model, virtual_penalty,
                sum(
                    virtual_cost_coeff * p_virtual[b, t]^2 for b in axes(p_virtual, 1),
                    t in 1:T
                )
            )
        elseif isa(virtual_cost_coeff, Dict)
            @expression(model, virtual_penalty,
                sum(
                    virtual_cost_coeff[b] * p_virtual[b, t]^2 for b in axes(p_virtual, 1),
                    t in 1:T
                )
            )
        else
            error("Unsupported virtual_cost_coeff format")
        end
        penalty_expr = model[:virtual_penalty]
    end

    @objective(model, Min,
        model[:generation_cost] + penalty_expr
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
