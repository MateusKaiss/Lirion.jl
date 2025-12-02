function compute_generation_cost(model::JuMP.Model, obj::ObjectiveData)
    T = get(model.ext, :time_periods, 1)
    active_generators = findall(obj.generator_status .== 1)

    @expression(model, generation_cost,
        sum(
            begin
                coefs = obj.generator_costs[t][g]
                sum(
                    coefs[i] *
                    (obj.sbase * model[:p_generation][g, t])^(length(coefs) - i)
                    for i in eachindex(coefs)
                )
            end
            for t in 1:T, g in active_generators
            if haskey(obj.generator_costs, t) && haskey(obj.generator_costs[t], g)
        )
    )

    model[:generation_cost] = generation_cost
    return nothing
end

# function compute_curtailment_cost(model::JuMP.Model, obj::ObjectiveData)
#     T = get(model.ext, :time_periods, 1)
#     renewable_generators = findall(
#         g -> obj.plant_type[g] == "R" && obj.generator_status[g] == 1,
#         1:obj.number_of_generators,
#     )
#
#     @expression(model, curtailment_cost,
#         sum(
#             obj.curtailment_cost[g] * model[:curtailment][g, t]
#             for g in renewable_generators, t in 1:T
#         )
#     )
#
#     model[:curtailment_cost] = curtailment_cost
#     return nothing
# end

# =======================================================================
# Dispatch by model type (DC / AC / SOCP)
# =======================================================================
function objective!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{DC},
    ::Type{SubproblemBuild},
)
    obj = sys.objective
    compute_generation_cost(model, obj)

    @objective(model, Min,
        model[:generation_cost]
    )
end

function objective!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{AC},
    ::Type{SubproblemBuild},
)
    obj = sys.objective
    compute_generation_cost(model, obj)

    @objective(model, Min,
        model[:generation_cost]
    )
end

function objective!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{SOCP},
    ::Type{SubproblemBuild},
)
    obj = sys.objective
    compute_generation_cost(model, obj)

    @objective(model, Min,
        model[:generation_cost]
    )
end