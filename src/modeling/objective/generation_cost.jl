function add_generation_cost!(model::JuMP.Model, obj::ObjectiveData)
    T = get(model.ext, :time_periods, 1)
    active_generators = findall(obj.generator_status .== 1)
    costs = obj.generator_costs

    @expression(model, generation_cost,
        sum(
            begin
                coefs = costs[t][g]
                sum(
                    coefs[i] *
                    (obj.sbase * model[:p_generation][g, t])^(length(coefs) - i)
                    for i in eachindex(coefs)
                )
            end
            for t in 1:T, g in active_generators
            if t <= length(costs) && g <= length(costs[t]) && !isempty(costs[t][g])
        )
    )

    model[:generation_cost] = generation_cost

    return nothing
end
