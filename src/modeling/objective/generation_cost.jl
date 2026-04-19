function add_generation_cost!(model::JuMP.Model, obj::ObjectiveData)
    T = get(model.ext, :time_periods, 1)
    active_generators = findall(obj.generator_status .== 1)
    costs = obj.generator_costs

    @expression(model, fuel_cost,
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

    # Curtailment penalty (only if curtailment variables exist in the
    # model). We penalize post-per-unit dispatch (sbase scaling) using the
    # per-generator curtailment_cost read from the raw gen table.
    if haskey(model, :curtailment) && !isempty(obj.curtailment_cost)
        curt = model[:curtailment]
        curt_costs = obj.curtailment_cost
        @expression(model, curtailment_cost_expr,
            sum(
                curt_costs[g] * obj.sbase * curt[g, t]
                for g in axes(curt, 1), t in axes(curt, 2)
                if g <= length(curt_costs) && curt_costs[g] != 0.0
            )
        )
        model[:curtailment_cost_expr] = curtailment_cost_expr
        @expression(model, generation_cost, fuel_cost + curtailment_cost_expr)
    else
        @expression(model, generation_cost, fuel_cost)
    end

    model[:generation_cost] = generation_cost
    return nothing
end
