function add_curtailment_constraints!(model::JuMP.Model, gen::GenData)
    T = get(model.ext, :time_periods, 1)
    pmax = gen.pmax

    @constraint(model, [g in gen.curtailable, t in 1:T],
        model[:p_generation][g, t] + model[:curtailment][g, t] == pmax[g]
    )

    return nothing
end
