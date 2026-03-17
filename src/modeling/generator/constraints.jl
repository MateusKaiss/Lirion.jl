function add_curtailment_constraints!(model::JuMP.Model, gen::GenData)
    T = get(model.ext, :time_periods, 1)
    G = gen.n

    curtailable = [g for g in 1:G if has_curtailment(gen, g)]

    @constraint(model, [g in curtailable, t in 1:T],
        model[:p_generation][g, t] + model[:curtailment][g, t] == gen.pmax[g]
    )

    return nothing
end
