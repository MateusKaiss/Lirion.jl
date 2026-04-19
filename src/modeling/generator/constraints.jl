function add_curtailment_constraints!(model::JuMP.Model, gen::GenData)
    T = get(model.ext, :time_periods, 1)
    pmax_ts = gen.pmax_ts

    # Renewable dispatch envelope: p_gen + curt == available capacity at t.
    # When pmax_ts[g, t] == 0 this forces both p_gen and curt to 0.
    @constraint(model, [g in gen.curtailable, t in 1:T],
        model[:p_generation][g, t] + model[:curtailment][g, t] == pmax_ts[g, t]
    )

    return nothing
end
