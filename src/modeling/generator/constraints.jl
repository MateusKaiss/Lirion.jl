# ============================================================
# Generator constraints
# ============================================================

function add_curtailment_constraints!(model::JuMP.Model, gen::GenData)
    T = get(model.ext, :time_periods, 1)
    pmax_ts = gen.pmax_ts

    @constraint(model, [g in gen.curtailable, t in 1:T],
        model[:p_generation][g, t] + model[:curtailment][g, t] == pmax_ts[g, t]
    )

    return nothing
end
