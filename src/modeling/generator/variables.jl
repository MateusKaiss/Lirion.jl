@inline has_curtailment(gen::GenData, g::Int) =
    gen.existing[g] && gen.plant_type[g] == "R"

# =====================================================================
# VARIABLES
# =====================================================================

function add_active_generation_variables!(model::JuMP.Model, gen::GenData)
    T = get(model.ext, :time_periods, 1)
    G = gen.n
    pg_start = gen.pg_start
    pg_fallback = gen.pg
    pmin = gen.pmin
    pmax = gen.pmax

    @variable(
        model,
        p_generation[g = 1:G, t = 1:T; gen.existing[g]],
        lower_bound = pmin[g],
        upper_bound = pmax[g],
        start = isnan(pg_start[g]) ? pg_fallback[g] : pg_start[g],
    )

    model[:p_generation] = p_generation
    return nothing
end

function add_reactive_generation_variables!(model::JuMP.Model, gen::GenData)
    T = get(model.ext, :time_periods, 1)
    G = gen.n
    qg_start = gen.qg_start
    qg_fallback = gen.qg
    qmin = gen.qmin
    qmax = gen.qmax

    @variable(
        model,
        q_generation[g = 1:G, t = 1:T; gen.existing[g]],
        lower_bound = qmin[g],
        upper_bound = qmax[g],
        start = isnan(qg_start[g]) ? qg_fallback[g] : qg_start[g],
    )

    model[:q_generation] = q_generation
    return nothing
end

function add_curtailment_variables!(model::JuMP.Model, gen::GenData)
    T = get(model.ext, :time_periods, 1)

    @variable(model, curtailment[g in gen.curtailable, t in 1:T] >= 0.0)

    model[:curtailment] = curtailment
    return nothing
end
