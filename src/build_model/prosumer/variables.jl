# =====================================================================
# Prosumers (usa ProsumerData em SystemData)
# =====================================================================

function create_prosumer_variables!(model::JuMP.Model, pros::ProsumerData)
    T = get(model.ext, :time_periods, 1)
    np = pros.number_of_prosumers

    @variable(
        model,
        p_prosumer_injected[p in 1:np, t in 1:T],
    )

    @variable(
        model,
        q_prosumer_injected[p in 1:np, t in 1:T],
    )

    @variable(
        model,
        pv_generation_used[p in 1:np, t in 1:T],
        lower_bound = 0.0,
    )

    @variable(
        model,
        pv_curtailment[p in 1:np, t in 1:T],
        lower_bound = 0.0,
    )

    # JuMP já registra por nome, mas se quiser explitar:
    model[:p_prosumer_injected] = p_prosumer_injected
    model[:q_prosumer_injected] = q_prosumer_injected
    model[:pv_generation_used] = pv_generation_used
    model[:pv_curtailment] = pv_curtailment

    return nothing
end
