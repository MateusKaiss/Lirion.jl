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

function create_prosumer_constraints!(model::JuMP.Model, pros::ProsumerData)
    T = get(model.ext, :time_periods, 1)
    np = pros.number_of_prosumers

    # Pré-computa capacidade² para não elevar ao quadrado dentro do macro
    inverter_cap_squared = pros.inverter_capacity .^ 2

    # Ativo: p_inj = p_used - demanda
    @constraint(model, [p in 1:np, t in 1:T],
        model[:p_prosumer_injected][p, t] ==
        model[:pv_generation_used][p, t] - pros.demand[p]
    )

    # Reativo: q_inj = q_gen - q_dem
    @constraint(model, [p in 1:np, t in 1:T],
        model[:q_prosumer_injected][p, t] ==
        pros.q_generation[p] - pros.q_demand[p]
    )

    # Balanço da geração FV: usada + corte = geração disponível
    @constraint(model, [p in 1:np, t in 1:T],
        model[:pv_generation_used][p, t] +
        model[:pv_curtailment][p, t] ==
        pros.generation[p]
    )

    # Limite de capacidade do inversor
    @constraint(model, [p in 1:np, t in 1:T],
        model[:p_prosumer_injected][p, t]^2 +
        model[:q_prosumer_injected][p, t]^2 <=
        inverter_cap_squared[p]
    )

    return nothing
end

# =======================================================================
# Dispatch by model type (DC / AC / SOCP)
# =======================================================================
function prosumer!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{T},
    ::Type{SubproblemBuild},
) where {T <: OPFModel}
    pros = sys.prosumer

    if pros.number_of_prosumers == 0
        @info "No prosumers found. Skipping prosumer modeling."
        return nothing
    end

    create_prosumer_variables!(model, pros)
    create_prosumer_constraints!(model, pros)

    return nothing
end
