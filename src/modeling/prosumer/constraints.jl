# =======================================================================
# CONSTRAINTS
# =======================================================================

function add_prosumer_constraints!(model::JuMP.Model, pros::ProsumerData)
    T = get(model.ext, :time_periods, 1)
    np = pros.number_of_prosumers
    inverter_cap_sq = pros.inverter_capacity_sq

    # Active power injection
    @constraint(model, [p in 1:np, t in 1:T],
        model[:p_prosumer_injected][p, t] ==
        model[:pv_generation_used][p, t] - pros.demand[p]
    )

    # Reactive power injection
    @constraint(model, [p in 1:np, t in 1:T],
        model[:q_prosumer_injected][p, t] ==
        pros.q_generation[p] - pros.q_demand[p]
    )

    # PV balance
    @constraint(model, [p in 1:np, t in 1:T],
        model[:pv_generation_used][p, t] +
        model[:pv_curtailment][p, t] ==
        pros.generation[p]
    )

    # Inverter capacity (SOC constraint)
    @constraint(model, [p in 1:np, t in 1:T],
        model[:p_prosumer_injected][p, t]^2 +
        model[:q_prosumer_injected][p, t]^2 <=
        inverter_cap_sq[p]
    )

    return nothing
end
