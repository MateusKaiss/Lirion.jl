function get_prosumer_parameters(data)
    return (
        number_of_prosumers = length(data["bus"]),
        bus_index = data["bus"],
        generation = data["generation"],
        demand = data["demand"],
        q_generation = data["q_generation"],
        q_demand = data["q_demand"],
        cost_injection = data["cost_injection"],
        inverter_capacity = data["inverter_capacity"],
    )
end

function create_prosumer_variables!(model::JuMP.Model, params)
    T = get(model.ext, :time_periods, 1)

    @variable(model, p_prosumer_injected[1:params.number_of_prosumers, 1:T])
    @variable(model, q_prosumer_injected[1:params.number_of_prosumers, 1:T])
    @variable(model, pv_generation_used[1:params.number_of_prosumers, 1:T] >= 0)
    @variable(model, pv_curtailment[1:params.number_of_prosumers, 1:T] >= 0)
end

function create_prosumer_constraints!(model::JuMP.Model, params)
    T = get(model.ext, :time_periods, 1)

    @constraint(
        model,
        [p = 1:params.number_of_prosumers, t = 1:T],
        model[:p_prosumer_injected][p, t] ==
        model[:pv_generation_used][p, t] - params.demand[p]
    )

    @constraint(
        model,
        [p = 1:params.number_of_prosumers, t = 1:T],
        model[:q_prosumer_injected][p, t] == params.q_generation[p] - params.q_demand[p]
    )

    @constraint(
        model,
        [p = 1:params.number_of_prosumers, t = 1:T],
        model[:pv_generation_used][p, t] + model[:pv_curtailment][p, t] ==
        params.generation[p]
    )

    @constraint(
        model,
        [p = 1:params.number_of_prosumers, t = 1:T],
        model[:p_prosumer_injected][p, t]^2 + model[:q_prosumer_injected][p, t]^2 <=
        params.inverter_capacity[p]^2
    )
end

function prosumer(
    model::JuMP.Model,
    data,
    ::Type{T},
    ::Type{SubproblemBuild},
) where {T <: OPFModel}
    if !haskey(data, "prosumer") || isempty(data["prosumer"])
        @info "No prosumers found. Skipping prosumer modeling."
        return
    end

    prosumer_data = data["prosumer"]
    params = get_prosumer_parameters(prosumer_data)

    create_prosumer_variables!(model, params)
    create_prosumer_constraints!(model, params)

    return nothing
end
