function get_reactive_power_parameters(matpower_data)
    has_storage = haskey(matpower_data, "storage")

    return (
        number_of_buses = length(matpower_data["bus"]["bus_i"]),
        generator_bus_index = matpower_data["gen"]["gen_bus"],
        battery_bus_index = has_storage ? matpower_data["storage"]["storage_bus"] : Int[],
        demand_bus_index = matpower_data["bus"]["bus_i"],
        branch_bus_from = matpower_data["branch"]["f_bus"],
        branch_bus_to = matpower_data["branch"]["t_bus"],
        shunt_susceptance = matpower_data["bus"]["bs"],
        reactive_demand = matpower_data["bus"]["qd"],
    )
end

function compute_reactive_power_balance(model, matpower_data, bus_index)
    params = get_reactive_power_parameters(matpower_data)

    return sum(
        model[:q_generation][g] for g in eachindex(params.generator_bus_index)
        if params.generator_bus_index[g] == bus_index;
        init = 0.0,
    ) -
           sum(
        model[:q_battery_injected][bat] for bat in eachindex(params.battery_bus_index)
        if params.battery_bus_index[bat] == bus_index;
        init = 0.0,
    ) -
           sum(
        params.reactive_demand[d] for d in eachindex(params.demand_bus_index)
        if params.demand_bus_index[d] == bus_index;
        init = 0.0,
    ) -
           sum(
        model[:q_flow_branches][l] for l in eachindex(params.branch_bus_from)
        if params.branch_bus_from[l] == bus_index;
        init = 0.0,
    ) -
           sum(
        model[:q_reverse_flow_branches][l] for l in eachindex(params.branch_bus_to)
        if params.branch_bus_to[l] == bus_index;
        init = 0.0,
    )
end

function reactive_power_balance!(model, matpower_data, ::Type{DC})
    return nothing
end

function reactive_power_balance!(model, matpower_data, ::Type{AC})
    params = get_reactive_power_parameters(matpower_data)

    for b in 1:params.number_of_buses
        reactive_balance = compute_reactive_power_balance(model, matpower_data, b)
        @constraint(
            model,
            reactive_balance - params.shunt_susceptance[b] * model[:v_magnitude][b]^2 == 0
        )
    end
end

function reactive_power_balance!(model, matpower_data, ::Type{SOCP})
    params = get_reactive_power_parameters(matpower_data)

    for b in 1:params.number_of_buses
        reactive_balance = compute_reactive_power_balance(model, matpower_data, b)
        @constraint(
            model,
            reactive_balance - params.shunt_susceptance[b] * model[:v_magnitude_sq][b] == 0
        )
    end
end