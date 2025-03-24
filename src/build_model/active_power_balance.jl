function get_power_balance_parameters(matpower_data)
    has_storage = haskey(matpower_data, "storage")

    return (
        number_of_buses = length(matpower_data["bus"]["bus_i"]),
        generator_bus_index = matpower_data["gen"]["gen_bus"],
        battery_bus_index = has_storage ? matpower_data["storage"]["bus"] : Int[],
        demand_bus_index = matpower_data["bus"]["bus_i"],
        branch_bus_from = matpower_data["branch"]["f_bus"],
        branch_bus_to = matpower_data["branch"]["t_bus"],
        shunt_conductance = matpower_data["bus"]["gs"],
        demand = matpower_data["bus"]["pd"],
        existing = matpower_data["gen"]["gen_status"],
    )
end

function compute_power_balance(model, matpower_data, bus_index)
    params = get_power_balance_parameters(matpower_data)

    return sum(
        model[:p_generation][g] for g in eachindex(params.generator_bus_index)
        if params.generator_bus_index[g] == bus_index && params.existing[g] == 1;
        init = 0.0
    ) +
    sum(
        model[:p_battery_injected][bat] for bat in eachindex(params.battery_bus_index)
        if params.battery_bus_index[bat] == bus_index;
        init = 0.0
    ) -
    sum(
        params.demand[d] for d in eachindex(params.demand_bus_index)
        if params.demand_bus_index[d] == bus_index;
        init = 0.0
    ) -
    sum(
        model[:p_flow_branches][l] for l in eachindex(params.branch_bus_from)
        if params.branch_bus_from[l] == bus_index;
        init = 0.0
    ) -
    sum(
        model[:p_reverse_flow_branches][l] for l in eachindex(params.branch_bus_to)
        if params.branch_bus_to[l] == bus_index;
        init = 0.0
    )
end

function active_power_balance!(model, matpower_data, ::Type{DC}, ::Type{SubproblemBuild})
    params = get_power_balance_parameters(matpower_data)

    for b in 1:params.number_of_buses
        power_balance = compute_power_balance(model, matpower_data, b)
        @constraint(model, power_balance + params.shunt_conductance[b] == 0)
    end
end

function active_power_balance!(model, matpower_data, ::Type{AC}, ::Type{SubproblemBuild})
    params = get_power_balance_parameters(matpower_data)

    for b in 1:params.number_of_buses
        power_balance = compute_power_balance(model, matpower_data, b)
        @constraint(model, power_balance + params.shunt_conductance[b] * model[:v_magnitude][b]^2 == 0)
    end
end

function active_power_balance!(model, matpower_data, ::Type{SOCP}, ::Type{SubproblemBuild})
    params = get_power_balance_parameters(matpower_data)

    for b in 1:params.number_of_buses
        power_balance = compute_power_balance(model, matpower_data, b)
        @constraint(model, power_balance + params.shunt_conductance[b] * model[:v_magnitude_sq][b] == 0)
    end
end