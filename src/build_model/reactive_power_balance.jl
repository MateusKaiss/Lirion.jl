function get_reactive_power_parameters(bus_data, gen_data, battery_data, demand_data, branch_data, transformer_data)
    return (
        number_of_buses = length(bus_data["bus_i"]),
        generator_bus_index = gen_data["gen_bus"],
        battery_bus_index = battery_data["storage_bus"],
        demand_bus_index = demand_data["storage_bus"],
        branch_bus_from = branch_data["f_bus"],
        branch_bus_to = branch_data["t_bus"],
        transformer_bus_from = transformer_data["f_bus"],
        transformer_bus_to = transformer_data["t_bus"],
        shunt_susceptance = bus_data["bs"],
        reactive_demand = demand_data["qd"],
    )
end

function compute_reactive_power_balance(model, bus_data, gen_data, battery_data, demand_data, branch_data, transformer_data, bus_index)
    params = get_reactive_power_parameters(bus_data, gen_data, battery_data, demand_data, branch_data, transformer_data)

    return sum(
        (
            model[:q_generation][g] for g in eachindex(params.generator_bus_index) if
            params.generator_bus_index[g] == bus_index
        );
        init = 0,
    ) -
    sum(
        (
            model[:q_battery_injected][bat] for
            bat in eachindex(params.battery_bus_index) if
            params.battery_bus_index[bat] == bus_index
        );
        init = 0,
    ) -
    sum(
        (
            params.reactive_demand[d] for d in eachindex(params.demand_bus_index) if
            params.demand_bus_index[d] == bus_index
        );
        init = 0,
    ) -
    sum(
        (
            model[:q_flow_branches][l] for
            l in eachindex(params.branch_bus_from) if params.branch_bus_from[l] == bus_index
        );
        init = 0,
    ) -
    sum(
        (
            model[:q_reverse_flow_branches][l] for
            l in eachindex(params.branch_bus_to) if params.branch_bus_to[l] == bus_index
        );
        init = 0,
    ) 
end

function reactive_power_balance!(model::JuMP.Model, bus_data, gen_data, battery_data, demand_data, branch_data, transformer_data, ::Type{DC})
    return nothing
end

function reactive_power_balance!(model::JuMP.Model, bus_data, gen_data, battery_data, demand_data, branch_data, transformer_data, ::Type{AC})
    params = get_reactive_power_parameters(bus_data, gen_data, battery_data, demand_data, branch_data, transformer_data)

    for b in 1:params.number_of_buses
        reactive_power_balance = compute_reactive_power_balance(model, bus_data, gen_data, battery_data, demand_data, branch_data, transformer_data, b)
        @constraint(
            model,
            reactive_power_balance -
            params.shunt_susceptance[b] * model[:v_magnitude][b]^2 == 0
        )
    end
end

function reactive_power_balance!(model::JuMP.Model, bus_data, gen_data, battery_data, demand_data, branch_data, transformer_data, ::Type{SOCP})
    params = get_reactive_power_parameters(bus_data, gen_data, battery_data, demand_data, branch_data, transformer_data)

    for b in 1:params.number_of_buses
        reactive_power_balance = compute_reactive_power_balance(model, bus_data, gen_data, battery_data, demand_data, branch_data, transformer_data, b)
        @constraint(
            model,
            reactive_power_balance -
            params.shunt_susceptance[b] * model[:v_magnitude_sq][b] == 0
        )
    end
end