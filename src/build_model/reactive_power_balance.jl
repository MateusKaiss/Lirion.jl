function get_reactive_power_parameters(matpower_data)
    has_storage = haskey(matpower_data, "storage")
    has_prosumer = haskey(matpower_data, "prosumer")

    return (
        number_of_buses = length(matpower_data["bus"]["bus_i"]),
        generator_bus_index = matpower_data["gen"]["gen_bus"],
        battery_bus_index = has_storage ? matpower_data["storage"]["storage_bus"] : Int[],
        prosumer_bus_index = has_prosumer ? matpower_data["prosumer"]["bus"] : Int[],
        demand_bus_index = matpower_data["bus"]["bus_i"],
        branch_bus_from = matpower_data["branch"]["f_bus"],
        branch_bus_to = matpower_data["branch"]["t_bus"],
        shunt_susceptance = matpower_data["bus"]["bs"],
        has_timeseries = haskey(matpower_data, "bus_timeseries"),
    )
end

function get_reactive_demand_at_time(data::Dict{String, Any}, bus_index::Int, t::Int)
    if haskey(data, "bus_timeseries")
        df = data["bus_timeseries"]
        return sum(
            row.qd for row in eachrow(df)
            if row.bus == bus_index && row.time == t;
            init = 0.0,
        )
    else
        return sum(
            data["bus"]["qd"][d] for d in eachindex(data["bus"]["bus_i"])
            if data["bus"]["bus_i"][d] == bus_index;
            init = 0.0,
        )
    end
end

function compute_reactive_power_balance(model, matpower_data, bus_index, t)
    params = get_reactive_power_parameters(matpower_data)

    return sum(
        model[:q_generation][g, t] for g in eachindex(params.generator_bus_index)
        if params.generator_bus_index[g] == bus_index;
        init = 0.0,
    ) -
           sum(
        model[:q_battery_injected][bat, t] for bat in eachindex(params.battery_bus_index)
        if params.battery_bus_index[bat] == bus_index;
        init = 0.0,
    ) +
           sum(
        model[:q_prosumer_injected][p, t] for p in eachindex(params.prosumer_bus_index)
        if params.prosumer_bus_index[p] == bus_index;
        init = 0.0,
    ) -
           get_reactive_demand_at_time(matpower_data, bus_index, t) -
           sum(
        model[:q_flow_branches][l, t] for l in eachindex(params.branch_bus_from)
        if params.branch_bus_from[l] == bus_index;
        init = 0.0,
    ) -
           sum(
        model[:q_reverse_flow_branches][l, t] for l in eachindex(params.branch_bus_to)
        if params.branch_bus_to[l] == bus_index;
        init = 0.0,
    )
end

function reactive_power_balance!(model, matpower_data, ::Type{DC}, ::Type{SubproblemBuild})
    return nothing
end

function reactive_power_balance!(model, matpower_data, ::Type{AC}, ::Type{SubproblemBuild})
    params = get_reactive_power_parameters(matpower_data)
    T = get(model.ext, :time_periods, 1)

    for b in 1:params.number_of_buses, t in 1:T
        reactive_balance = compute_reactive_power_balance(model, matpower_data, b, t)
        @constraint(
            model,
            reactive_balance + params.shunt_susceptance[b] * model[:v_magnitude][b, t]^2 ==
            0
        )
    end
end

function reactive_power_balance!(
    model,
    matpower_data,
    ::Type{SOCP},
    ::Type{SubproblemBuild},
)
    params = get_reactive_power_parameters(matpower_data)
    T = get(model.ext, :time_periods, 1)

    for b in 1:params.number_of_buses, t in 1:T
        reactive_balance = compute_reactive_power_balance(model, matpower_data, b, t)
        @constraint(
            model,
            reactive_balance + params.shunt_susceptance[b] * model[:v_magnitude_sq][b, t] ==
            0
        )
    end
end
