function compute_reactive_power_balance(model, matpower_data, params, bus_index, t)
    q_gen = sum(
        model[:q_generation][g, t]
        for g in get(params.bus_to_generators, bus_index, Int[])
        if params.generator_status[g] == 1; init = 0.0)

    q_bat = sum(
        model[:q_battery_injected][b, t]
        for b in get(params.bus_to_batteries, bus_index, Int[]); init = 0.0)

    q_pros = sum(
        model[:q_prosumer_injected][p, t]
        for p in get(params.bus_to_prosumers, bus_index, Int[]); init = 0.0)

    qd = if params.has_timeseries
        ts = matpower_data["bus_timeseries"]
        sum(
            ts["qd"][i] for i in eachindex(ts["bus"])
            if ts["bus"][i] == bus_index && ts["time"][i] == t; init = 0.0)
    else
        bus_list = matpower_data["bus"]["bus_i"]
        qd_list = matpower_data["bus"]["qd"]
        @inbounds sum(
            qd_list[i] for i in eachindex(bus_list)
            if bus_list[i] == bus_index; init = 0.0)
    end

    q_cons = if params.has_timeseries && haskey(matpower_data, "consumer_timeseries")
        tsc = matpower_data["consumer_timeseries"]
        sum(
            tsc["q_demand"][i] for i in eachindex(tsc["bus"])
            if tsc["bus"][i] == bus_index && tsc["time"][i] == t; init = 0.0)
    elseif haskey(matpower_data, "consumer")
        c = matpower_data["consumer"]
        sum(c["q_demand"][i] for i in eachindex(c["bus"])
            if c["bus"][i] == bus_index; init = 0.0)
    else
        0.0
    end

    sum_q_out = _sum_q_out(model, params, bus_index, t)

    return -q_gen + sum_q_out + qd + q_cons - q_bat + q_pros
end

function reactive_power_balance!(model, matpower_data, ::Type{DC}, ::Type{SubproblemBuild})
    return nothing
end

function reactive_power_balance!(model, matpower_data, ::Type{AC}, ::Type{SubproblemBuild})
    params = matpower_data["reactive_balance_params"]
    T = get(model.ext, :time_periods, 1)
    valid_bus_ids = filter(b -> !params.is_virtual[params.bus_idx_map[b]], params.bus_ids)

    @constraint(model, [b in valid_bus_ids, t in 1:T],
        compute_reactive_power_balance(model, matpower_data, params, b, t)
        +
        params.shunt_susceptance[params.bus_idx_map[b]] * model[:v_magnitude][b, t]^2 == 0
    )
end

function reactive_power_balance!(
    model,
    matpower_data,
    ::Type{SOCP},
    ::Type{SubproblemBuild},
)
    params = matpower_data["reactive_balance_params"]
    T = get(model.ext, :time_periods, 1)
    valid_bus_ids = filter(b -> !params.is_virtual[params.bus_idx_map[b]], params.bus_ids)

    @constraint(model, [b in valid_bus_ids, t in 1:T],
        compute_reactive_power_balance(model, matpower_data, params, b, t)
        -
        params.shunt_susceptance[params.bus_idx_map[b]] * model[:v_magnitude_sq][b, t] == 0
    )
end
