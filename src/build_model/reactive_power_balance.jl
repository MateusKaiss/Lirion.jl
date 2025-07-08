# ---------------------------------------------------------------------
# Reactive Balance Expression
# ---------------------------------------------------------------------
function compute_reactive_power_balance(model, matpower_data, params, bus_index, t)
    q_gen = sum(
        model[:q_generation][g, t] for
        g in get(params.bus_to_generators, bus_index, Int[])
        if params.generator_status[g] == 1;
        init = 0.0,
    )
    q_bat = sum(
        model[:q_battery_injected][b, t] for
        b in get(params.bus_to_batteries, bus_index, Int[]);
        init = 0.0,
    )
    q_pros = sum(
        model[:q_prosumer_injected][p, t] for
        p in get(params.bus_to_prosumers, bus_index, Int[]);
        init = 0.0,
    )
    q_out = sum(
        model[:q_flow_branches][l, t] for
        l in get(params.bus_to_branches_from, bus_index, Int[])
        if params.branch_status[l] == 1;
        init = 0.0,
    )
    q_in = sum(
        model[:q_reverse_flow_branches][l, t] for
        l in get(params.bus_to_branches_to, bus_index, Int[])
        if params.branch_status[l] == 1;
        init = 0.0,
    )

    qd = if params.has_timeseries
        ts = matpower_data["bus_timeseries"]
        sum(
            ts["qd"][i]
            for i in eachindex(ts["bus"])
            if ts["bus"][i] == bus_index && ts["time"][i] == t;
            init = 0.0,
        )
    else
        bus_list = matpower_data["bus"]["bus_i"]
        qd_list = matpower_data["bus"]["qd"]
        @inbounds sum(
            qd_list[i] for i in eachindex(bus_list) if bus_list[i] == bus_index;
            init = 0.0,
        )
    end

    q_cons = if params.has_timeseries && haskey(matpower_data, "consumer_timeseries")
        tsc = matpower_data["consumer_timeseries"]
        sum(
            tsc["q_demand"][i]
            for i in eachindex(tsc["bus"])
            if tsc["bus"][i] == bus_index && tsc["time"][i] == t;
            init = 0.0,
        )
    elseif haskey(matpower_data, "consumer")
        c = matpower_data["consumer"]
        sum(
            c["q_demand"][i]
            for i in eachindex(c["bus"])
            if c["bus"][i] == bus_index;
            init = 0.0,
        )
    else
        0.0
    end
    return q_gen - q_bat + q_pros - qd - q_cons - q_out - q_in
end

# ---------------------------------------------------------------------
# Reactive Balance Constraints by Model Type
# ---------------------------------------------------------------------
function reactive_power_balance!(model, matpower_data, ::Type{DC}, ::Type{SubproblemBuild})
    return nothing
end

function reactive_power_balance!(model, matpower_data, ::Type{AC}, ::Type{SubproblemBuild})
    params = matpower_data["reactive_balance_params"]
    T = get(model.ext, :time_periods, 1)

    @inbounds for (idx, b) in enumerate(params.bus_ids), t in 1:T
        balance = compute_reactive_power_balance(model, matpower_data, params, b, t)
        @constraint(
            model,
            balance + params.shunt_susceptance[idx] * model[:v_magnitude][b, t]^2 == 0
        )
    end
end

function reactive_power_balance!(
    model,
    matpower_data,
    ::Type{SOCP},
    ::Type{SubproblemBuild},
)
    params = matpower_data["reactive_balance_params"]
    T = get(model.ext, :time_periods, 1)

    @inbounds for (idx, b) in enumerate(params.bus_ids), t in 1:T
        balance = compute_reactive_power_balance(model, matpower_data, params, b, t)
        @constraint(
            model,
            balance + params.shunt_susceptance[idx] * model[:v_magnitude_sq][b, t] == 0
        )
    end
end
