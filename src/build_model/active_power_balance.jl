function compute_power_balance(model, matpower_data, params, bus_index, t)
    p_gen = sum(
        model[:p_generation][g, t]
        for g in get(params.bus_to_generators, bus_index, Int[])
        if params.existing[g] == 1;
        init = 0.0,
    )
    p_bat = sum(
        model[:p_battery_injected][b, t] for
        b in get(params.bus_to_batteries, bus_index, Int[]);
        init = 0.0,
    )

    p_pros = sum(
        model[:p_prosumer_injected][p, t] for
        p in get(params.bus_to_prosumers, bus_index, Int[]);
        init = 0.0,
    )
    p_out = sum(
        model[:p_flow_branches][l, t] for
        l in get(params.bus_to_branches_from, bus_index, Int[])
        if params.branch_status[l] == 1;
        init = 0.0,
    )
    p_in = sum(
        model[:p_reverse_flow_branches][l, t] for
        l in get(params.bus_to_branches_to, bus_index, Int[])
        if params.branch_status[l] == 1;
        init = 0.0,
    )

    pd = if params.has_timeseries
        ts = matpower_data["bus_timeseries"]
        sum(
            ts["pd"][i]
            for i in eachindex(ts["bus"])
            if ts["bus"][i] == bus_index && ts["time"][i] == t;
            init = 0.0,
        )
    else
        bus_list = matpower_data["bus"]["bus_i"]
        pd_list = matpower_data["bus"]["pd"]
        @inbounds sum(
            pd_list[i] for i in eachindex(bus_list) if bus_list[i] == bus_index;
            init = 0.0,
        )
    end

    p_cons = if params.has_timeseries
        tsc = matpower_data["consumer_timeseries"]
        sum(
            tsc["demand"][i]
            for i in eachindex(tsc["bus"])
            if tsc["bus"][i] == bus_index && tsc["time"][i] == t;
            init = 0.0,
        )
    else
        c = matpower_data["consumer"]
        sum(
            c["demand"][i]
            for i in eachindex(c["bus"])
            if c["bus"][i] == bus_index;
            init = 0.0,
        )
    end

    return p_gen - p_bat - p_pros - pd - p_cons - p_out - p_in
end

# ---------------------------------------------------------------------
# Dispatching Active Power Balance Variables by Model Type
# ---------------------------------------------------------------------
function active_power_balance!(model, matpower_data, ::Type{DC}, ::Type{SubproblemBuild})
    params = matpower_data["power_balance_params"]
    T = get(model.ext, :time_periods, 1)

    @constraint(model, [b in params.bus_ids, t in 1:T],
        compute_power_balance(model, matpower_data, params, b, t) -
        params.shunt_conductance[params.bus_idx_map[b]] == 0
    )
end

function active_power_balance!(model, matpower_data, ::Type{AC}, ::Type{SubproblemBuild})
    params = matpower_data["power_balance_params"]
    T = get(model.ext, :time_periods, 1)

    @constraint(model, [b in params.bus_ids, t in 1:T],
        compute_power_balance(model, matpower_data, params, b, t) -
        params.shunt_conductance[params.bus_idx_map[b]] * model[:v_magnitude][b, t]^2 == 0
    )
end

function active_power_balance!(model, matpower_data, ::Type{SOCP}, ::Type{SubproblemBuild})
    params = matpower_data["power_balance_params"]
    T = get(model.ext, :time_periods, 1)

    @constraint(model, [b in params.bus_ids, t in 1:T],
        compute_power_balance(model, matpower_data, params, b, t) -
        params.shunt_conductance[params.bus_idx_map[b]] * model[:v_magnitude_sq][b, t] == 0
    )
end
