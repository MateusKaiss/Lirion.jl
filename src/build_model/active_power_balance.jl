function _sum_p_out(model, params, b, t)
    p_fr = model[:p_flow_branches]          # from → to   
    p_to = model[:p_reverse_flow_branches]  # to   → from 
    s = 0.0
    for l in get(params.bus_to_branches_from, b, Int[])    # ramos com from = b
        if params.branch_status[l] == 1
            s += p_fr[l, t]
        end
    end
    for l in get(params.bus_to_branches_to, b, Int[])      # ramos com to = b
        if params.branch_status[l] == 1
            s += p_to[l, t]
        end
    end
    return s
end

function _sum_q_out(model, params, b, t)
    q_fr = model[:q_flow_branches]          # 2D
    q_to = model[:q_reverse_flow_branches]  # 2D
    s = 0.0
    for l in get(params.bus_to_branches_from, b, Int[])
        if params.branch_status[l] == 1
            s += q_fr[l, t]
        end
    end
    for l in get(params.bus_to_branches_to, b, Int[])
        if params.branch_status[l] == 1
            s += q_to[l, t]
        end
    end
    return s
end

function compute_power_balance(model, matpower_data, params, bus_index, t)
    # geração
    p_gen = sum(
        model[:p_generation][g, t]
        for g in get(params.bus_to_generators, bus_index, Int[])
        if params.existing[g] == 1;
        init = 0.0,
    )
    # bateria / prosumer (injeções positivas, depois entram com sinal -)
    p_bat = sum(
        model[:p_battery_injected][b, t]
        for b in get(params.bus_to_batteries, bus_index, Int[]); init = 0.0)
    p_pros = sum(
        model[:p_prosumer_injected][p, t]
        for p in get(params.bus_to_prosumers, bus_index, Int[]); init = 0.0)

    # virtual (se houver)
    p_virt =
        (haskey(model.ext, :p_virtual) &&
         bus_index in axes(model.ext[:p_virtual], 1)) ?
        model.ext[:p_virtual][bus_index, t] : 0.0

    # demandas
    pd = if params.has_timeseries
        ts = matpower_data["bus_timeseries"]
        sum(
            ts["pd"][i] for i in eachindex(ts["bus"])
            if ts["bus"][i] == bus_index && ts["time"][i] == t; init = 0.0)
    else
        bus_list = matpower_data["bus"]["bus_i"]
        pd_list = matpower_data["bus"]["pd"]
        @inbounds sum(
            pd_list[i] for i in eachindex(bus_list)
            if bus_list[i] == bus_index; init = 0.0)
    end

    p_cons = if params.has_timeseries && haskey(matpower_data, "consumer_timeseries")
        tsc = matpower_data["consumer_timeseries"]
        sum(
            tsc["demand"][i] for i in eachindex(tsc["bus"])
            if tsc["bus"][i] == bus_index && tsc["time"][i] == t; init = 0.0)
    elseif haskey(matpower_data, "consumer")
        c = matpower_data["consumer"]
        sum(c["demand"][i] for i in eachindex(c["bus"])
            if c["bus"][i] == bus_index; init = 0.0)
    else
        0.0
    end

    sum_p_out = _sum_p_out(model, params, bus_index, t)

    return -p_gen + sum_p_out + pd + p_cons - p_bat - p_pros + p_virt
end

function active_power_balance!(model, matpower_data, ::Type{DC}, ::Type{SubproblemBuild})
    params = matpower_data["power_balance_params"]
    T = get(model.ext, :time_periods, 1)
    valid_bus_ids = filter(b -> !params.is_virtual[params.bus_idx_map[b]], params.bus_ids)

    @constraint(model, [b in valid_bus_ids, t in 1:T],
        compute_power_balance(model, matpower_data, params, b, t)
        -
        params.shunt_conductance[params.bus_idx_map[b]] == 0
    )
end

function active_power_balance!(model, matpower_data, ::Type{AC}, ::Type{SubproblemBuild})
    params = matpower_data["power_balance_params"]
    T = get(model.ext, :time_periods, 1)
    valid_bus_ids = filter(b -> !params.is_virtual[params.bus_idx_map[b]], params.bus_ids)

    @constraint(model, [b in valid_bus_ids, t in 1:T],
        compute_power_balance(model, matpower_data, params, b, t)
        -
        params.shunt_conductance[params.bus_idx_map[b]] * model[:v_magnitude][b, t]^2 == 0
    )
end

function active_power_balance!(model, matpower_data, ::Type{SOCP}, ::Type{SubproblemBuild})
    params = matpower_data["power_balance_params"]
    T = get(model.ext, :time_periods, 1)
    valid_bus_ids = filter(b -> !params.is_virtual[params.bus_idx_map[b]], params.bus_ids)

    @constraint(model, [b in valid_bus_ids, t in 1:T],
        compute_power_balance(model, matpower_data, params, b, t)
        -
        params.shunt_conductance[params.bus_idx_map[b]] * model[:v_magnitude_sq][b, t] == 0
    )
end
