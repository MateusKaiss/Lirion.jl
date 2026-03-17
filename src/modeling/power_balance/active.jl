function compute_active_power_balance(
    model::JuMP.Model,
    sys::SystemData,
    bus_id::Int,
    t::Int,
)
    pf = sys.powerflow
    gen = sys.gen
    bus = sys.bus
    bus_ts = sys.bus_ts
    consumer_ts = sys.consumer_ts

    bus_idx = pf.bus_idx_map[bus_id]

    p_gen = sum(
        model[:p_generation][g, t]
        for g in pf.bus_to_generators[bus_idx]
        if gen.existing[g];
        init = 0.0,
    )

    p_bat =
        haskey(model, :p_battery_injected) ?
        sum(
            model[:p_battery_injected][b, t]
            for b in pf.bus_to_batteries[bus_idx];
            init = 0.0,
        ) : 0.0

    p_pros =
        haskey(model, :p_prosumer_injected) ?
        sum(
            model[:p_prosumer_injected][p, t]
            for p in pf.bus_to_prosumers[bus_idx];
            init = 0.0,
        ) : 0.0

    p_out =
        haskey(model, :p_flow_branches) ?
        sum(
            model[:p_flow_branches][l, t]
            for l in pf.bus_to_branches_from[bus_idx]
            if pf.branch_status[l] == 1;
            init = 0.0,
        ) : 0.0

    p_in =
        haskey(model, :p_reverse_flow_branches) ?
        sum(
            model[:p_reverse_flow_branches][l, t]
            for l in pf.bus_to_branches_to[bus_idx]
            if pf.branch_status[l] == 1;
            init = 0.0,
        ) : 0.0

    pd = 0.0
    if bus_ts.has_timeseries && !isempty(bus_ts.time)
        @inbounds for i in eachindex(bus_ts.bus)
            if bus_ts.bus[i] == bus_id && bus_ts.time[i] == t
                pd += bus_ts.pd[i]
            end
        end
    else
        pd = bus.pd[bus_idx]
    end

    p_cons = 0.0
    if consumer_ts.has_timeseries && !isempty(consumer_ts.time)
        @inbounds for i in eachindex(consumer_ts.bus)
            if consumer_ts.bus[i] == bus_id && consumer_ts.time[i] == t
                p_cons += consumer_ts.demand[i]
            end
        end
    end

    return p_gen - p_bat - p_pros - pd - p_cons - p_out - p_in
end
