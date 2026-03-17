function compute_reactive_power_balance(
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

    q_gen =
        haskey(model, :q_generation) ?
        sum(
            model[:q_generation][g, t]
            for g in pf.bus_to_generators[bus_idx]
            if gen.existing[g];
            init = 0.0,
        ) : 0.0

    q_bat =
        haskey(model, :q_battery_injected) ?
        sum(
            model[:q_battery_injected][b, t]
            for b in pf.bus_to_batteries[bus_idx];
            init = 0.0,
        ) : 0.0

    q_pros =
        haskey(model, :q_prosumer_injected) ?
        sum(
            model[:q_prosumer_injected][p, t]
            for p in pf.bus_to_prosumers[bus_idx];
            init = 0.0,
        ) : 0.0

    q_out =
        haskey(model, :q_flow_branches) ?
        sum(
            model[:q_flow_branches][l, t]
            for l in pf.bus_to_branches_from[bus_idx]
            if pf.branch_status[l] == 1;
            init = 0.0,
        ) : 0.0

    q_in =
        haskey(model, :q_reverse_flow_branches) ?
        sum(
            model[:q_reverse_flow_branches][l, t]
            for l in pf.bus_to_branches_to[bus_idx]
            if pf.branch_status[l] == 1;
            init = 0.0,
        ) : 0.0

    qd = 0.0
    if bus_ts.has_timeseries && !isempty(bus_ts.time)
        @inbounds for i in eachindex(bus_ts.bus)
            if bus_ts.bus[i] == bus_id && bus_ts.time[i] == t
                qd += bus_ts.qd[i]
            end
        end
    else
        qd = bus.qd[bus_idx]
    end

    q_cons = 0.0
    if consumer_ts.has_timeseries && !isempty(consumer_ts.time)
        @inbounds for i in eachindex(consumer_ts.bus)
            if consumer_ts.bus[i] == bus_id && consumer_ts.time[i] == t
                q_cons += consumer_ts.q_demand[i]
            end
        end
    end

    return q_gen - q_bat + q_pros - qd - q_cons - q_out - q_in
end
