function compute_reactive_power_balance(
    model::JuMP.Model,
    sys::SystemData,
    bus_id::Int,
    t::Int,
)
    pf = sys.powerflow
    gen = sys.gen
    bus = sys.bus

    bus_idx = pf.bus_idx_map[bus_id]

    q_gen =
        haskey(model, :q_generation) ?
        sum(
            model[:q_generation][g, t]
            for g in pf.bus_to_generators[bus_idx]
            if gen.existing[g];
            init = 0.0,
        ) : 0.0

    q_stor =
        haskey(model, :q_storage_injected) ?
        sum(
            model[:q_storage_injected][b, t]
            for b in pf.bus_to_storages[bus_idx];
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

    qd = pf.has_bus_ts ? pf.qd_ts[bus_idx, t] : bus.qd[bus_idx]

    return q_gen - q_stor - qd - q_out - q_in
end
