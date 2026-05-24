# ============================================================
# Active power balance
# ============================================================

function compute_active_power_balance(
    model::JuMP.Model,
    sys::SystemData,
    bus_id::Int,
    t::Int,
)
    pf = sys.powerflow
    gen = sys.gen
    bus = sys.bus

    bus_idx = pf.bus_idx_map[bus_id]

    p_gen = sum(
        model[:p_generation][g, t]
        for g in pf.bus_to_generators[bus_idx]
        if gen.existing[g];
        init = 0.0,
    )

    p_stor =
        haskey(model, :p_storage_injected) ?
        sum(
            model[:p_storage_injected][b, t]
            for b in pf.bus_to_storages[bus_idx];
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

    pd = pf.has_bus_ts ? pf.pd_ts[bus_idx, t] : bus.pd[bus_idx]

    return p_gen - p_stor - pd - p_out - p_in
end
