function compute_power_balance(
    model::JuMP.Model,
    sys::SystemData,
    bus_id::Int,
    t::Int,
)
    pf          = sys.powerflow
    gen         = sys.gen
    bus         = sys.bus
    bus_ts      = sys.bus_ts
    consumer_ts = sys.consumer_ts

    # --------------------------------------------------------------
    # Geração Ativa
    # --------------------------------------------------------------
    p_gen = sum(
        model[:p_generation][g, t]
        for g in get(pf.bus_to_generators, bus_id, Int[])
        if gen.existing[g];
        init = 0.0,
    )

    # --------------------------------------------------------------
    # Baterias
    # --------------------------------------------------------------
    p_bat = if haskey(model, :p_battery_injected)
        sum(
            model[:p_battery_injected][b, t] for
            b in get(pf.bus_to_batteries, bus_id, Int[]);
            init = 0.0,
        )
    else
        0.0
    end

    # --------------------------------------------------------------
    # Prosumers
    # --------------------------------------------------------------
    p_pros = if haskey(model, :p_prosumer_injected)
        sum(
            model[:p_prosumer_injected][p, t] for
            p in get(pf.bus_to_prosumers, bus_id, Int[]);
            init = 0.0,
        )
    else
        0.0
    end

    # --------------------------------------------------------------
    # Fluxos de ramo
    # --------------------------------------------------------------
    p_out = if haskey(model, :p_flow_branches)
        sum(
            model[:p_flow_branches][l, t] for
            l in get(pf.bus_to_branches_from, bus_id, Int[])
            if pf.branch_status[l] == 1;
            init = 0.0,
        )
    else
        0.0
    end

    p_in = if haskey(model, :p_reverse_flow_branches)
        sum(
            model[:p_reverse_flow_branches][l, t] for
            l in get(pf.bus_to_branches_to, bus_id, Int[])
            if pf.branch_status[l] == 1;
            init = 0.0,
        )
    else
        0.0
    end

    # --------------------------------------------------------------
    # Carga de barra (pd): timeseries OU estática
    # --------------------------------------------------------------
    pd = 0.0
    if bus_ts.has_timeseries && !isempty(bus_ts.time)
        @inbounds for i in eachindex(bus_ts.bus)
            if bus_ts.bus[i] == bus_id && bus_ts.time[i] == t
                pd += bus_ts.pd[i]
            end
        end
    else
        idx = pf.bus_idx_map[bus_id]
        pd  = bus.pd[idx]
    end

    # --------------------------------------------------------------
    # Carga adicional de "consumer" (ConsumerTimeseriesData)
    # --------------------------------------------------------------
    p_cons = 0.0
    if consumer_ts.has_timeseries && !isempty(consumer_ts.time)
        @inbounds for i in eachindex(consumer_ts.bus)
            if consumer_ts.bus[i] == bus_id && consumer_ts.time[i] == t
                p_cons += consumer_ts.demand[i]
            end
        end
    end

    if haskey(model.ext, :p_virtual) &&
       bus_id in axes(model.ext[:p_virtual], 1) &&
       t      in axes(model.ext[:p_virtual], 2)
        p_virt = model.ext[:p_virtual][bus_id, t]
    else
        p_virt = 0.0
    end

    return p_gen - p_bat - p_pros - pd - p_cons - p_out - p_in
end

# =======================================================================
# Dispatch by model type (DC / AC / SOCP)
# =======================================================================
function active_power_balance!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{DC},
    ::Type{SubproblemBuild},
)
    pf = sys.powerflow
    T  = get(model.ext, :time_periods, 1)

    valid_bus_ids = [
        b for b in pf.bus_ids if !pf.is_virtual[pf.bus_idx_map[b]]
    ]

    @constraint(model, [b in valid_bus_ids, t in 1:T],
        compute_power_balance(model, sys, b, t) -
        pf.shunt_conductance[pf.bus_idx_map[b]] == 0
    )

    return nothing
end

function active_power_balance!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{AC},
    ::Type{SubproblemBuild},
)
    pf = sys.powerflow
    T  = get(model.ext, :time_periods, 1)

    valid_bus_ids = [
        b for b in pf.bus_ids if !pf.is_virtual[pf.bus_idx_map[b]]
    ]

    @constraint(model, [b in valid_bus_ids, t in 1:T],
        compute_power_balance(model, sys, b, t) -
        pf.shunt_conductance[pf.bus_idx_map[b]] *
        model[:v_magnitude][b, t]^2 == 0
    )

    return nothing
end

function active_power_balance!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{SOCP},
    ::Type{SubproblemBuild},
)
    pf = sys.powerflow
    T  = get(model.ext, :time_periods, 1)

    valid_bus_ids = [
        b for b in pf.bus_ids if !pf.is_virtual[pf.bus_idx_map[b]]
    ]

    @constraint(model, [b in valid_bus_ids, t in 1:T],
        compute_power_balance(model, sys, b, t) -
        pf.shunt_conductance[pf.bus_idx_map[b]] *
        model[:v_magnitude_sq][b, t] == 0
    )

    return nothing
end