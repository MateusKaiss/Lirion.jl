
# =======================================================================
# Dispatch by model type (DC / AC / SOCP)
# =======================================================================

function active_power_balance!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{DC},
)
    pf = sys.powerflow
    T = get(model.ext, :time_periods, 1)

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
)
    pf = sys.powerflow
    T = get(model.ext, :time_periods, 1)

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
)
    pf = sys.powerflow
    T = get(model.ext, :time_periods, 1)

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

# =======================================================================
# Dispatch by model type (DC / AC / SOCP)
# =======================================================================

function reactive_power_balance!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{DC},
)
    return nothing
end

function reactive_power_balance!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{AC},
)
    pf = sys.powerflow
    T = get(model.ext, :time_periods, 1)

    valid_bus_ids = [
        b for b in pf.bus_ids if !pf.is_virtual[pf.bus_idx_map[b]]
    ]

    @inbounds for b in valid_bus_ids, t in 1:T
        idx = pf.bus_idx_map[b]
        balance = compute_reactive_power_balance(model, sys, b, t)
        @constraint(model,
            balance +
            pf.shunt_susceptance[idx] * model[:v_magnitude][b, t]^2 == 0
        )
    end

    return nothing
end

function reactive_power_balance!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{SOCP},
)
    pf = sys.powerflow
    T = get(model.ext, :time_periods, 1)

    valid_bus_ids = [
        b for b in pf.bus_ids if !pf.is_virtual[pf.bus_idx_map[b]]
    ]

    @inbounds for b in valid_bus_ids, t in 1:T
        idx = pf.bus_idx_map[b]
        balance = compute_reactive_power_balance(model, sys, b, t)
        @constraint(model,
            balance +
            pf.shunt_susceptance[idx] * model[:v_magnitude_sq][b, t] == 0
        )
    end

    return nothing
end
