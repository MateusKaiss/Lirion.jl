function active_power_balance!(
    model::JuMP.Model,
    sys::SystemData,
    model_type::Type{<:OPFModel},
)
    add_active_power_balance_constraints!(model, sys, model_type)
    return nothing
end

function reactive_power_balance!(
    model::JuMP.Model,
    sys::SystemData,
    model_type::Type{<:OPFModel},
)
    add_reactive_power_balance_constraints!(model, sys, model_type)
    return nothing
end

# =======================================================================
# ACTIVE POWER
# =======================================================================

function add_active_power_balance_constraints!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{DC},
)
    pf = sys.powerflow
    T = get(model.ext, :time_periods, 1)
    valid_bus_ids = pf.non_virtual_bus_ids
    gs = pf.shunt_conductance
    bus_idx_map = pf.bus_idx_map

    @constraint(model, [b in valid_bus_ids, t in 1:T],
        compute_active_power_balance(model, sys, b, t) -
        gs[bus_idx_map[b]] == 0
    )

    return nothing
end

function add_active_power_balance_constraints!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{AC},
)
    pf = sys.powerflow
    T = get(model.ext, :time_periods, 1)
    valid_bus_ids = pf.non_virtual_bus_ids
    gs = pf.shunt_conductance
    bus_idx_map = pf.bus_idx_map

    @constraint(model, [b in valid_bus_ids, t in 1:T],
        compute_active_power_balance(model, sys, b, t) -
        gs[bus_idx_map[b]] * model[:v_magnitude][b, t]^2 == 0
    )

    return nothing
end

function add_active_power_balance_constraints!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{ACR},
)
    pf = sys.powerflow
    T = get(model.ext, :time_periods, 1)
    valid_bus_ids = pf.non_virtual_bus_ids
    gs = pf.shunt_conductance
    bus_idx_map = pf.bus_idx_map

    @constraint(model, [b in valid_bus_ids, t in 1:T],
        compute_active_power_balance(model, sys, b, t) -
        gs[bus_idx_map[b]] *
        (model[:v_real][b, t]^2 + model[:v_imag][b, t]^2) == 0
    )

    return nothing
end

function add_active_power_balance_constraints!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{SOCWR},
)
    pf = sys.powerflow
    T = get(model.ext, :time_periods, 1)
    valid_bus_ids = pf.non_virtual_bus_ids
    gs = pf.shunt_conductance
    bus_idx_map = pf.bus_idx_map

    @constraint(model, [b in valid_bus_ids, t in 1:T],
        compute_active_power_balance(model, sys, b, t) -
        gs[bus_idx_map[b]] * model[:v_magnitude_sq][b, t] == 0
    )

    return nothing
end

# =======================================================================
# REACTIVE POWER
# =======================================================================

function add_reactive_power_balance_constraints!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{DC},
)
    return nothing
end

function add_reactive_power_balance_constraints!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{AC},
)
    pf = sys.powerflow
    T = get(model.ext, :time_periods, 1)
    valid_bus_ids = pf.non_virtual_bus_ids
    bs = pf.shunt_susceptance
    bus_idx_map = pf.bus_idx_map

    @constraint(model, [b in valid_bus_ids, t in 1:T],
        compute_reactive_power_balance(model, sys, b, t) +
        bs[bus_idx_map[b]] * model[:v_magnitude][b, t]^2 == 0
    )

    return nothing
end

function add_reactive_power_balance_constraints!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{ACR},
)
    pf = sys.powerflow
    T = get(model.ext, :time_periods, 1)
    valid_bus_ids = pf.non_virtual_bus_ids
    bs = pf.shunt_susceptance
    bus_idx_map = pf.bus_idx_map

    @constraint(model, [b in valid_bus_ids, t in 1:T],
        compute_reactive_power_balance(model, sys, b, t) +
        bs[bus_idx_map[b]] *
        (model[:v_real][b, t]^2 + model[:v_imag][b, t]^2) == 0
    )

    return nothing
end

function add_reactive_power_balance_constraints!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{SOCWR},
)
    pf = sys.powerflow
    T = get(model.ext, :time_periods, 1)
    valid_bus_ids = pf.non_virtual_bus_ids
    bs = pf.shunt_susceptance
    bus_idx_map = pf.bus_idx_map

    @constraint(model, [b in valid_bus_ids, t in 1:T],
        compute_reactive_power_balance(model, sys, b, t) +
        bs[bus_idx_map[b]] * model[:v_magnitude_sq][b, t] == 0
    )

    return nothing
end
