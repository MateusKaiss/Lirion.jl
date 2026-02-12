@inline function has_reference_bus(bus::BusData)::Bool
    return bus.reference_bus !== nothing
end

@inline function initial_bus_angle(bus::BusData)::Float64
    ref = bus.reference_bus
    ref === nothing && return 0.0
    idx = bus.id_to_idx[ref]
    va0 = bus.vastart[idx]

    return (va0 === nothing || isnan(va0)) ? 0.0 : va0
end

@inline function normalize_initial_angle(va_b::Union{Nothing, Float64}, theta_ref::Float64)
    if va_b === nothing || isnan(va_b)
        return 0.0
    else
        return va_b - theta_ref
    end
end

function create_voltage_angle_variables!(model::JuMP.Model, bus::BusData)
    T = get(model.ext, :time_periods, 1)
    bus_ids = bus.id
    va_start = bus.vastart

    @variable(model, theta[b in bus_ids, t in 1:T])

    if has_reference_bus(bus)
        ref = bus.reference_bus::Int
        @constraint(model, [t in 1:T], theta[ref, t] == 0.0)
    end

    theta_start = initial_bus_angle(bus)

    @inbounds for b in bus_ids
        idx = bus.id_to_idx[b]
        val = normalize_initial_angle(va_start[idx], theta_start)
        for t in 1:T
            JuMP.set_start_value(theta[b, t], val)
        end
    end

    return nothing
end

function create_voltage_magnitude_variables!(model::JuMP.Model, bus::BusData)
    T = get(model.ext, :time_periods, 1)
    bus_ids = bus.id

    @variable(
        model,
        v_magnitude[b in bus_ids, t in 1:T],
        lower_bound = bus.vmin[bus.id_to_idx[b]],
        upper_bound = bus.vmax[bus.id_to_idx[b]],
        start = something(bus.vstart[bus.id_to_idx[b]], 1.0),
    )

    return nothing
end

function create_socp_bus_variables!(model::JuMP.Model, bus::BusData)
    T = get(model.ext, :time_periods, 1)
    bus_ids = bus.id

    @variable(
        model,
        v_magnitude_sq[b in bus_ids, t in 1:T],
        lower_bound = bus.vmin[bus.id_to_idx[b]]^2,
        upper_bound = bus.vmax[bus.id_to_idx[b]]^2,
        start = (something(bus.vstart[bus.id_to_idx[b]], 1.0))^2,
    )

    return nothing
end

# =======================================================================
# Dispatch by model type (DC / AC / SOCP)
# =======================================================================
function buses!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{DC},
)
    bus = sys.bus
    create_voltage_angle_variables!(model, bus)
    return nothing
end

function buses!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{AC},
)
    bus = sys.bus
    create_voltage_angle_variables!(model, bus)
    create_voltage_magnitude_variables!(model, bus)
    return nothing
end

function buses!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{SOCP},
)
    bus = sys.bus
    create_voltage_angle_variables!(model, bus)
    create_socp_bus_variables!(model, bus)
    return nothing
end
