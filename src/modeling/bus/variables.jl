@inline function initial_bus_angle(bus::BusData)::Float64
    ref = bus.reference_bus
    ref == 0 && return 0.0
    idx = bus.id_to_idx[ref]
    idx == 0 && return 0.0
    va0 = bus.va_start[idx]
    return isnan(va0) ? 0.0 : va0
end

@inline function normalize_initial_angle(va_b::Float64, theta_ref::Float64)
    return isnan(va_b) ? 0.0 : va_b - theta_ref
end

function add_voltage_angle_variables!(model::JuMP.Model, bus::BusData)
    T = get(model.ext, :time_periods, 1)
    bus_ids = bus.id
    va_start = bus.va_start
    id_to_idx = bus.id_to_idx

    @variable(model, theta[b in bus_ids, t in 1:T])

    if has_reference_bus(bus)
        ref = bus.reference_bus
        @constraint(model, [t in 1:T], theta[ref, t] == 0.0)
    end

    theta_start = initial_bus_angle(bus)

    @inbounds for b in bus_ids
        idx = id_to_idx[b]
        val = normalize_initial_angle(va_start[idx], theta_start)
        for t in 1:T
            JuMP.set_start_value(theta[b, t], val)
        end
    end

    return nothing
end

function add_voltage_magnitude_variables!(model::JuMP.Model, bus::BusData)
    T = get(model.ext, :time_periods, 1)
    bus_ids = bus.id
    id_to_idx = bus.id_to_idx
    v_min = bus.v_min
    v_max = bus.v_max
    vm_start = bus.vm_start

    @variable(
        model,
        v_magnitude[b in bus_ids, t in 1:T],
        lower_bound = v_min[id_to_idx[b]],
        upper_bound = v_max[id_to_idx[b]],
        start = isnan(vm_start[id_to_idx[b]]) ? 1.0 : vm_start[id_to_idx[b]],
    )

    return nothing
end

function add_socp_bus_variables!(model::JuMP.Model, bus::BusData)
    T = get(model.ext, :time_periods, 1)
    bus_ids = bus.id
    id_to_idx = bus.id_to_idx
    v_min = bus.v_min
    v_max = bus.v_max
    vm_start = bus.vm_start

    @variable(
        model,
        v_magnitude_sq[b in bus_ids, t in 1:T],
        lower_bound = v_min[id_to_idx[b]]^2,
        upper_bound = v_max[id_to_idx[b]]^2,
        start = (isnan(vm_start[id_to_idx[b]]) ? 1.0 : vm_start[id_to_idx[b]])^2,
    )

    return nothing
end
