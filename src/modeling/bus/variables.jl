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

    @variable(model, v_magnitude[b in bus_ids, t in 1:T])

    @inbounds for b in bus_ids
        idx = id_to_idx[b]
        lo = v_min[idx]
        hi = v_max[idx]
        start_val = isnan(vm_start[idx]) ? 1.0 : vm_start[idx]
        for t in 1:T
            v = v_magnitude[b, t]
            JuMP.set_lower_bound(v, lo)
            JuMP.set_upper_bound(v, hi)
            JuMP.set_start_value(v, start_val)
        end
    end

    return nothing
end

function add_rectangular_bus_variables!(model::JuMP.Model, bus::BusData)
    T = get(model.ext, :time_periods, 1)
    bus_ids = bus.id
    id_to_idx = bus.id_to_idx
    v_min = bus.v_min
    v_max = bus.v_max
    vm_start = bus.vm_start
    va_start = bus.va_start

    @variable(model, v_real[b in bus_ids, t in 1:T])
    @variable(model, v_imag[b in bus_ids, t in 1:T])

    if has_reference_bus(bus)
        ref = bus.reference_bus
        @constraint(model, [t in 1:T], v_imag[ref, t] == 0.0)
    end

    theta_start = initial_bus_angle(bus)

    @inbounds for b in bus_ids
        idx = id_to_idx[b]
        vmax = v_max[idx]
        vmin_sq = v_min[idx]^2
        vmax_sq = vmax^2

        vm0 = isnan(vm_start[idx]) ? 1.0 : vm_start[idx]
        va0 = isnan(va_start[idx]) ? 0.0 : (va_start[idx] - theta_start)
        vr0 = vm0 * cos(va0)
        vi0 = vm0 * sin(va0)

        for t in 1:T
            vr_var = v_real[b, t]
            vi_var = v_imag[b, t]
            JuMP.set_lower_bound(vr_var, -vmax)
            JuMP.set_upper_bound(vr_var, vmax)
            JuMP.set_lower_bound(vi_var, -vmax)
            JuMP.set_upper_bound(vi_var, vmax)
            JuMP.set_start_value(vr_var, vr0)
            JuMP.set_start_value(vi_var, vi0)
        end

        for t in 1:T
            @constraint(model, v_real[b, t]^2 + v_imag[b, t]^2 >= vmin_sq)
            @constraint(model, v_real[b, t]^2 + v_imag[b, t]^2 <= vmax_sq)
        end
    end

    return nothing
end

function add_socwr_bus_variables!(model::JuMP.Model, bus::BusData)
    T = get(model.ext, :time_periods, 1)
    bus_ids = bus.id
    id_to_idx = bus.id_to_idx
    v_min = bus.v_min
    v_max = bus.v_max
    vm_start = bus.vm_start

    @variable(model, v_magnitude_sq[b in bus_ids, t in 1:T])

    @inbounds for b in bus_ids
        idx = id_to_idx[b]
        lo = v_min[idx]^2
        hi = v_max[idx]^2
        vm0 = isnan(vm_start[idx]) ? 1.0 : vm_start[idx]
        start_val = vm0 * vm0
        for t in 1:T
            v = v_magnitude_sq[b, t]
            JuMP.set_lower_bound(v, lo)
            JuMP.set_upper_bound(v, hi)
            JuMP.set_start_value(v, start_val)
        end
    end

    return nothing
end
