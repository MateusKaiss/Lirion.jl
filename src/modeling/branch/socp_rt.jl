@inline function _rt_product_bounds(bus::BusData, br::BranchData)
    angmin = br.angmin
    angmax = br.angmax
    bf = br.from
    bt = br.to
    L = br.n

    R_min = fill(-Inf, L)
    R_max = fill(Inf, L)
    T_min = fill(-Inf, L)
    T_max = fill(Inf, L)

    @inbounds for l in 1:L
        bus_i = bf[l]
        bus_j = bt[l]
        idx_i = bus.id_to_idx[bus_i]
        idx_j = bus.id_to_idx[bus_j]

        vm_fr_min = bus.vmin[idx_i]
        vm_fr_max = bus.vmax[idx_i]
        vm_to_min = bus.vmin[idx_j]
        vm_to_max = bus.vmax[idx_j]

        amin = angmin[l]
        amax = angmax[l]

        if amin >= 0
            R_max[l] = min(R_max[l], vm_fr_max * vm_to_max * cos(amin))
            R_min[l] = max(R_min[l], vm_fr_min * vm_to_min * cos(amax))
            T_max[l] = min(T_max[l], vm_fr_max * vm_to_max * sin(amax))
            T_min[l] = max(T_min[l], vm_fr_min * vm_to_min * sin(amin))
        end
        if amax <= 0
            R_max[l] = min(R_max[l], vm_fr_max * vm_to_max * cos(amax))
            R_min[l] = max(R_min[l], vm_fr_min * vm_to_min * cos(amin))
            T_max[l] = min(T_max[l], vm_fr_min * vm_to_max * sin(amax))
            T_min[l] = max(T_min[l], vm_fr_max * vm_to_max * sin(amin))
        end
        if amin < 0 && amax > 0
            R_max[l] = min(R_max[l], vm_fr_max * vm_to_max * 1.0)
            R_min[l] = max(R_min[l], vm_fr_min * vm_to_min * min(cos(amin), cos(amax)))
            T_max[l] = min(T_max[l], vm_fr_max * vm_to_max * sin(amax))
            T_min[l] = max(T_min[l], vm_fr_max * vm_to_max * sin(amin))
        end
    end

    return R_min, R_max, T_min, T_max
end

function create_rt_variables_with_bounds!(
    model::JuMP.Model,
    bus::BusData,
    br::BranchData;
    bounded::Bool = true,
)
    T = get(model.ext, :time_periods, 1)
    active = br.active

    @variable(model, R[l in active, t in 1:T], base_name = "R")
    @variable(model, T_var[l in active, t in 1:T], base_name = "T_var")

    vm_start = bus.vm_start
    va_start = bus.va_start
    shf = br.tap_angle

    @inbounds for l in active, t in 1:T
        bus_i = br.from[l]
        bus_j = br.to[l]

        idx_i = bus.id_to_idx[bus_i]
        idx_j = bus.id_to_idx[bus_j]

        Vm_i = vm_start[idx_i]
        Vm_j = vm_start[idx_j]
        Va_i = va_start[idx_i]
        Va_j = va_start[idx_j]

        if !(
            Vm_i === nothing || Vm_j === nothing ||
            Va_i === nothing || Va_j === nothing
        ) && !any(isnan.((Vm_i, Vm_j, Va_i, Va_j, shf[l])))
            set_start_value(R[l, t], Vm_i * Vm_j * cos(Va_i - Va_j - shf[l]))
            set_start_value(T_var[l, t], Vm_i * Vm_j * sin(Va_i - Va_j - shf[l]))
        end
    end

    if bounded
        R_min, R_max, T_min, T_max = _rt_product_bounds(bus, br)
        @inbounds for l in active, t in 1:T
            isfinite(R_min[l]) && set_lower_bound(R[l, t], R_min[l])
            isfinite(R_max[l]) && set_upper_bound(R[l, t], R_max[l])
            isfinite(T_min[l]) && set_lower_bound(T_var[l, t], T_min[l])
            isfinite(T_max[l]) && set_upper_bound(T_var[l, t], T_max[l])
        end
    end
    return nothing
end

@inline function _tan_clip(x; eps = 1e-6)
    b = (pi / 2) - eps
    x = min(max(x, -b), b)
    return tan(x)
end

function add_rt_tangent_angle_cuts!(model::JuMP.Model, br::BranchData)
    T = get(model.ext, :time_periods, 1)
    active = br.active

    R = model[:R]
    T_var = model[:T_var]

    @constraint(model, [l in active, t in 1:T; !isnan(br.angmax[l])],
        T_var[l, t] <= _tan_clip(br.angmax[l]) * R[l, t]
    )
    @constraint(model, [l in active, t in 1:T; !isnan(br.angmin[l])],
        T_var[l, t] >= _tan_clip(br.angmin[l]) * R[l, t]
    )
    return nothing
end
