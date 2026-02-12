function create_dc_branch_expressions!(model::JuMP.Model, br::BranchData)
    T = get(model.ext, :time_periods, 1)
    active = br.active

    @expression(model, p_flow_branches[i in active, t in 1:T],
        br.b[i] * (
            model[:theta][br.from[i], t] -
            model[:theta][br.to[i], t]
        )
    )

    @expression(model, p_reverse_flow_branches[i in active, t in 1:T],
        -model[:p_flow_branches][i, t]
    )

    @expression(model, losses_branches[i in active, t in 1:T],
        model[:p_flow_branches][i, t]^2 * br.r[i]
    )
    return nothing
end

function create_dc_branch_constraints!(model::JuMP.Model, br::BranchData)
    T = get(model.ext, :time_periods, 1)
    active = br.active

    @constraint(model, capacity[i in active, t in 1:T],
        -br.rate_a[i] <= model[:p_flow_branches][i, t] <= br.rate_a[i]
    )

    @constraint(model,
        branch_angle_lb[i in active, t in 1:T; !isnan(br.angmin[i])],
        br.angmin[i] <=
        model[:theta][br.from[i], t] - model[:theta][br.to[i], t]
    )

    @constraint(model,
        branch_angle_ub[i in active, t in 1:T; !isnan(br.angmax[i])],
        model[:theta][br.from[i], t] - model[:theta][br.to[i], t] <=
        br.angmax[i]
    )

    return nothing
end

function create_ac_branch_expressions!(model::JuMP.Model, br::BranchData)
    T = get(model.ext, :time_periods, 1)
    active = br.active

    @expression(model, p_flow_branches[l in active, t in 1:T],
        br.g[l] * (br.inv_tap[l] * model[:v_magnitude][br.from[l], t])^2 -
        br.inv_tap[l] * model[:v_magnitude][br.from[l], t] *
        model[:v_magnitude][br.to[l], t] *
        (
            br.g[l] * cos(
                model[:theta][br.from[l], t] -
                model[:theta][br.to[l], t] - br.tap_angle[l],
            ) +
            br.b[l] * sin(
                model[:theta][br.from[l], t] -
                model[:theta][br.to[l], t] - br.tap_angle[l],
            )
        )
    )

    @expression(model, p_reverse_flow_branches[l in active, t in 1:T],
        br.g[l] * model[:v_magnitude][br.to[l], t]^2 -
        br.inv_tap[l] * model[:v_magnitude][br.to[l], t] *
        model[:v_magnitude][br.from[l], t] *
        (
            br.g[l] * cos(
                model[:theta][br.from[l], t] -
                model[:theta][br.to[l], t] - br.tap_angle[l],
            ) -
            br.b[l] * sin(
                model[:theta][br.from[l], t] -
                model[:theta][br.to[l], t] - br.tap_angle[l],
            )
        )
    )

    @expression(model, q_flow_branches[l in active, t in 1:T],
        -(br.b[l] + br.b_sh[l] / 2) *
        (br.inv_tap[l] * model[:v_magnitude][br.from[l], t])^2 -
        br.inv_tap[l] * model[:v_magnitude][br.from[l], t] *
        model[:v_magnitude][br.to[l], t] *
        (
            -br.b[l] * cos(
                model[:theta][br.from[l], t] -
                model[:theta][br.to[l], t] - br.tap_angle[l],
            ) +
            br.g[l] * sin(
                model[:theta][br.from[l], t] -
                model[:theta][br.to[l], t] - br.tap_angle[l],
            )
        )
    )

    @expression(model, q_reverse_flow_branches[l in active, t in 1:T],
        -(br.b[l] + br.b_sh[l] / 2) *
        model[:v_magnitude][br.to[l], t]^2 -
        br.inv_tap[l] * model[:v_magnitude][br.to[l], t] *
        model[:v_magnitude][br.from[l], t] *
        (
            -br.b[l] * cos(
                model[:theta][br.from[l], t] -
                model[:theta][br.to[l], t] - br.tap_angle[l],
            ) -
            br.g[l] * sin(
                model[:theta][br.from[l], t] -
                model[:theta][br.to[l], t] - br.tap_angle[l],
            )
        )
    )
    return nothing
end

function create_ac_branch_constraints!(model::JuMP.Model, br::BranchData)
    T = get(model.ext, :time_periods, 1)
    active = br.active

    @constraint(model, [l in active, t in 1:T],
        model[:p_flow_branches][l, t]^2 +
        model[:q_flow_branches][l, t]^2 <= br.rate_a[l]^2
    )
    @constraint(model, [l in active, t in 1:T],
        model[:p_reverse_flow_branches][l, t]^2 +
        model[:q_reverse_flow_branches][l, t]^2 <= br.rate_a[l]^2
    )

    @constraint(model,
        branch_angle_lb[i in active, t in 1:T; !isnan(br.angmin[i])],
        br.angmin[i] <=
        model[:theta][br.from[i], t] - model[:theta][br.to[i], t]
    )

    @constraint(model,
        branch_angle_ub[i in active, t in 1:T; !isnan(br.angmax[i])],
        model[:theta][br.from[i], t] - model[:theta][br.to[i], t] <=
        br.angmax[i]
    )

    return nothing
end

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

    vm_start = bus.vstart
    va_start = bus.vastart
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

function create_socp_branch_expressions!(model::JuMP.Model, br::BranchData)
    T = get(model.ext, :time_periods, 1)
    active = br.active

    inv_tap = br.inv_tap
    R = model[:R]
    Tv = model[:T_var]

    tr = br.tr
    ti = br.ti

    @expression(model, p_flow_branches[l in active, t in 1:T],
        br.g[l] * inv_tap[l]^2 * model[:v_magnitude_sq][br.from[l], t] +
        (-br.g[l] * tr[l] + br.b[l] * ti[l]) * inv_tap[l]^2 * R[l, t] +
        (-br.b[l] * tr[l] - br.g[l] * ti[l]) * inv_tap[l]^2 * Tv[l, t]
    )
    @expression(model, q_flow_branches[l in active, t in 1:T],
        -(br.b[l] + br.b_sh[l] / 2) * inv_tap[l]^2 *
        model[:v_magnitude_sq][br.from[l], t] -
        (-br.b[l] * tr[l] - br.g[l] * ti[l]) * inv_tap[l]^2 * R[l, t] +
        (-br.g[l] * tr[l] + br.b[l] * ti[l]) * inv_tap[l]^2 * Tv[l, t]
    )

    @expression(model, p_reverse_flow_branches[l in active, t in 1:T],
        br.g[l] * model[:v_magnitude_sq][br.to[l], t] +
        (-br.g[l] * tr[l] - br.b[l] * ti[l]) * inv_tap[l]^2 * R[l, t] +
        (-br.b[l] * tr[l] + br.g[l] * ti[l]) * inv_tap[l]^2 * (-Tv[l, t])
    )
    @expression(model, q_reverse_flow_branches[l in active, t in 1:T],
        -(br.b[l] + br.b_sh[l] / 2) *
        model[:v_magnitude_sq][br.to[l], t] -
        (-br.b[l] * tr[l] + br.g[l] * ti[l]) * inv_tap[l]^2 * R[l, t] +
        (-br.g[l] * tr[l] - br.b[l] * ti[l]) * inv_tap[l]^2 * (-Tv[l, t])
    )

    return nothing
end

function create_socp_branch_constraints!(model::JuMP.Model, br::BranchData)
    T = get(model.ext, :time_periods, 1)
    active = br.active

    @constraint(model, branch_angle_lb[i in active, t in 1:T; !isnan(br.angmin[i])],
        br.angmin[i] <=
        model[:theta][br.from[i], t] - model[:theta][br.to[i], t]
    )
    @constraint(model, branch_angle_ub[i in active, t in 1:T; !isnan(br.angmax[i])],
        model[:theta][br.from[i], t] - model[:theta][br.to[i], t] <=
        br.angmax[i]
    )

    @constraint(model, [l in active, t in 1:T],
        [
            model[:v_magnitude_sq][br.from[l], t] / sqrt(2),
            model[:v_magnitude_sq][br.to[l], t] / sqrt(2),
            model[:R][l, t],
            model[:T_var][l, t],
        ] in JuMP.RotatedSecondOrderCone()
    )

    @constraint(model, [l in active, t in 1:T],
        [
            br.rate_a[l],
            model[:p_flow_branches][l, t],
            model[:q_flow_branches][l, t],
        ] in JuMP.SecondOrderCone()
    )
    @constraint(model, [l in active, t in 1:T],
        [
            br.rate_a[l],
            model[:p_reverse_flow_branches][l, t],
            model[:q_reverse_flow_branches][l, t],
        ] in JuMP.SecondOrderCone()
    )

    add_rt_tangent_angle_cuts!(model, br)

    return nothing
end

# =======================================================================
# Dispatch by model type (DC / AC / SOCP)
# =======================================================================
function branches!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{DC},
)
    br = sys.branch
    create_dc_branch_expressions!(model, br)
    create_dc_branch_constraints!(model, br)
    return nothing
end

function branches!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{AC},
)
    br = sys.branch
    create_ac_branch_expressions!(model, br)
    create_ac_branch_constraints!(model, br)
    return nothing
end

function branches!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{SOCP},
)
    bus = sys.bus
    br = sys.branch

    create_rt_variables_with_bounds!(model, bus, br)
    create_socp_branch_expressions!(model, br)
    create_socp_branch_constraints!(model, br)

    return nothing
end
