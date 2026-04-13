function add_socwr_branch_expressions!(model::JuMP.Model, br::BranchData)
    T = get(model.ext, :time_periods, 1)
    active = br.active

    inv_tap = br.inv_tap
    R = model[:R]
    Tv = model[:T_var]

    tap_real = br.tap_real
    tap_imag = br.tap_imag

    @expression(model, p_flow_branches[l in active, t in 1:T],
        br.g[l] * inv_tap[l]^2 * model[:v_magnitude_sq][br.from[l], t] +
        (-br.g[l] * tap_real[l] + br.b[l] * tap_imag[l]) * inv_tap[l]^2 * R[l, t] +
        (-br.b[l] * tap_real[l] - br.g[l] * tap_imag[l]) * inv_tap[l]^2 * Tv[l, t]
    )

    @expression(model, q_flow_branches[l in active, t in 1:T],
        -(br.b[l] + br.b_sh[l] / 2) * inv_tap[l]^2 *
        model[:v_magnitude_sq][br.from[l], t] -
        (-br.b[l] * tap_real[l] - br.g[l] * tap_imag[l]) * inv_tap[l]^2 * R[l, t] +
        (-br.g[l] * tap_real[l] + br.b[l] * tap_imag[l]) * inv_tap[l]^2 * Tv[l, t]
    )

    @expression(model, p_reverse_flow_branches[l in active, t in 1:T],
        br.g[l] * model[:v_magnitude_sq][br.to[l], t] +
        (-br.g[l] * tap_real[l] - br.b[l] * tap_imag[l]) * inv_tap[l]^2 * R[l, t] +
        (-br.b[l] * tap_real[l] + br.g[l] * tap_imag[l]) * inv_tap[l]^2 * (-Tv[l, t])
    )

    @expression(model, q_reverse_flow_branches[l in active, t in 1:T],
        -(br.b[l] + br.b_sh[l] / 2) *
        model[:v_magnitude_sq][br.to[l], t] -
        (-br.b[l] * tap_real[l] + br.g[l] * tap_imag[l]) * inv_tap[l]^2 * R[l, t] +
        (-br.g[l] * tap_real[l] - br.b[l] * tap_imag[l]) * inv_tap[l]^2 * (-Tv[l, t])
    )

    return nothing
end

function add_socwr_branch_constraints!(model::JuMP.Model, br::BranchData)
    T = get(model.ext, :time_periods, 1)
    active = br.active

    @constraint(model, branch_angle_lb[i in active, t in 1:T; !isnan(br.ang_min[i])],
        br.ang_min[i] <=
        model[:theta][br.from[i], t] - model[:theta][br.to[i], t]
    )

    @constraint(model, branch_angle_ub[i in active, t in 1:T; !isnan(br.ang_max[i])],
        model[:theta][br.from[i], t] - model[:theta][br.to[i], t] <=
        br.ang_max[i]
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
