
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
