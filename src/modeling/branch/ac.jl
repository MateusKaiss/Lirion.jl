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
