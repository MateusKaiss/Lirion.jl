# ============================================================
# AC (polar V, θ)
# ============================================================

function add_ac_branch_expressions!(model::JuMP.Model, br::BranchData)
    T = get(model.ext, :time_periods, 1)
    active = br.active

    theta = model[:theta]
    v_magnitude = model[:v_magnitude]

    @expression(model, _ac_dtheta[l in active, t in 1:T],
        theta[br.from[l], t] - theta[br.to[l], t] - br.tap_angle[l],
    )

    @expression(model, _ac_cos_dtheta[l in active, t in 1:T], cos(_ac_dtheta[l, t]))
    @expression(model, _ac_sin_dtheta[l in active, t in 1:T], sin(_ac_dtheta[l, t]))

    @expression(model, p_flow_branches[l in active, t in 1:T],
        (br.g[l] + br.g_fr[l]) * (br.inv_tap[l] * v_magnitude[br.from[l], t])^2 -
        br.inv_tap[l] * v_magnitude[br.from[l], t] * v_magnitude[br.to[l], t] *
        (br.g[l] * _ac_cos_dtheta[l, t] + br.b[l] * _ac_sin_dtheta[l, t])
    )

    @expression(model, p_reverse_flow_branches[l in active, t in 1:T],
        (br.g[l] + br.g_to[l]) * v_magnitude[br.to[l], t]^2 -
        br.inv_tap[l] * v_magnitude[br.to[l], t] * v_magnitude[br.from[l], t] *
        (br.g[l] * _ac_cos_dtheta[l, t] - br.b[l] * _ac_sin_dtheta[l, t])
    )

    @expression(model, q_flow_branches[l in active, t in 1:T],
        -(br.b[l] + br.b_fr[l]) *
        (br.inv_tap[l] * v_magnitude[br.from[l], t])^2 -
        br.inv_tap[l] * v_magnitude[br.from[l], t] * v_magnitude[br.to[l], t] *
        (-br.b[l] * _ac_cos_dtheta[l, t] + br.g[l] * _ac_sin_dtheta[l, t])
    )

    @expression(model, q_reverse_flow_branches[l in active, t in 1:T],
        -(br.b[l] + br.b_to[l]) * v_magnitude[br.to[l], t]^2 -
        br.inv_tap[l] * v_magnitude[br.to[l], t] * v_magnitude[br.from[l], t] *
        (-br.b[l] * _ac_cos_dtheta[l, t] - br.g[l] * _ac_sin_dtheta[l, t])
    )

    return nothing
end

function add_ac_branch_constraints!(model::JuMP.Model, br::BranchData)
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
        branch_angle_lb[i in active, t in 1:T; !isnan(br.ang_min[i])],
        br.ang_min[i] <=
        model[:theta][br.from[i], t] - model[:theta][br.to[i], t]
    )

    @constraint(model,
        branch_angle_ub[i in active, t in 1:T; !isnan(br.ang_max[i])],
        model[:theta][br.from[i], t] - model[:theta][br.to[i], t] <=
        br.ang_max[i]
    )

    return nothing
end
