function add_dc_branch_expressions!(model::JuMP.Model, br::BranchData)
    T = get(model.ext, :time_periods, 1)
    active = br.active

    @expression(model, p_flow_branches[i in active, t in 1:T],
        -br.b[i] * (
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

function add_dc_branch_constraints!(model::JuMP.Model, br::BranchData)
    T = get(model.ext, :time_periods, 1)
    active = br.active

    @constraint(model, capacity[i in active, t in 1:T],
        -br.rate_a[i] <= model[:p_flow_branches][i, t] <= br.rate_a[i]
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
