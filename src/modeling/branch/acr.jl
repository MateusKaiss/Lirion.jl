function add_acr_branch_expressions!(model::JuMP.Model, br::BranchData)
    T = get(model.ext, :time_periods, 1)
    active = br.active

    vr = model[:v_real]
    vi = model[:v_imag]

    @expression(model, _acr_v_sq_fr[l in active, t in 1:T],
        vr[br.from[l], t]^2 + vi[br.from[l], t]^2,
    )
    @expression(model, _acr_v_sq_to[l in active, t in 1:T],
        vr[br.to[l], t]^2 + vi[br.to[l], t]^2,
    )
    @expression(model, _acr_cos_mix[l in active, t in 1:T],
        vr[br.from[l], t] * vr[br.to[l], t] +
        vi[br.from[l], t] * vi[br.to[l], t],
    )
    @expression(model, _acr_sin_mix[l in active, t in 1:T],
        vi[br.from[l], t] * vr[br.to[l], t] -
        vr[br.from[l], t] * vi[br.to[l], t],
    )

    @expression(model, p_flow_branches[l in active, t in 1:T],
        br.g[l] * br.inv_tap[l]^2 * _acr_v_sq_fr[l, t] +
        (-br.g[l] * br.tap_real[l] + br.b[l] * br.tap_imag[l]) *
        br.inv_tap[l]^2 * _acr_cos_mix[l, t] +
        (-br.b[l] * br.tap_real[l] - br.g[l] * br.tap_imag[l]) *
        br.inv_tap[l]^2 * _acr_sin_mix[l, t],
    )

    @expression(model, q_flow_branches[l in active, t in 1:T],
        -(br.b[l] + br.b_sh[l] / 2) * br.inv_tap[l]^2 * _acr_v_sq_fr[l, t] -
        (-br.b[l] * br.tap_real[l] - br.g[l] * br.tap_imag[l]) *
        br.inv_tap[l]^2 * _acr_cos_mix[l, t] +
        (-br.g[l] * br.tap_real[l] + br.b[l] * br.tap_imag[l]) *
        br.inv_tap[l]^2 * _acr_sin_mix[l, t],
    )

    @expression(model, p_reverse_flow_branches[l in active, t in 1:T],
        br.g[l] * _acr_v_sq_to[l, t] +
        (-br.g[l] * br.tap_real[l] - br.b[l] * br.tap_imag[l]) *
        br.inv_tap[l]^2 * _acr_cos_mix[l, t] -
        (-br.b[l] * br.tap_real[l] + br.g[l] * br.tap_imag[l]) *
        br.inv_tap[l]^2 * _acr_sin_mix[l, t],
    )

    @expression(model, q_reverse_flow_branches[l in active, t in 1:T],
        -(br.b[l] + br.b_sh[l] / 2) * _acr_v_sq_to[l, t] -
        (-br.b[l] * br.tap_real[l] + br.g[l] * br.tap_imag[l]) *
        br.inv_tap[l]^2 * _acr_cos_mix[l, t] -
        (-br.g[l] * br.tap_real[l] - br.b[l] * br.tap_imag[l]) *
        br.inv_tap[l]^2 * _acr_sin_mix[l, t],
    )

    return nothing
end

function add_acr_branch_constraints!(model::JuMP.Model, br::BranchData)
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
        branch_angle_lb[l in active, t in 1:T; !isnan(br.ang_min[l])],
        _tan_clip(br.ang_min[l]) * model[:_acr_cos_mix][l, t] <=
        model[:_acr_sin_mix][l, t]
    )

    @constraint(model,
        branch_angle_ub[l in active, t in 1:T; !isnan(br.ang_max[l])],
        model[:_acr_sin_mix][l, t] <=
        _tan_clip(br.ang_max[l]) * model[:_acr_cos_mix][l, t]
    )

    return nothing
end
