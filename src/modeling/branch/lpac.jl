# ============================================================
# Variables
# ============================================================

function add_lpac_branch_variables!(model::JuMP.Model, br::BranchData)
    T = get(model.ext, :time_periods, 1)
    active = br.active
    theta = model[:theta]

    # cs variable per branch — cosine of angle difference.
    @variable(model, cs[l in active, t in 1:T], start = 1.0)

    @inbounds for l in active
        angmin = br.ang_min[l]
        angmax = br.ang_max[l]

        # Compute bounds on cos(Δθ) from angle limits.
        if !isnan(angmin) && angmin >= 0
            cs_lo = cos(angmax)
            cs_hi = cos(angmin)
        elseif !isnan(angmax) && angmax <= 0
            cs_lo = cos(angmin)
            cs_hi = cos(angmax)
        elseif !isnan(angmin) && !isnan(angmax) && angmin < 0 && angmax > 0
            cs_lo = min(cos(angmin), cos(angmax))
            cs_hi = 1.0
        else
            cs_lo = -1.0
            cs_hi = 1.0
        end

        for t in 1:T
            JuMP.set_lower_bound(cs[l, t], cs_lo)
            JuMP.set_upper_bound(cs[l, t], cs_hi)
        end

        # Quadratic upper bound: cs ≤ 1 − (1−cos(vad_max))/vad_max² · Δθ²
        vad_max = max(
            isnan(angmin) ? 0.0 : abs(angmin),
            isnan(angmax) ? 0.0 : abs(angmax),
        )
        if vad_max > 0
            c_coeff = (1.0 - cos(vad_max)) / vad_max^2
            for t in 1:T
                dtheta = theta[br.from[l], t] - theta[br.to[l], t]
                @constraint(model, cs[l, t] <= 1.0 - c_coeff * dtheta^2)
            end
        end
    end

    return nothing
end

# ============================================================
# Expressions
# ============================================================

function add_lpac_branch_expressions!(model::JuMP.Model, br::BranchData)
    T = get(model.ext, :time_periods, 1)
    active = br.active

    theta = model[:theta]
    phi = model[:phi]
    cs = model[:cs]

    # From-side P: (g + g_fr)/tm² · (1+2ϕ_fr) + (-g·tr+b·ti)/tm² · (cs+ϕ_fr+ϕ_to) + (-b·tr-g·ti)/tm² · Δθ
    @expression(model, p_flow_branches[l in active, t in 1:T],
        (br.g[l] + br.g_fr[l]) * br.inv_tap[l]^2 * (1.0 + 2.0 * phi[br.from[l], t]) +
        (-br.g[l] * br.tap_real[l] + br.b[l] * br.tap_imag[l]) *
        br.inv_tap[l]^2 * (cs[l, t] + phi[br.from[l], t] + phi[br.to[l], t]) +
        (-br.b[l] * br.tap_real[l] - br.g[l] * br.tap_imag[l]) *
        br.inv_tap[l]^2 * (theta[br.from[l], t] - theta[br.to[l], t]),
    )

    # From-side Q: -(b + b_fr)/tm² · (1+2ϕ_fr) - (-b·tr-g·ti)/tm² · (cs+ϕ_fr+ϕ_to) + (-g·tr+b·ti)/tm² · Δθ
    @expression(model, q_flow_branches[l in active, t in 1:T],
        -(br.b[l] + br.b_fr[l]) * br.inv_tap[l]^2 *
        (1.0 + 2.0 * phi[br.from[l], t]) -
        (-br.b[l] * br.tap_real[l] - br.g[l] * br.tap_imag[l]) *
        br.inv_tap[l]^2 * (cs[l, t] + phi[br.from[l], t] + phi[br.to[l], t]) +
        (-br.g[l] * br.tap_real[l] + br.b[l] * br.tap_imag[l]) *
        br.inv_tap[l]^2 * (theta[br.from[l], t] - theta[br.to[l], t]),
    )

    # To-side P: (g + g_to) · (1+2ϕ_to) + (-g·tr-b·ti)/tm² · (cs+ϕ_fr+ϕ_to) - (-b·tr+g·ti)/tm² · Δθ
    @expression(model, p_reverse_flow_branches[l in active, t in 1:T],
        (br.g[l] + br.g_to[l]) * (1.0 + 2.0 * phi[br.to[l], t]) +
        (-br.g[l] * br.tap_real[l] - br.b[l] * br.tap_imag[l]) *
        br.inv_tap[l]^2 * (cs[l, t] + phi[br.from[l], t] + phi[br.to[l], t]) -
        (-br.b[l] * br.tap_real[l] + br.g[l] * br.tap_imag[l]) *
        br.inv_tap[l]^2 * (theta[br.from[l], t] - theta[br.to[l], t]),
    )

    # To-side Q: -(b + b_to) · (1+2ϕ_to) - (-b·tr+g·ti)/tm² · (cs+ϕ_fr+ϕ_to) - (-g·tr-b·ti)/tm² · Δθ
    @expression(model, q_reverse_flow_branches[l in active, t in 1:T],
        -(br.b[l] + br.b_to[l]) * (1.0 + 2.0 * phi[br.to[l], t]) -
        (-br.b[l] * br.tap_real[l] + br.g[l] * br.tap_imag[l]) *
        br.inv_tap[l]^2 * (cs[l, t] + phi[br.from[l], t] + phi[br.to[l], t]) -
        (-br.g[l] * br.tap_real[l] - br.b[l] * br.tap_imag[l]) *
        br.inv_tap[l]^2 * (theta[br.from[l], t] - theta[br.to[l], t]),
    )

    return nothing
end

# ============================================================
# Constraints
# ============================================================

function add_lpac_branch_constraints!(model::JuMP.Model, br::BranchData)
    T = get(model.ext, :time_periods, 1)
    active = br.active

    # Explicit flow variables for thermal limits (strongly convex in 2 vars).
    @variable(model, _lpac_p_fr[l in active, t in 1:T])
    @variable(model, _lpac_q_fr[l in active, t in 1:T])
    @variable(model, _lpac_p_to[l in active, t in 1:T])
    @variable(model, _lpac_q_to[l in active, t in 1:T])

    @inbounds for l in active
        ra = br.rate_a[l]
        for t in 1:T
            @constraint(model, _lpac_p_fr[l, t] == model[:p_flow_branches][l, t])
            @constraint(model, _lpac_q_fr[l, t] == model[:q_flow_branches][l, t])
            @constraint(model, _lpac_p_to[l, t] == model[:p_reverse_flow_branches][l, t])
            @constraint(model, _lpac_q_to[l, t] == model[:q_reverse_flow_branches][l, t])

            JuMP.set_lower_bound(_lpac_p_fr[l, t], -ra)
            JuMP.set_upper_bound(_lpac_p_fr[l, t], ra)
            JuMP.set_lower_bound(_lpac_q_fr[l, t], -ra)
            JuMP.set_upper_bound(_lpac_q_fr[l, t], ra)
            JuMP.set_lower_bound(_lpac_p_to[l, t], -ra)
            JuMP.set_upper_bound(_lpac_p_to[l, t], ra)
            JuMP.set_lower_bound(_lpac_q_to[l, t], -ra)
            JuMP.set_upper_bound(_lpac_q_to[l, t], ra)
        end
    end

    @constraint(model, [l in active, t in 1:T],
        _lpac_p_fr[l, t]^2 + _lpac_q_fr[l, t]^2 <= br.rate_a[l]^2
    )
    @constraint(model, [l in active, t in 1:T],
        _lpac_p_to[l, t]^2 + _lpac_q_to[l, t]^2 <= br.rate_a[l]^2
    )

    @constraint(model,
        branch_angle_lb[l in active, t in 1:T; !isnan(br.ang_min[l])],
        br.ang_min[l] <=
        model[:theta][br.from[l], t] - model[:theta][br.to[l], t]
    )

    @constraint(model,
        branch_angle_ub[l in active, t in 1:T; !isnan(br.ang_max[l])],
        model[:theta][br.from[l], t] - model[:theta][br.to[l], t] <=
        br.ang_max[l]
    )

    return nothing
end
