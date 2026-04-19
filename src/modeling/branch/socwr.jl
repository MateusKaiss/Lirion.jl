function add_socwr_branch_expressions!(model::JuMP.Model, br::BranchData)
    T = get(model.ext, :time_periods, 1)
    active = br.active

    inv_tap = br.inv_tap
    R = model[:R]
    Tv = model[:T_var]

    tap_real = br.tap_real
    tap_imag = br.tap_imag

    @expression(model, p_flow_branches[l in active, t in 1:T],
        (br.g[l] + br.g_fr[l]) * inv_tap[l]^2 * model[:v_magnitude_sq][br.from[l], t] +
        (-br.g[l] * tap_real[l] + br.b[l] * tap_imag[l]) * inv_tap[l]^2 * R[l, t] +
        (-br.b[l] * tap_real[l] - br.g[l] * tap_imag[l]) * inv_tap[l]^2 * Tv[l, t]
    )

    @expression(model, q_flow_branches[l in active, t in 1:T],
        -(br.b[l] + br.b_fr[l]) * inv_tap[l]^2 *
        model[:v_magnitude_sq][br.from[l], t] -
        (-br.b[l] * tap_real[l] - br.g[l] * tap_imag[l]) * inv_tap[l]^2 * R[l, t] +
        (-br.g[l] * tap_real[l] + br.b[l] * tap_imag[l]) * inv_tap[l]^2 * Tv[l, t]
    )

    @expression(model, p_reverse_flow_branches[l in active, t in 1:T],
        (br.g[l] + br.g_to[l]) * model[:v_magnitude_sq][br.to[l], t] +
        (-br.g[l] * tap_real[l] - br.b[l] * tap_imag[l]) * inv_tap[l]^2 * R[l, t] +
        (-br.b[l] * tap_real[l] + br.g[l] * tap_imag[l]) * inv_tap[l]^2 * (-Tv[l, t])
    )

    @expression(model, q_reverse_flow_branches[l in active, t in 1:T],
        -(br.b[l] + br.b_to[l]) *
        model[:v_magnitude_sq][br.to[l], t] -
        (-br.b[l] * tap_real[l] + br.g[l] * tap_imag[l]) * inv_tap[l]^2 * R[l, t] +
        (-br.g[l] * tap_real[l] - br.b[l] * tap_imag[l]) * inv_tap[l]^2 * (-Tv[l, t])
    )

    return nothing
end

function add_socwr_branch_constraints!(
    model::JuMP.Model, bus::BusData, br::BranchData,
)
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

    # Lifted Nonlinear Cuts (Coffrin-Hijazi-Van Hentenryck 2015,
    # arXiv:1512.04644). Strengthen the W-space relaxation so Lirion
    # matches the PowerModels SOCWR tightness. Same cut derivation used
    # in QC and SDPWRM.
    add_socwr_lnc_cuts!(model, bus, br)

    return nothing
end

function add_socwr_lnc_cuts!(model::JuMP.Model, bus::BusData, br::BranchData)
    T = get(model.ext, :time_periods, 1)
    active = br.active
    id_to_idx = bus.id_to_idx
    v_min = bus.v_min
    v_max = bus.v_max

    for l in active
        amin = br.ang_min[l]
        amax = br.ang_max[l]
        (isnan(amin) || isnan(amax)) && continue
        amin >= -pi / 2 && amin <= pi / 2 || continue
        amax >= -pi / 2 && amax <= pi / 2 || continue
        amin < amax || continue

        idx_i = id_to_idx[br.from[l]]
        idx_j = id_to_idx[br.to[l]]

        vf_lb, vf_ub = v_min[idx_i], v_max[idx_i]
        vt_lb, vt_ub = v_min[idx_j], v_max[idx_j]

        phi = (amax + amin) / 2
        d = (amax - amin) / 2

        sf = vf_lb + vf_ub
        st = vt_lb + vt_ub

        for t in 1:T
            wf = model[:v_magnitude_sq][br.from[l], t]
            wt = model[:v_magnitude_sq][br.to[l], t]
            wr = model[:R][l, t]
            wi = model[:T_var][l, t]

            JuMP.@constraint(model,
                sf * st * (cos(phi) * wr + sin(phi) * wi) -
                vt_ub * cos(d) * st * wf -
                vf_ub * cos(d) * sf * wt >=
                vf_ub * vt_ub * cos(d) * (vf_lb * vt_lb - vf_ub * vt_ub)
            )
            JuMP.@constraint(model,
                sf * st * (cos(phi) * wr + sin(phi) * wi) -
                vt_lb * cos(d) * st * wf -
                vf_lb * cos(d) * sf * wt >=
                -vf_lb * vt_lb * cos(d) * (vf_lb * vt_lb - vf_ub * vt_ub)
            )
        end
    end

    return nothing
end
