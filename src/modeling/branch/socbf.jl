# =======================================================================
# SOC-BF (Branch Flow SOC relaxation
# =======================================================================

function add_socbf_branch_variables!(model::JuMP.Model, bus::BusData, br::BranchData)
    T = get(model.ext, :time_periods, 1)
    active = br.active
    id_to_idx = bus.id_to_idx

    @variable(model, p_flow_branches[l in active, t in 1:T])
    @variable(model, q_flow_branches[l in active, t in 1:T])
    @variable(model, p_reverse_flow_branches[l in active, t in 1:T])
    @variable(model, q_reverse_flow_branches[l in active, t in 1:T])
    @variable(model, ccm[l in active, t in 1:T], lower_bound = 0.0)

    @inbounds for l in active, t in 1:T
        ra = br.rate_a[l]

        set_lower_bound(p_flow_branches[l, t], -ra)
        set_upper_bound(p_flow_branches[l, t], ra)
        set_lower_bound(q_flow_branches[l, t], -ra)
        set_upper_bound(q_flow_branches[l, t], ra)
        set_lower_bound(p_reverse_flow_branches[l, t], -ra)
        set_upper_bound(p_reverse_flow_branches[l, t], ra)
        set_lower_bound(q_reverse_flow_branches[l, t], -ra)
        set_upper_bound(q_reverse_flow_branches[l, t], ra)

        # ccm upper bound: (rate_a * tap / v_min_from)^2
        idx_fr = id_to_idx[br.from[l]]
        v_min_fr = bus.v_min[idx_fr]
        ccm_ub = (ra * br.tap_position[l] / v_min_fr)^2
        if isfinite(ccm_ub)
            set_upper_bound(ccm[l, t], ccm_ub)
        end
    end

    return nothing
end

function add_socbf_branch_constraints!(model::JuMP.Model, bus::BusData, br::BranchData)
    T = get(model.ext, :time_periods, 1)
    active = br.active

    w = model[:v_magnitude_sq]
    p_fr = model[:p_flow_branches]
    q_fr = model[:q_flow_branches]
    p_to = model[:p_reverse_flow_branches]
    q_to = model[:q_reverse_flow_branches]
    ccm = model[:ccm]

    r = br.r
    x = br.x
    g_fr = br.g_fr
    b_fr = br.b_fr
    g_to = br.g_to
    b_to = br.b_to
    inv_t = br.inv_tap
    tp = br.tap_real
    ti = br.tap_imag

    # --- Active power losses ---
    @constraint(model, socbf_p_loss[l in active, t in 1:T],
        p_fr[l, t] + p_to[l, t] ==
        r[l] * (
            ccm[l, t] +
            (g_fr[l]^2 + b_fr[l]^2) * inv_t[l]^2 * w[br.from[l], t] -
            2 * (g_fr[l] * p_fr[l, t] - b_fr[l] * q_fr[l, t])
        ) +
        g_fr[l] * inv_t[l]^2 * w[br.from[l], t] +
        g_to[l] * w[br.to[l], t]
    )

    # --- Reactive power losses ---
    @constraint(model, socbf_q_loss[l in active, t in 1:T],
        q_fr[l, t] + q_to[l, t] ==
        x[l] * (
            ccm[l, t] +
            (g_fr[l]^2 + b_fr[l]^2) * inv_t[l]^2 * w[br.from[l], t] -
            2 * (g_fr[l] * p_fr[l, t] - b_fr[l] * q_fr[l, t])
        ) -
        b_fr[l] * inv_t[l]^2 * w[br.from[l], t] -
        b_to[l] * w[br.to[l], t]
    )

    # --- Voltage magnitude difference (Ohm's law) ---
    @constraint(model, socbf_voltage_drop[l in active, t in 1:T],
        (1 + 2 * (r[l] * g_fr[l] - x[l] * b_fr[l])) *
        inv_t[l]^2 * w[br.from[l], t] - w[br.to[l], t] ==
        2 * (r[l] * p_fr[l, t] + x[l] * q_fr[l, t]) -
        (r[l]^2 + x[l]^2) * (
            ccm[l, t] +
            (g_fr[l]^2 + b_fr[l]^2) * inv_t[l]^2 * w[br.from[l], t] -
            2 * (g_fr[l] * p_fr[l, t] - b_fr[l] * q_fr[l, t])
        )
    )

    # --- SOC relaxation: w_fr/tm² · ccm ≥ p_fr² + q_fr² ---
    @constraint(model, socbf_soc[l in active, t in 1:T],
        [
            inv_t[l]^2 * w[br.from[l], t],
            ccm[l, t] / 2,
            p_fr[l, t],
            q_fr[l, t],
        ] in JuMP.RotatedSecondOrderCone()
    )

    # --- Thermal limits (from end) ---
    @constraint(model, socbf_thermal_fr[l in active, t in 1:T],
        [
            br.rate_a[l],
            p_fr[l, t],
            q_fr[l, t],
        ] in JuMP.SecondOrderCone()
    )

    # --- Thermal limits (to end) ---
    @constraint(model, socbf_thermal_to[l in active, t in 1:T],
        [
            br.rate_a[l],
            p_to[l, t],
            q_to[l, t],
        ] in JuMP.SecondOrderCone()
    )

    # --- Angle difference constraints (BF form, no theta variable) ---
    _add_socbf_angle_constraints!(model, br)

    return nothing
end

function _add_socbf_angle_constraints!(model::JuMP.Model, br::BranchData)
    T = get(model.ext, :time_periods, 1)
    active = br.active

    w = model[:v_magnitude_sq]
    p_fr = model[:p_flow_branches]
    q_fr = model[:q_flow_branches]

    r = br.r
    x = br.x
    g_fr = br.g_fr
    b_fr = br.b_fr
    inv_t = br.inv_tap
    tr = br.tap_real
    ti = br.tap_imag

    @inbounds for l in active
        amin = br.ang_min[l]
        amax = br.ang_max[l]
        (isnan(amin) && isnan(amax)) && continue

        tzr = r[l] * tr[l] + x[l] * ti[l]
        tzi = r[l] * ti[l] - x[l] * tr[l]

        for t in 1:T
            w_fr_s = inv_t[l]^2 * w[br.from[l], t]

            # Re{(V_fr/tm)·conj(I_s·z)} components
            re_term =
                (tr[l] + tzr * g_fr[l] + tzi * b_fr[l]) * w_fr_s -
                tzr * p_fr[l, t] + tzi * q_fr[l, t]
            im_term =
                (ti[l] + tzi * g_fr[l] - tzr * b_fr[l]) * w_fr_s -
                tzi * p_fr[l, t] - tzr * q_fr[l, t]

            if !isnan(amin)
                @constraint(model, _tan_clip(amin) * re_term <= im_term)
            end
            if !isnan(amax)
                @constraint(model, _tan_clip(amax) * re_term >= im_term)
            end
        end
    end

    return nothing
end
