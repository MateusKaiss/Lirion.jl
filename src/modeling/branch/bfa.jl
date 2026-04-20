# =======================================================================
# BFA (Branch Flow Approximation — LinDistFlow)
# =======================================================================

function add_bfa_branch_variables!(model::JuMP.Model, bus::BusData, br::BranchData)
    T = get(model.ext, :time_periods, 1)
    active = br.active

    @variable(model, p_flow_branches[l in active, t in 1:T])
    @variable(model, q_flow_branches[l in active, t in 1:T])
    @variable(model, p_reverse_flow_branches[l in active, t in 1:T])
    @variable(model, q_reverse_flow_branches[l in active, t in 1:T])

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
    end

    return nothing
end

function add_bfa_branch_constraints!(model::JuMP.Model, bus::BusData, br::BranchData)
    T = get(model.ext, :time_periods, 1)
    active = br.active

    w = model[:v_magnitude_sq]
    p_fr = model[:p_flow_branches]
    q_fr = model[:q_flow_branches]
    p_to = model[:p_reverse_flow_branches]
    q_to = model[:q_reverse_flow_branches]

    r = br.r
    x = br.x
    g_fr = br.g_fr
    b_fr = br.b_fr
    g_to = br.g_to
    b_to = br.b_to
    inv_t = br.inv_tap

    # --- Active power losses (shunt only) ---
    @constraint(model, bfa_p_loss[l in active, t in 1:T],
        p_fr[l, t] + p_to[l, t] ==
        g_fr[l] * inv_t[l]^2 * w[br.from[l], t] +
        g_to[l] * w[br.to[l], t]
    )

    # --- Reactive power losses (shunt only) ---
    @constraint(model, bfa_q_loss[l in active, t in 1:T],
        q_fr[l, t] + q_to[l, t] ==
        -b_fr[l] * inv_t[l]^2 * w[br.from[l], t] -
        b_to[l] * w[br.to[l], t]
    )

    # --- Linear voltage magnitude difference ---
    @constraint(model, bfa_voltage_drop[l in active, t in 1:T],
        inv_t[l]^2 * w[br.from[l], t] - w[br.to[l], t] ==
        2 * (r[l] * p_fr[l, t] + x[l] * q_fr[l, t])
    )


    # --- Angle difference cuts (BF form, reused from SOCBF) ---
    _add_socbf_angle_constraints!(model, br)

    return nothing
end
