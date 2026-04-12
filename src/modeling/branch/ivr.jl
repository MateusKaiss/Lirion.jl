function add_ivr_branch_variables!(model::JuMP.Model, br::BranchData)
    T = get(model.ext, :time_periods, 1)
    active = br.active

    @variable(model, cr_fr[l in active, t in 1:T])
    @variable(model, ci_fr[l in active, t in 1:T])
    @variable(model, cr_to[l in active, t in 1:T])
    @variable(model, ci_to[l in active, t in 1:T])

    return nothing
end

function add_ivr_branch_constraints!(
    model::JuMP.Model,
    bus::BusData,
    br::BranchData,
)
    T = get(model.ext, :time_periods, 1)
    active = br.active

    vr = model[:v_real]
    vi = model[:v_imag]
    cr_fr = model[:cr_fr]
    ci_fr = model[:ci_fr]
    cr_to = model[:cr_to]
    ci_to = model[:ci_to]

    # Ohm's law: I_from = Y_ff·V_from + Y_ft·V_to  (linear)
    #            I_to   = Y_tf·V_from + Y_tt·V_to
    #
    # Y_ff = (g+g_fr + j(b+b_fr)) / tm²
    # Y_ft = -(g+jb)(tr+jti) / tm²  →  Y_ft_r = -(g·tr - b·ti)/tm²
    #                                    Y_ft_i = -(g·ti + b·tr)/tm²
    # Y_tf = -(g+jb)(tr-jti) / tm²  →  Y_tf_r = -(g·tr + b·ti)/tm²
    #                                    Y_tf_i = -(b·tr - g·ti)/tm² = (g·ti - b·tr)/tm²
    # Y_tt = g+g_to + j(b+b_to)

    @inbounds for l in active
        g = br.g[l]
        b_ln = br.b[l]
        g_fr_v = br.g_fr[l]
        b_fr_v = br.b_fr[l]
        g_to_v = br.g_to[l]
        b_to_v = br.b_to[l]
        tr = br.tap_real[l]
        ti = br.tap_imag[l]
        itm2 = br.inv_tap[l]^2

        Y_ff_r = (g + g_fr_v) * itm2
        Y_ff_i = (b_ln + b_fr_v) * itm2
        Y_ft_r = -(g * tr - b_ln * ti) * itm2
        Y_ft_i = -(g * ti + b_ln * tr) * itm2
        Y_tf_r = -(g * tr + b_ln * ti) * itm2
        Y_tf_i = (g * ti - b_ln * tr) * itm2
        Y_tt_r = g + g_to_v
        Y_tt_i = b_ln + b_to_v

        f = br.from[l]
        t = br.to[l]

        for tt in 1:T
            # cr_fr = Y_ff_r·vr_fr - Y_ff_i·vi_fr + Y_ft_r·vr_to - Y_ft_i·vi_to
            @constraint(model,
                cr_fr[l, tt] ==
                Y_ff_r * vr[f, tt] - Y_ff_i * vi[f, tt] +
                Y_ft_r * vr[t, tt] - Y_ft_i * vi[t, tt]
            )
            # ci_fr = Y_ff_i·vr_fr + Y_ff_r·vi_fr + Y_ft_i·vr_to + Y_ft_r·vi_to
            @constraint(model,
                ci_fr[l, tt] ==
                Y_ff_i * vr[f, tt] + Y_ff_r * vi[f, tt] +
                Y_ft_i * vr[t, tt] + Y_ft_r * vi[t, tt]
            )
            # cr_to = Y_tf_r·vr_fr - Y_tf_i·vi_fr + Y_tt_r·vr_to - Y_tt_i·vi_to
            @constraint(model,
                cr_to[l, tt] ==
                Y_tf_r * vr[f, tt] - Y_tf_i * vi[f, tt] +
                Y_tt_r * vr[t, tt] - Y_tt_i * vi[t, tt]
            )
            # ci_to = Y_tf_i·vr_fr + Y_tf_r·vi_fr + Y_tt_i·vr_to + Y_tt_r·vi_to
            @constraint(model,
                ci_to[l, tt] ==
                Y_tf_i * vr[f, tt] + Y_tf_r * vi[f, tt] +
                Y_tt_i * vr[t, tt] + Y_tt_r * vi[t, tt]
            )
        end
    end

    # Thermal limits: |V|²·|I|² ≤ rate² (quartic)
    @inbounds for l in active
        ra2 = br.rate_a[l]^2
        f = br.from[l]
        t = br.to[l]
        for tt in 1:T
            @constraint(model,
                (vr[f, tt]^2 + vi[f, tt]^2) *
                (cr_fr[l, tt]^2 + ci_fr[l, tt]^2) <= ra2
            )
            @constraint(model,
                (vr[t, tt]^2 + vi[t, tt]^2) *
                (cr_to[l, tt]^2 + ci_to[l, tt]^2) <= ra2
            )
        end
    end

    # Angle limits via voltage product (same as ACR)
    @constraint(model,
        branch_angle_lb[l in active, t in 1:T; !isnan(br.ang_min[l])],
        tan(br.ang_min[l]) *
        (vr[br.from[l], t] * vr[br.to[l], t] + vi[br.from[l], t] * vi[br.to[l], t]) <=
        (vi[br.from[l], t] * vr[br.to[l], t] - vr[br.from[l], t] * vi[br.to[l], t])
    )

    @constraint(model,
        branch_angle_ub[l in active, t in 1:T; !isnan(br.ang_max[l])],
        (vi[br.from[l], t] * vr[br.to[l], t] - vr[br.from[l], t] * vi[br.to[l], t]) <=
        tan(br.ang_max[l]) *
        (vr[br.from[l], t] * vr[br.to[l], t] + vi[br.from[l], t] * vi[br.to[l], t])
    )

    return nothing
end
