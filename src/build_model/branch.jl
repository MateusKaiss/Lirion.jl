function compute_dc_branch_expressions(model, params)
    T = get(model.ext, :time_periods, 1)

    @expression(model, p_flow_branches[i in params.active, t in 1:T],
        params.b[i] *
        (model[:theta][params.bus_from[i], t] - model[:theta][params.bus_to[i], t])
    )

    @expression(model, p_reverse_flow_branches[i in params.active, t in 1:T],
        -model[:p_flow_branches][i, t]
    )

    @expression(model, losses_branches[i in params.active, t in 1:T],
        model[:p_flow_branches][i, t]^2 * (params.resistance[i])
    )
end

function create_dc_branch_constraints(model, params)
    T = get(model.ext, :time_periods, 1)

    @constraint(model, capacity[i in params.active, t in 1:T],
        -params.capacity[i] <= model[:p_flow_branches][i, t] <= params.capacity[i]
    )

    @constraint(model,
        branch_angle_lb[i in params.active, t in 1:T; !isnan(params.min_ang[i])],
        params.min_ang[i] <=
        model[:theta][params.bus_from[i], t] - model[:theta][params.bus_to[i], t]
    )

    @constraint(model,
        branch_angle_ub[i in params.active, t in 1:T; !isnan(params.max_ang[i])],
        model[:theta][params.bus_from[i], t] - model[:theta][params.bus_to[i], t] <=
        params.max_ang[i]
    )
end

# ---------------------------------------------------------------------
# AC Branch Modeling
# ---------------------------------------------------------------------
function compute_ac_branch_expressions(model, params)
    T = get(model.ext, :time_periods, 1)

    @expression(model, p_flow_branches[l in params.active, t in 1:T],
        (
            params.g[l] *
            (params.inv_tap[l] * model[:v_magnitude][params.bus_from[l], t])^2 -
            params.inv_tap[l] * model[:v_magnitude][params.bus_from[l], t] *
            model[:v_magnitude][params.bus_to[l], t] *
            (
                params.g[l] * cos(
                    model[:theta][params.bus_from[l], t] -
                    model[:theta][params.bus_to[l], t] - params.tap_angle[l],
                ) +
                params.b[l] * sin(
                    model[:theta][params.bus_from[l], t] -
                    model[:theta][params.bus_to[l], t] - params.tap_angle[l],
                )
            )
        )
    )

    @expression(model, p_reverse_flow_branches[l in params.active, t in 1:T],
        (
            params.g[l] * model[:v_magnitude][params.bus_to[l], t]^2 -
            params.inv_tap[l] * model[:v_magnitude][params.bus_to[l], t] *
            model[:v_magnitude][params.bus_from[l], t] *
            (
                params.g[l] * cos(
                    model[:theta][params.bus_from[l], t] -
                    model[:theta][params.bus_to[l], t] - params.tap_angle[l],
                ) -
                params.b[l] * sin(
                    model[:theta][params.bus_from[l], t] -
                    model[:theta][params.bus_to[l], t] - params.tap_angle[l],
                )
            )
        )
    )

    @expression(model, q_flow_branches[l in params.active, t in 1:T],
        (
            -(params.b[l] + params.b_shunt[l] / 2) *
            (params.inv_tap[l] * model[:v_magnitude][params.bus_from[l], t])^2 -
            params.inv_tap[l] * model[:v_magnitude][params.bus_from[l], t] *
            model[:v_magnitude][params.bus_to[l], t] *
            (
                -params.b[l] * cos(
                    model[:theta][params.bus_from[l], t] -
                    model[:theta][params.bus_to[l], t] - params.tap_angle[l],
                ) +
                params.g[l] * sin(
                    model[:theta][params.bus_from[l], t] -
                    model[:theta][params.bus_to[l], t] - params.tap_angle[l],
                )
            )
        )
    )

    @expression(model, q_reverse_flow_branches[l in params.active, t in 1:T],
        (
            -(params.b[l] + params.b_shunt[l] / 2) *
            model[:v_magnitude][params.bus_to[l], t]^2 -
            params.inv_tap[l] * model[:v_magnitude][params.bus_to[l], t] *
            model[:v_magnitude][params.bus_from[l], t] *
            (
                -params.b[l] * cos(
                    model[:theta][params.bus_from[l], t] -
                    model[:theta][params.bus_to[l], t] - params.tap_angle[l],
                ) -
                params.g[l] * sin(
                    model[:theta][params.bus_from[l], t] -
                    model[:theta][params.bus_to[l], t] - params.tap_angle[l],
                )
            )
        )
    )
end

function create_ac_branch_constraints(model, params)
    T = get(model.ext, :time_periods, 1)

    @variable(model, s_branches[l in params.active, t = 1:T] >= 0)
    @variable(model, s_reverse_branches[l in params.active, t = 1:T] >= 0)

    @constraint(model, [l in params.active, t in 1:T],
        s_branches[l, t]^2 >=
        model[:p_flow_branches][l, t]^2 + model[:q_flow_branches][l, t]^2
    )

    @constraint(model, [l in params.active, t in 1:T],
        s_branches[l, t] <= params.capacity[l]
    )

    @constraint(model, [l in params.active, t in 1:T],
        s_reverse_branches[l, t]^2 >=
        model[:p_reverse_flow_branches][l, t]^2 + model[:q_reverse_flow_branches][l, t]^2
    )

    @constraint(model, [l in params.active, t in 1:T],
        s_reverse_branches[l, t] <= params.capacity[l]
    )

    @constraint(model,
        branch_angle_lb[i in params.active, t in 1:T; !isnan(params.min_ang[i])],
        params.min_ang[i] <=
        model[:theta][params.bus_from[i], t] - model[:theta][params.bus_to[i], t]
    )

    @constraint(model,
        branch_angle_ub[i in params.active, t in 1:T; !isnan(params.max_ang[i])],
        model[:theta][params.bus_from[i], t] - model[:theta][params.bus_to[i], t] <=
        params.max_ang[i]
    )
    return nothing
end


buspair(i, j) = (max(i, j), min(i, j))

function compute_socp_branch_expressions(model, params)
    T = get(model.ext, :time_periods, 1)
    active = findall(i -> params.existing[i] == 1, 1:params.number_of_branches)

    w     = model[:v_magnitude_sq]
    tap   = params.tap_position
    shift = params.tap_angle
    tr    = tap .* cos.(shift)
    ti    = tap .* sin.(shift)
    it2   = [1.0 / (isnan(t) || t == 0.0 ? 1.0 : t)^2 for t in tap]  # 1/|tap|^2

    # ---------- wr, wi (por par de barras) ----------
    @variable(model,
        wr[p in Set(buspair(params.bus_from[l], params.bus_to[l]) for l in active), t in 1:T],
        base_name = "0_wr",
    )
    @variable(model,
        wi[p in Set(buspair(params.bus_from[l], params.bus_to[l]) for l in active), t in 1:T],
        base_name = "0_wi",
    )

    # bounds wr/wi como no PM
    for l in active, t in 1:T
        i = params.bus_from[l]; j = params.bus_to[l]
        pkey = buspair(i, j)

        vmin_i = params.vmin[i]; vmax_i = params.vmax[i]
        vmin_j = params.vmin[j]; vmax_j = params.vmax[j]
        admax  = max(abs(params.min_ang[l]), abs(params.max_ang[l]))

        wr_lb = vmin_i * vmin_j * cos(admax)
        wr_ub = vmax_i * vmax_j
        JuMP.set_lower_bound(wr[pkey, t], wr_lb)
        JuMP.set_upper_bound(wr[pkey, t], wr_ub)

        wi_abs = vmax_i * vmax_j * sin(admax)
        JuMP.set_lower_bound(wi[pkey, t], -wi_abs)
        JuMP.set_upper_bound(wi[pkey, t],  wi_abs)
    end

    # ---------- fluxos p,q em DOIS sentidos ----------
    @variable(model, p_fr[(l,t) in Iterators.product(active, 1:T)], base_name = "0_p_fr")
    @variable(model, q_fr[(l,t) in Iterators.product(active, 1:T)], base_name = "0_q_fr")
    @variable(model, p_to[(l,t) in Iterators.product(active, 1:T)], base_name = "0_p_to")
    @variable(model, q_to[(l,t) in Iterators.product(active, 1:T)], base_name = "0_q_to")

    for l in active, t in 1:T
        cap = params.capacity[l]
        set_lower_bound(p_fr[(l,t)], -cap); set_upper_bound(p_fr[(l,t)],  cap)
        set_lower_bound(q_fr[(l,t)], -cap); set_upper_bound(q_fr[(l,t)],  cap)
        set_lower_bound(p_to[(l,t)], -cap); set_upper_bound(p_to[(l,t)],  cap)
        set_lower_bound(q_to[(l,t)], -cap); set_upper_bound(q_to[(l,t)],  cap)
    end

    model[:wr] = wr
    model[:wi] = wi
    model[:p_flow_branches]         = p_fr    # from → to
    model[:q_flow_branches]         = q_fr
    model[:p_reverse_flow_branches] = p_to    # to → from
    model[:q_reverse_flow_branches] = q_to

    @inbounds for l in active, t in 1:T
        i = params.bus_from[l]; j = params.bus_to[l]
        pkey = buspair(i, j)

        g = params.g[l]
        b = params.b[l]

        b_sh_total = abs(params.b_shunt[l])   
        b_fr = b_sh_total/2
        b_to = b_sh_total/2

        # p (from → to)
        c_w_fr_p  = g * it2[l]
        c_wr_fr_p = ( g*tr[l] -  b*ti[l]) * it2[l]
        c_wi_fr_p = ( b*tr[l] +  g*ti[l]) * it2[l]

        # q (from → to)
        c_w_fr_q  = (b + b_fr) * it2[l]
        c_wr_fr_q = -(b*tr[l] + g*ti[l]) * it2[l]
        c_wi_fr_q =  (g*tr[l] - b*ti[l]) * it2[l]

        # p (to → from)
        c_w_to_p  =  g
        c_wr_to_p = ( g*tr[l] +  b*ti[l]) * it2[l]
        c_wi_to_p = -(b*tr[l] -  g*ti[l]) * it2[l]

        # q (to → from)
        c_w_to_q  = (b + b_to)
        c_wr_to_q = -(b*tr[l] - g*ti[l]) * it2[l]
        c_wi_to_q = -(g*tr[l] + b*ti[l]) * it2[l]

        @constraint(model, -c_w_fr_p*w[i,t] + c_wr_fr_p*wr[pkey,t] + c_wi_fr_p*wi[pkey,t] + p_fr[(l,t)] == 0)
        @constraint(model, c_w_fr_q*w[i,t] + c_wr_fr_q*wr[pkey,t] + c_wi_fr_q*wi[pkey,t] + q_fr[(l,t)] == 0)

        @constraint(model, -c_w_to_p*w[j,t] + c_wr_to_p*wr[pkey,t] + c_wi_to_p*wi[pkey,t] + p_to[(l,t)] == 0)
        @constraint(model, c_w_to_q*w[j,t] + c_wr_to_q*wr[pkey,t] + c_wi_to_q*wi[pkey,t] + q_to[(l,t)] == 0)
    end

    return nothing
end


function create_socp_branch_constraints(model, params)
    T = get(model.ext, :time_periods, 1)
    active = findall(i -> params.existing[i] == 1, 1:params.number_of_branches)

    w  = model[:v_magnitude_sq]
    wr = model[:wr]
    wi = model[:wi]
    p_fr = model[:p_flow_branches]
    q_fr = model[:q_flow_branches]
    p_to = model[:p_reverse_flow_branches]
    q_to = model[:q_reverse_flow_branches]

    @constraint(model, [l in active, t in 1:T],
        [
            w[params.bus_from[l], t] / sqrt(2),
            w[params.bus_to[l],   t] / sqrt(2),
            wr[buspair(params.bus_from[l], params.bus_to[l]), t],
            wi[buspair(params.bus_from[l], params.bus_to[l]), t],
        ] in JuMP.RotatedSecondOrderCone()
    )

    @constraint(model, [l in active, t in 1:T],
        [params.capacity[l], p_fr[(l,t)], q_fr[(l,t)]] in JuMP.SecondOrderCone()
    )
    @constraint(model, [l in active, t in 1:T],
        [params.capacity[l], p_to[(l,t)], q_to[(l,t)]] in JuMP.SecondOrderCone()
    )

    for l in active, t in 1:T
        i = params.bus_from[l]; j = params.bus_to[l]
        pkey = buspair(i, j)

        admax = max(abs(params.min_ang[l]), abs(params.max_ang[l]))
        τ = tan(admax)

        @constraint(model,  τ*wr[pkey,t] + wi[pkey,t] >= 0)
        @constraint(model, -τ*wr[pkey,t] + wi[pkey,t] <= 0)
    end

    return nothing
end

# ------------------------------------------------------------------------
# Dispatching Branch Variables, Expressions and Constraints by Model Type
# ------------------------------------------------------------------------
function branches!(model::JuMP.Model, matpower_data, ::Type{DC}, ::Type{SubproblemBuild})
    params = matpower_data["branch_params"]
    compute_dc_branch_expressions(model, params)
    create_dc_branch_constraints(model, params)
    return nothing
end

function branches!(model::JuMP.Model, matpower_data, ::Type{AC}, ::Type{SubproblemBuild})
    params = matpower_data["branch_params"]
    compute_ac_branch_expressions(model, params)
    create_ac_branch_constraints(model, params)
    return nothing
end

function branches!(model::JuMP.Model, matpower_data, ::Type{SOCP}, ::Type{SubproblemBuild})
    params = matpower_data["branch_params"]
    compute_socp_branch_expressions(model, params)
    create_socp_branch_constraints(model, params)
    return nothing
end