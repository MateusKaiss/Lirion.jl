# ---------------------------------------------------------------------
# DC Branch Modeling
# ---------------------------------------------------------------------
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

    @constraint(model, [i in params.active, t in 1:T],
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

# ---------------------------------------------------------------------
# Creating SOCP Branch Variables, Expressions and Constraints
# ---------------------------------------------------------------------
function compute_socp_branch_expressions(model, params)
    T = get(model.ext, :time_periods, 1)
    active = findall(i -> params.existing[i] == 1, 1:params.number_of_branches)
    inv_tap = [1.0 / (isnan(t) || t == 0.0 ? 1.0 : t) for t in params.tap_position]

    @variable(model, R[l in active, t = 1:T])
    @variable(model, T_var[l in active, t = 1:T])

    @inbounds for l in active, t in 1:T
        i_from = params.bus_idx_map[params.bus_from[l]]
        i_to = params.bus_idx_map[params.bus_to[l]]

        Vm_from = params.initial_voltage_magnitude[i_from]
        Vm_to = params.initial_voltage_magnitude[i_to]
        Va_from = params.initial_voltage_angle[i_from]
        Va_to = params.initial_voltage_angle[i_to]

        if !any(isnan.([Vm_from, Vm_to, Va_from, Va_to]))
            shift = params.tap_angle
            set_start_value(R[l, t], Vm_from * Vm_to * cos(Va_from - Va_to - shift[l]))
            set_start_value(T_var[l, t], Vm_from * Vm_to * sin(Va_from - Va_to - shift[l]))
        end
    end

    tap = params.tap_position
    shift = params.tap_angle
    tr = tap .* cos.(shift)
    ti = tap .* sin.(shift)

    @expression(model, p_flow_branches[l in active, t in 1:T],
        (
            params.g[l] * inv_tap[l]^2 * model[:v_magnitude_sq][params.bus_from[l], t] +
            (-params.g[l] * tr[l] + params.b[l] * ti[l]) * inv_tap[l]^2 * R[l, t] +
            (-params.b[l] * tr[l] - params.g[l] * ti[l]) * inv_tap[l]^2 * T_var[l, t]
        )
    )

    @expression(model, p_reverse_flow_branches[l in active, t in 1:T],
        (
            params.g[l] * model[:v_magnitude_sq][params.bus_to[l], t] +
            (-params.g[l] * tr[l] - params.b[l] * ti[l]) * inv_tap[l]^2 * R[l, t] +
            (-params.b[l] * tr[l] + params.g[l] * ti[l]) * inv_tap[l]^2 * -T_var[l, t]
        )
    )

    @expression(model, q_flow_branches[l in active, t in 1:T],
        (
            -(params.b[l] + params.b_shunt[l] / 2) * inv_tap[l]^2 *
            model[:v_magnitude_sq][params.bus_from[l], t] -
            (-params.b[l] * tr[l] - params.g[l] * ti[l]) * inv_tap[l]^2 * R[l, t] +
            (-params.g[l] * tr[l] + params.b[l] * ti[l]) * inv_tap[l]^2 * T_var[l, t]
        )
    )

    @expression(model, q_reverse_flow_branches[l in active, t in 1:T],
        (
            -(params.b[l] + params.b_shunt[l] / 2) *
            model[:v_magnitude_sq][params.bus_to[l], t] -
            (-params.b[l] * tr[l] + params.g[l] * ti[l]) * inv_tap[l]^2 * R[l, t] +
            (-params.g[l] * tr[l] - params.b[l] * ti[l]) * inv_tap[l]^2 * -T_var[l, t]
        )
    )
end

function create_socp_branch_constraints(model, params)
    T = get(model.ext, :time_periods, 1)
    active = findall(i -> params.existing[i] == 1, 1:params.number_of_branches)

    @constraint(model, branch_angle_lb[i in active, t in 1:T; !isnan(params.min_ang[i])],
        params.min_ang[i] <=
        model[:theta][params.bus_from[i], t] - model[:theta][params.bus_to[i], t]
    )

    @constraint(model, branch_angle_ub[i in active, t in 1:T; !isnan(params.max_ang[i])],
        model[:theta][params.bus_from[i], t] - model[:theta][params.bus_to[i], t] <=
        params.max_ang[i]
    )

    @constraint(model, [l in active, t in 1:T],
        [
            model[:v_magnitude_sq][params.bus_from[l], t] / sqrt(2),
            model[:v_magnitude_sq][params.bus_to[l], t] / sqrt(2),
            model[:R][l, t],
            model[:T_var][l, t],
        ] in JuMP.RotatedSecondOrderCone()
    )

    @constraint(model, [l in active, t in 1:T],
        [
            params.capacity[l],
            model[:p_flow_branches][l, t],
            model[:q_flow_branches][l, t],
        ] in JuMP.SecondOrderCone()
    )

    @constraint(model, [l in active, t in 1:T],
        [
            params.capacity[l],
            model[:p_reverse_flow_branches][l, t],
            model[:q_reverse_flow_branches][l, t],
        ] in JuMP.SecondOrderCone()
    )

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
