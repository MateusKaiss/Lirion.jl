function get_branch_parameters(branch_data, matpower_data)
    return (
        number_of_branches = length(branch_data["f_bus"]),
        s_base = matpower_data["baseMVA"],
        resistance = branch_data["br_r"],
        reactance = branch_data["br_x"],
        bus_from = branch_data["f_bus"],
        bus_to = branch_data["t_bus"],
        capacity = branch_data["rate_a"],
        min_ang = branch_data["angmin"],
        max_ang = branch_data["angmax"],
        existing = branch_data["br_status"],
        b_shunt = branch_data["br_b"],
        tap_position = branch_data["tap"],
        initial_voltage_magnitude = matpower_data["bus"]["vm"],
        initial_voltage_angle = matpower_data["bus"]["va"],
    )
end

# ---------------------------------------------------------------------
# Creating DC Branch Variables, Expressions and Constraints
# ---------------------------------------------------------------------
function compute_dc_branch_expressions(model, params)
    br_b = -params.reactance ./ (params.resistance .^ 2 .+ params.reactance .^ 2)

    @expression(
        model,
        p_flow_branches[i = 1:params.number_of_branches],
        params.existing[i] == 1 ?
        params.s_base * br_b[i] *
        (model[:theta][params.bus_from[i]] - model[:theta][params.bus_to[i]]) : 0
    )

    @expression(
        model,
        p_reverse_flow_branches[i = 1:params.number_of_branches],
        -p_flow_branches[i]
    )

    @expression(
        model,
        losses_branches[i = 1:params.number_of_branches],
        params.existing[i] == 1 ?
        (p_flow_branches[i]^2 * (params.resistance[i] / params.s_base)) : 0
    )
end

function create_dc_branch_constraints(model, params)
    active_branches = findall(i -> params.existing[i] == 1, 1:params.number_of_branches)

    @constraint(
        model,
        [i in active_branches],
        -params.capacity[i] <= model[:p_flow_branches][i] <= params.capacity[i]
    )
    @constraint(
        model,
        branch_angle_lb[i in active_branches],
        params.min_ang[i] <=
        model[:theta][params.bus_from[i]] - model[:theta][params.bus_to[i]]
    )
    @constraint(
        model,
        branch_angle_ub[i in active_branches],
        model[:theta][params.bus_from[i]] - model[:theta][params.bus_to[i]] <=
        params.max_ang[i]
    )
end

# ---------------------------------------------------------------------
# Creating AC Branch Variables, Expressions and Constraints
# ---------------------------------------------------------------------
function compute_ac_branch_expressions(model, params)
    g = params.resistance ./ (params.resistance .^ 2 .+ params.reactance .^ 2)
    b = -params.reactance ./ (params.resistance .^ 2 .+ params.reactance .^ 2)

    tap_position = copy(params.tap_position)
    for i in eachindex(tap_position)
        if isnan(tap_position[i]) || tap_position[i] == 0.0
            tap_position[i] = 1.0
        end
    end
    inv_tap = 1.0 ./ tap_position

    @expression(
        model,
        p_flow_branches[l = 1:params.number_of_branches],
        params.existing[l] == 1 ?
        params.s_base *
        (
            g[l] * (inv_tap[l] * model[:v_magnitude][params.bus_from[l]])^2 -
            inv_tap[l] * model[:v_magnitude][params.bus_from[l]] *
            model[:v_magnitude][params.bus_to[l]] *
            (
                g[l] *
                cos(model[:theta][params.bus_from[l]] - model[:theta][params.bus_to[l]]) +
                b[l] *
                sin(model[:theta][params.bus_from[l]] - model[:theta][params.bus_to[l]])
            )
        ) : 0
    )

    @expression(
        model,
        p_reverse_flow_branches[l = 1:params.number_of_branches],
        params.existing[l] == 1 ?
        params.s_base *
        (
            g[l] * model[:v_magnitude][params.bus_to[l]]^2 -
            inv_tap[l] * model[:v_magnitude][params.bus_to[l]] *
            model[:v_magnitude][params.bus_from[l]] *
            (
                g[l] *
                cos(model[:theta][params.bus_from[l]] - model[:theta][params.bus_to[l]]) -
                b[l] *
                sin(model[:theta][params.bus_from[l]] - model[:theta][params.bus_to[l]])
            )
        ) : 0
    )

    @expression(
        model,
        q_flow_branches[l = 1:params.number_of_branches],
        params.existing[l] == 1 ?
        params.s_base *
        (
            -(b[l] + params.b_shunt[l] / 2) *
            (inv_tap[l] * model[:v_magnitude][params.bus_from[l]])^2 -
            inv_tap[l] * model[:v_magnitude][params.bus_from[l]] *
            model[:v_magnitude][params.bus_to[l]] *
            (
                -b[l] *
                cos(model[:theta][params.bus_from[l]] - model[:theta][params.bus_to[l]]) +
                g[l] *
                sin(model[:theta][params.bus_from[l]] - model[:theta][params.bus_to[l]])
            )
        ) : 0
    )

    @expression(
        model,
        q_reverse_flow_branches[l = 1:params.number_of_branches],
        params.existing[l] == 1 ?
        params.s_base *
        (
            -(b[l] + params.b_shunt[l] / 2) * model[:v_magnitude][params.bus_to[l]]^2 -
            inv_tap[l] * model[:v_magnitude][params.bus_to[l]] *
            model[:v_magnitude][params.bus_from[l]] *
            (
                -b[l] *
                cos(model[:theta][params.bus_from[l]] - model[:theta][params.bus_to[l]]) -
                g[l] *
                sin(model[:theta][params.bus_from[l]] - model[:theta][params.bus_to[l]])
            )
        ) : 0
    )
end

function create_ac_branch_constraints(model, params)
    active_branches = findall(i -> params.existing[i] == 1, 1:params.number_of_branches)

    @variable(model, s_branches[1:params.number_of_branches] >= 0)
    @variable(model, s_reverse_branches[1:params.number_of_branches] >= 0)

    @constraint(
        model,
        [l in active_branches],
        s_branches[l]^2 >= model[:p_flow_branches][l]^2 + model[:q_flow_branches][l]^2
    )

    @constraint(
        model,
        [l in active_branches],
        s_branches[l] <= params.capacity[l]
    )

    @constraint(
        model,
        [l in active_branches],
        s_reverse_branches[l]^2 >=
        model[:p_reverse_flow_branches][l]^2 + model[:q_reverse_flow_branches][l]^2
    )

    @constraint(
        model,
        [l in active_branches],
        s_reverse_branches[l] <= params.capacity[l]
    )

    @constraint(
        model,
        branch_angle_lb[i in active_branches; !isnan(params.min_ang[i])],
        params.min_ang[i] <=
        model[:theta][params.bus_from[i]] - model[:theta][params.bus_to[i]]
    )

    @constraint(
        model,
        branch_angle_ub[i in active_branches; !isnan(params.max_ang[i])],
        model[:theta][params.bus_from[i]] - model[:theta][params.bus_to[i]] <=
        params.max_ang[i]
    )
end

# ---------------------------------------------------------------------
# Creating SOCP Branch Variables, Expressions and Constraints
# ---------------------------------------------------------------------
function compute_socp_branch_expressions(model, params)
    g = params.resistance ./ (params.resistance .^ 2 .+ params.reactance .^ 2)
    b = -params.reactance ./ (params.resistance .^ 2 .+ params.reactance .^ 2)

    tap_position = copy(params.tap_position)
    for i in eachindex(tap_position)
        if isnan(tap_position[i]) || tap_position[i] == 0.0
            tap_position[i] = 1.0
        end
    end
    inv_tap = 1.0 ./ tap_position

    active_branches = findall(i -> params.existing[i] == 1, 1:params.number_of_branches)

    @variable(model, R[1:params.number_of_branches])
    @variable(model, T[1:params.number_of_branches])

    for i in active_branches
        V_from = params.initial_voltage_magnitude[params.bus_from[i]]
        V_to = params.initial_voltage_magnitude[params.bus_to[i]]
        theta_from = params.initial_voltage_angle[params.bus_from[i]]
        theta_to = params.initial_voltage_angle[params.bus_to[i]]

        if isnan(V_from) || isnan(V_to) || isnan(theta_from) || isnan(theta_to)
            continue
        end

        set_start_value(R[i], V_from * V_to * cos(theta_from - theta_to))
        set_start_value(T[i], V_from * V_to * sin(theta_from - theta_to))
    end

    @expression(
        model,
        p_flow_branches[l = 1:params.number_of_branches],
        params.existing[l] == 1 ?
        params.s_base * (
            g[l] * (inv_tap[l]^2 * model[:v_magnitude_sq][params.bus_from[l]]) -
            inv_tap[l] * g[l] * R[l] - inv_tap[l] * b[l] * T[l]
        ) : 0
    )

    @expression(
        model,
        p_reverse_flow_branches[l = 1:params.number_of_branches],
        params.existing[l] == 1 ?
        params.s_base * (
            g[l] * model[:v_magnitude_sq][params.bus_to[l]] -
            inv_tap[l] * g[l] * R[l] + inv_tap[l] * b[l] * T[l]
        ) : 0
    )

    @expression(
        model,
        q_flow_branches[l = 1:params.number_of_branches],
        params.existing[l] == 1 ?
        params.s_base * (
            -(b[l] + (params.b_shunt[l] / 2)) *
            (inv_tap[l]^2 * model[:v_magnitude_sq][params.bus_from[l]]) +
            inv_tap[l] * b[l] * R[l] - inv_tap[l] * g[l] * T[l]) : 0
    )

    @expression(
        model,
        q_reverse_flow_branches[l = 1:params.number_of_branches],
        params.existing[l] == 1 ?
        params.s_base * (
            -(b[l] + (params.b_shunt[l] / 2)) * model[:v_magnitude_sq][params.bus_to[l]] +
            inv_tap[l] * b[l] * R[l] + inv_tap[l] * g[l] * T[l]
        ) : 0
    )
end

function create_socp_branch_constraints(model, params)
    active_branches = findall(i -> params.existing[i] == 1, 1:params.number_of_branches)

    @constraint(
        model,
        branch_angle_lb[i in active_branches; !isnan(params.min_ang[i])],
        params.min_ang[i] <=
        model[:theta][params.bus_from[i]] - model[:theta][params.bus_to[i]]
    )

    @constraint(
        model,
        branch_angle_ub[i in active_branches; !isnan(params.max_ang[i])],
        model[:theta][params.bus_from[i]] - model[:theta][params.bus_to[i]] <=
        params.max_ang[i]
    )

    @constraint(
        model,
        [m in active_branches],
        [
            model[:v_magnitude_sq][params.bus_from[m]] / sqrt(2),
            model[:v_magnitude_sq][params.bus_to[m]] / sqrt(2),
            model[:R][m],
            model[:T][m],
        ] in JuMP.RotatedSecondOrderCone()
    )

    @constraint(
        model,
        [l in active_branches],
        [params.capacity[l], model[:p_flow_branches][l], model[:q_flow_branches][l]] in
        JuMP.SecondOrderCone()
    )

    @constraint(
        model,
        [l in active_branches],
        [
            params.capacity[l],
            model[:p_reverse_flow_branches][l],
            model[:q_reverse_flow_branches][l],
        ] in
        JuMP.SecondOrderCone()
    )
end

# ------------------------------------------------------------------------
# Dispatching Branch Variables, Expressions and Constraints by Model Type
# ------------------------------------------------------------------------
function branches!(model::JuMP.Model, matpower_data, ::Type{DC}, ::Type{SubproblemBuild})
    branch_data = matpower_data["branch"]
    params = get_branch_parameters(branch_data, matpower_data)
    compute_dc_branch_expressions(model, params)
    create_dc_branch_constraints(model, params)
    return nothing
end

function branches!(model::JuMP.Model, matpower_data, ::Type{AC}, ::Type{SubproblemBuild})
    branch_data = matpower_data["branch"]
    params = get_branch_parameters(branch_data, matpower_data)
    compute_ac_branch_expressions(model, params)
    create_ac_branch_constraints(model, params)
    return nothing
end

function branches!(model::JuMP.Model, matpower_data, ::Type{SOCP}, ::Type{SubproblemBuild})
    branch_data = matpower_data["branch"]
    params = get_branch_parameters(branch_data, matpower_data)
    compute_socp_branch_expressions(model, params)
    create_socp_branch_constraints(model, params)
    return nothing
end
