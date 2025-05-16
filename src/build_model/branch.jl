# ---------------------------------------------------------------------
# Getting Branch Parameters (multi-temporal compatible, no mudança)
# ---------------------------------------------------------------------
function get_branch_parameters(branch_data, matpower_data)
    bus_ids = matpower_data["bus"]["bus_i"]
    bus_idx_map = Dict(b => i for (i, b) in enumerate(bus_ids))

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
        bus_ids = bus_ids,
        bus_idx_map = bus_idx_map,
    )
end

# ---------------------------------------------------------------------
# DC Branch Expressions (multi-temporal)
# ---------------------------------------------------------------------
function compute_dc_branch_expressions(model, params)
    T = get(model.ext, :time_periods, 1)  # default = 1 (single-period)
    br_b = -params.reactance ./ (params.resistance .^ 2 .+ params.reactance .^ 2)

    @expression(
        model,
        p_flow_branches[i = 1:params.number_of_branches, t = 1:T],
        params.existing[i] == 1 ?
        params.s_base * br_b[i] *
        (model[:theta][params.bus_from[i], t] - model[:theta][params.bus_to[i], t]) : 0
    )

    @expression(
        model,
        p_reverse_flow_branches[i = 1:params.number_of_branches, t = 1:T],
        -model[:p_flow_branches][i, t]
    )

    @expression(
        model,
        losses_branches[i = 1:params.number_of_branches, t = 1:T],
        params.existing[i] == 1 ?
        (model[:p_flow_branches][i, t]^2 * (params.resistance[i] / params.s_base)) : 0
    )
end

# ---------------------------------------------------------------------
# DC Branch Constraints (multi-temporal)
# ---------------------------------------------------------------------
function create_dc_branch_constraints(model, params)
    T = get(model.ext, :time_periods, 1)
    active_branches = findall(i -> params.existing[i] == 1, 1:params.number_of_branches)

    @constraint(
        model,
        [i in active_branches, t = 1:T],
        -params.capacity[i] <= model[:p_flow_branches][i, t] <= params.capacity[i]
    )

    @constraint(
        model,
        branch_angle_lb[i in active_branches, t = 1:T],
        params.min_ang[i] <=
        model[:theta][params.bus_from[i], t] - model[:theta][params.bus_to[i], t]
    )

    @constraint(
        model,
        branch_angle_ub[i in active_branches, t = 1:T],
        model[:theta][params.bus_from[i], t] - model[:theta][params.bus_to[i], t] <=
        params.max_ang[i]
    )
end

# ---------------------------------------------------------------------
# Creating AC Branch Variables, Expressions and Constraints
# ---------------------------------------------------------------------
function compute_ac_branch_expressions(model, params)
    T = get(model.ext, :time_periods, 1)
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
        p_flow_branches[l = 1:params.number_of_branches, t = 1:T],
        params.existing[l] == 1 ?
        params.s_base *
        (
            g[l] * (inv_tap[l] * model[:v_magnitude][params.bus_from[l], t])^2 -
            inv_tap[l] * model[:v_magnitude][params.bus_from[l], t] *
            model[:v_magnitude][params.bus_to[l], t] *
            (
                g[l] *
                cos(
                    model[:theta][params.bus_from[l], t] -
                    model[:theta][params.bus_to[l], t],
                ) +
                b[l] *
                sin(
                    model[:theta][params.bus_from[l], t] -
                    model[:theta][params.bus_to[l], t],
                )
            )
        ) : 0
    )

    @expression(
        model,
        p_reverse_flow_branches[l = 1:params.number_of_branches, t = 1:T],
        params.existing[l] == 1 ?
        params.s_base *
        (
            g[l] * model[:v_magnitude][params.bus_to[l], t]^2 -
            inv_tap[l] * model[:v_magnitude][params.bus_to[l], t] *
            model[:v_magnitude][params.bus_from[l], t] *
            (
                g[l] *
                cos(
                    model[:theta][params.bus_from[l], t] -
                    model[:theta][params.bus_to[l], t],
                ) -
                b[l] *
                sin(
                    model[:theta][params.bus_from[l], t] -
                    model[:theta][params.bus_to[l], t],
                )
            )
        ) : 0
    )

    @expression(
        model,
        q_flow_branches[l = 1:params.number_of_branches, t = 1:T],
        params.existing[l] == 1 ?
        params.s_base *
        (
            -(b[l] + params.b_shunt[l] / 2) *
            (inv_tap[l] * model[:v_magnitude][params.bus_from[l], t])^2 -
            inv_tap[l] * model[:v_magnitude][params.bus_from[l], t] *
            model[:v_magnitude][params.bus_to[l], t] *
            (
                -b[l] *
                cos(
                    model[:theta][params.bus_from[l], t] -
                    model[:theta][params.bus_to[l], t],
                ) +
                g[l] *
                sin(
                    model[:theta][params.bus_from[l], t] -
                    model[:theta][params.bus_to[l], t],
                )
            )
        ) : 0
    )

    @expression(
        model,
        q_reverse_flow_branches[l = 1:params.number_of_branches, t = 1:T],
        params.existing[l] == 1 ?
        params.s_base *
        (
            -(b[l] + params.b_shunt[l] / 2) * model[:v_magnitude][params.bus_to[l], t]^2 -
            inv_tap[l] * model[:v_magnitude][params.bus_to[l], t] *
            model[:v_magnitude][params.bus_from[l], t] *
            (
                -b[l] *
                cos(
                    model[:theta][params.bus_from[l], t] -
                    model[:theta][params.bus_to[l], t],
                ) -
                g[l] *
                sin(
                    model[:theta][params.bus_from[l], t] -
                    model[:theta][params.bus_to[l], t],
                )
            )
        ) : 0
    )
end

function create_ac_branch_constraints(model, params)
    T = get(model.ext, :time_periods, 1)
    active_branches = findall(i -> params.existing[i] == 1, 1:params.number_of_branches)

    @variable(model, s_branches[1:params.number_of_branches, t = 1:T] >= 0)
    @variable(model, s_reverse_branches[1:params.number_of_branches, t = 1:T] >= 0)

    @constraint(
        model,
        [l in active_branches, t = 1:T],
        s_branches[l, t]^2 >=
        model[:p_flow_branches][l, t]^2 + model[:q_flow_branches][l, t]^2
    )

    @constraint(
        model,
        [l in active_branches, t = 1:T],
        s_branches[l, t] <= params.capacity[l]
    )

    @constraint(
        model,
        [l in active_branches, t = 1:T],
        s_reverse_branches[l, t]^2 >=
        model[:p_reverse_flow_branches][l, t]^2 + model[:q_reverse_flow_branches][l, t]^2
    )

    @constraint(
        model,
        [l in active_branches, t = 1:T],
        s_reverse_branches[l, t] <= params.capacity[l]
    )

    @constraint(
        model,
        branch_angle_lb[i in active_branches, t = 1:T; !isnan(params.min_ang[i])],
        params.min_ang[i] <=
        model[:theta][params.bus_from[i], t] - model[:theta][params.bus_to[i], t]
    )

    @constraint(
        model,
        branch_angle_ub[i in active_branches, t = 1:T; !isnan(params.max_ang[i])],
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

    @variable(model, R[l in active_branches, t = 1:T])
    @variable(model, T_var[l in active_branches, t = 1:T])

    for t in 1:T, l in active_branches
        b_from = params.bus_from[l]
        b_to = params.bus_to[l]

        idx_from = params.bus_idx_map[b_from]
        idx_to = params.bus_idx_map[b_to]

        V_from = params.initial_voltage_magnitude[idx_from]
        V_to = params.initial_voltage_magnitude[idx_to]
        theta_from = params.initial_voltage_angle[idx_from]
        theta_to = params.initial_voltage_angle[idx_to]

        if !any(isnan.([V_from, V_to, theta_from, theta_to]))
            set_start_value(R[l, t], V_from * V_to * cos(theta_from - theta_to))
            set_start_value(T_var[l, t], V_from * V_to * sin(theta_from - theta_to))
        end
    end

    @expression(
        model,
        p_flow_branches[l = 1:params.number_of_branches, t = 1:T],
        params.existing[l] == 1 ?
        params.s_base * (
            g[l] * (inv_tap[l]^2 * model[:v_magnitude_sq][params.bus_from[l], t]) -
            inv_tap[l] * g[l] * R[l, t] - inv_tap[l] * b[l] * T_var[l, t]
        ) : 0
    )

    @expression(
        model,
        p_reverse_flow_branches[l = 1:params.number_of_branches, t = 1:T],
        params.existing[l] == 1 ?
        params.s_base * (
            g[l] * model[:v_magnitude_sq][params.bus_to[l], t] -
            inv_tap[l] * g[l] * R[l, t] + inv_tap[l] * b[l] * T_var[l, t]
        ) : 0
    )

    @expression(
        model,
        q_flow_branches[l = 1:params.number_of_branches, t = 1:T],
        params.existing[l] == 1 ?
        params.s_base * (
            -(b[l] + (params.b_shunt[l] / 2)) *
            (inv_tap[l]^2 * model[:v_magnitude_sq][params.bus_from[l], t]) +
            inv_tap[l] * b[l] * R[l, t] - inv_tap[l] * g[l] * T_var[l, t]
        ) : 0
    )

    @expression(
        model,
        q_reverse_flow_branches[l = 1:params.number_of_branches, t = 1:T],
        params.existing[l] == 1 ?
        params.s_base * (
            -(b[l] + (params.b_shunt[l] / 2)) *
            model[:v_magnitude_sq][params.bus_to[l], t] +
            inv_tap[l] * b[l] * R[l, t] + inv_tap[l] * g[l] * T_var[l, t]
        ) : 0
    )
end

function create_socp_branch_constraints(model, params)
    T = get(model.ext, :time_periods, 1)
    active_branches = findall(i -> params.existing[i] == 1, 1:params.number_of_branches)

    @constraint(
        model,
        branch_angle_lb[i in active_branches, t = 1:T; !isnan(params.min_ang[i])],
        params.min_ang[i] <=
        model[:theta][params.bus_from[i], t] - model[:theta][params.bus_to[i], t]
    )

    @constraint(
        model,
        branch_angle_ub[i in active_branches, t = 1:T; !isnan(params.max_ang[i])],
        model[:theta][params.bus_from[i], t] - model[:theta][params.bus_to[i], t] <=
        params.max_ang[i]
    )

    @constraint(
        model,
        [m in active_branches, t = 1:T],
        [
            model[:v_magnitude_sq][params.bus_from[m], t] / sqrt(2),
            model[:v_magnitude_sq][params.bus_to[m], t] / sqrt(2),
            model[:R][m, t],
            model[:T_var][m, t],
        ] in JuMP.RotatedSecondOrderCone()
    )

    # @constraint(
    #     model,
    #     [l in active_branches, t = 1:T],
    #     [
    #         (model[:v_magnitude_sq][params.bus_from[l], t] + model[:v_magnitude_sq][params.bus_to[l], t]) / 2,
    #         model[:R][l, t],
    #         model[:T_var][l, t]
    #     ] in JuMP.SecondOrderCone()
    # )

    @constraint(
        model,
        [l in active_branches, t = 1:T],
        [
            params.capacity[l],
            model[:p_flow_branches][l, t],
            model[:q_flow_branches][l, t],
        ] in
        JuMP.SecondOrderCone()
    )

    @constraint(
        model,
        [l in active_branches, t = 1:T],
        [
            params.capacity[l],
            model[:p_reverse_flow_branches][l, t],
            model[:q_reverse_flow_branches][l, t],
        ] in
        JuMP.SecondOrderCone()
    )

    return nothing
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
