# ---------------------------------------------------------------------
# Getting Branch Parameters (multi-temporal compatible, no mudança)
# ---------------------------------------------------------------------
function get_branch_parameters(branch_data, matpower_data)
    return (
        # Informações topológicas
        number_of_branches = length(branch_data["f_bus"]),
        bus_from = branch_data["f_bus"],
        bus_to = branch_data["t_bus"],
        bus_ids = matpower_data["bus"]["bus_i"],
        bus_idx_map = matpower_data["bus_idx_map"],

        # Impedâncias e tap
        resistance = branch_data["br_r"],
        reactance = branch_data["br_x"],
        b_shunt = branch_data["br_b"],
        tap_position = branch_data["tap"],

        # Limites e status
        capacity = branch_data["rate_a"],
        min_ang = branch_data["angmin"],
        max_ang = branch_data["angmax"],
        existing = branch_data["br_status"],

        # Dados do sistema base e condições iniciais
        s_base = matpower_data["baseMVA"],
        initial_voltage_magnitude = matpower_data["bus"]["vm"],
        initial_voltage_angle = matpower_data["bus"]["va"],
    )
end
function compute_dc_branch_expressions(model, params)
    T = get(model.ext, :time_periods, 1)
    active = findall(i -> params.existing[i] == 1, 1:params.number_of_branches)
    br_b = -params.reactance ./ (params.resistance .^ 2 .+ params.reactance .^ 2)

    @expression(model, p_flow_branches[i in active, t in 1:T],
        params.s_base * br_b[i] *
        (model[:theta][params.bus_from[i], t] - model[:theta][params.bus_to[i], t])
    )

    @expression(model, p_reverse_flow_branches[i in active, t in 1:T],
        -model[:p_flow_branches][i, t]
    )

    @expression(model, losses_branches[i in active, t in 1:T],
        model[:p_flow_branches][i, t]^2 * (params.resistance[i] / params.s_base)
    )
end

function create_dc_branch_constraints(model, params)
    T = get(model.ext, :time_periods, 1)
    active = findall(i -> params.existing[i] == 1, 1:params.number_of_branches)

    @constraint(model, [i in active, t in 1:T],
        -params.capacity[i] <= model[:p_flow_branches][i, t] <= params.capacity[i]
    )

    @constraint(model, branch_angle_lb[i in active, t in 1:T; !isnan(params.min_ang[i])],
        params.min_ang[i] <=
        model[:theta][params.bus_from[i], t] - model[:theta][params.bus_to[i], t]
    )

    @constraint(model, branch_angle_ub[i in active, t in 1:T; !isnan(params.max_ang[i])],
        model[:theta][params.bus_from[i], t] - model[:theta][params.bus_to[i], t] <=
        params.max_ang[i]
    )
end

function compute_ac_branch_expressions(model, params)
    T = get(model.ext, :time_periods, 1)
    active = findall(i -> params.existing[i] == 1, 1:params.number_of_branches)

    g = params.resistance ./ (params.resistance .^ 2 .+ params.reactance .^ 2)
    b = -params.reactance ./ (params.resistance .^ 2 .+ params.reactance .^ 2)

    # Tratar taps inválidos (0 ou NaN)
    tap = [isnan(t) || t == 0.0 ? 1.0 : t for t in params.tap_position]
    inv_tap = 1.0 ./ tap

    @expression(model, p_flow_branches[l in active, t in 1:T],
        params.s_base * (
            g[l] * (inv_tap[l] * model[:v_magnitude][params.bus_from[l], t])^2 -
            inv_tap[l] * model[:v_magnitude][params.bus_from[l], t] *
            model[:v_magnitude][params.bus_to[l], t] *
            (
                g[l] * cos(
                    model[:theta][params.bus_from[l], t] -
                    model[:theta][params.bus_to[l], t],
                ) +
                b[l] * sin(
                    model[:theta][params.bus_from[l], t] -
                    model[:theta][params.bus_to[l], t],
                )
            )
        )
    )

    @expression(model, p_reverse_flow_branches[l in active, t in 1:T],
        params.s_base * (
            g[l] * model[:v_magnitude][params.bus_to[l], t]^2 -
            inv_tap[l] * model[:v_magnitude][params.bus_to[l], t] *
            model[:v_magnitude][params.bus_from[l], t] *
            (
                g[l] * cos(
                    model[:theta][params.bus_from[l], t] -
                    model[:theta][params.bus_to[l], t],
                ) -
                b[l] * sin(
                    model[:theta][params.bus_from[l], t] -
                    model[:theta][params.bus_to[l], t],
                )
            )
        )
    )

    @expression(model, q_flow_branches[l in active, t in 1:T],
        params.s_base * (
            -(b[l] + params.b_shunt[l] / 2) *
            (inv_tap[l] * model[:v_magnitude][params.bus_from[l], t])^2 -
            inv_tap[l] * model[:v_magnitude][params.bus_from[l], t] *
            model[:v_magnitude][params.bus_to[l], t] *
            (
                -b[l] * cos(
                    model[:theta][params.bus_from[l], t] -
                    model[:theta][params.bus_to[l], t],
                ) +
                g[l] * sin(
                    model[:theta][params.bus_from[l], t] -
                    model[:theta][params.bus_to[l], t],
                )
            )
        )
    )

    @expression(model, q_reverse_flow_branches[l in active, t in 1:T],
        params.s_base * (
            -(b[l] + params.b_shunt[l] / 2) *
            model[:v_magnitude][params.bus_to[l], t]^2 -
            inv_tap[l] * model[:v_magnitude][params.bus_to[l], t] *
            model[:v_magnitude][params.bus_from[l], t] *
            (
                -b[l] * cos(
                    model[:theta][params.bus_from[l], t] -
                    model[:theta][params.bus_to[l], t],
                ) -
                g[l] * sin(
                    model[:theta][params.bus_from[l], t] -
                    model[:theta][params.bus_to[l], t],
                )
            )
        )
    )
end

function create_ac_branch_constraints(model, params)
    T = get(model.ext, :time_periods, 1)
    active_branches = findall(i -> params.existing[i] == 1, 1:params.number_of_branches)

    @variable(model, s_branches[l in active_branches, t = 1:T] >= 0)
    @variable(model, s_reverse_branches[l in active_branches, t = 1:T] >= 0)

    @constraint(model, [l in active_branches, t in 1:T],
        s_branches[l, t]^2 >=
        model[:p_flow_branches][l, t]^2 + model[:q_flow_branches][l, t]^2
    )

    @constraint(model, [l in active_branches, t in 1:T],
        s_branches[l, t] <= params.capacity[l]
    )

    @constraint(model, [l in active_branches, t in 1:T],
        s_reverse_branches[l, t]^2 >=
        model[:p_reverse_flow_branches][l, t]^2 + model[:q_reverse_flow_branches][l, t]^2
    )

    @constraint(model, [l in active_branches, t in 1:T],
        s_reverse_branches[l, t] <= params.capacity[l]
    )

    @constraint(model,
        branch_angle_lb[i in active_branches, t in 1:T; !isnan(params.min_ang[i])],
        params.min_ang[i] <=
        model[:theta][params.bus_from[i], t] - model[:theta][params.bus_to[i], t]
    )

    @constraint(model,
        branch_angle_ub[i in active_branches, t in 1:T; !isnan(params.max_ang[i])],
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

    tap = [isnan(t) || t == 0.0 ? 1.0 : t for t in params.tap_position]
    inv_tap = 1.0 ./ tap

    active = findall(i -> params.existing[i] == 1, 1:params.number_of_branches)

    @variable(model, R[l in active, t = 1:T])
    @variable(model, T_var[l in active, t = 1:T])

    for l in active, t in 1:T
        i_from = params.bus_idx_map[params.bus_from[l]]
        i_to = params.bus_idx_map[params.bus_to[l]]

        Vm_from = params.initial_voltage_magnitude[i_from]
        Vm_to = params.initial_voltage_magnitude[i_to]
        Va_from = params.initial_voltage_angle[i_from]
        Va_to = params.initial_voltage_angle[i_to]

        if !any(isnan.([Vm_from, Vm_to, Va_from, Va_to]))
            set_start_value(R[l, t], Vm_from * Vm_to * cos(Va_from - Va_to))
            set_start_value(T_var[l, t], Vm_from * Vm_to * sin(Va_from - Va_to))
        end
    end

    @expression(model, p_flow_branches[l = 1:params.number_of_branches, t = 1:T],
        params.existing[l] == 1 ?
        params.s_base * (
            g[l] * inv_tap[l]^2 * model[:v_magnitude_sq][params.bus_from[l], t] -
            inv_tap[l] * g[l] * R[l, t] - inv_tap[l] * b[l] * T_var[l, t]
        ) : 0
    )

    @expression(model, p_reverse_flow_branches[l = 1:params.number_of_branches, t = 1:T],
        params.existing[l] == 1 ?
        params.s_base * (
            g[l] * model[:v_magnitude_sq][params.bus_to[l], t] -
            inv_tap[l] * g[l] * R[l, t] + inv_tap[l] * b[l] * T_var[l, t]
        ) : 0
    )

    @expression(model, q_flow_branches[l = 1:params.number_of_branches, t = 1:T],
        params.existing[l] == 1 ?
        params.s_base * (
            -(b[l] + params.b_shunt[l] / 2) * inv_tap[l]^2 *
            model[:v_magnitude_sq][params.bus_from[l], t] +
            inv_tap[l] * b[l] * R[l, t] - inv_tap[l] * g[l] * T_var[l, t]
        ) : 0
    )

    @expression(model, q_reverse_flow_branches[l = 1:params.number_of_branches, t = 1:T],
        params.existing[l] == 1 ?
        params.s_base * (
            -(b[l] + params.b_shunt[l] / 2) * model[:v_magnitude_sq][params.bus_to[l], t] +
            inv_tap[l] * b[l] * R[l, t] + inv_tap[l] * g[l] * T_var[l, t]
        ) : 0
    )
end

function create_socp_branch_constraints(model, params)
    T = get(model.ext, :time_periods, 1)
    active = findall(i -> params.existing[i] == 1, 1:params.number_of_branches)

    @constraint(model, branch_angle_lb[i in active, t = 1:T; !isnan(params.min_ang[i])],
        params.min_ang[i] <=
        model[:theta][params.bus_from[i], t] - model[:theta][params.bus_to[i], t]
    )

    @constraint(model, branch_angle_ub[i in active, t = 1:T; !isnan(params.max_ang[i])],
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
