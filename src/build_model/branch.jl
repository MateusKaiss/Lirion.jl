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

# ---------------------------------------------------------------------
# Creating SOCP Branch Variables, Expressions and Constraints
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# Bounds R/T_var por ramo
# ---------------------------------------------------------------------
function _rt_product_bounds(bus_params, branch_params)
    vmin = bus_params.vmin
    vmax = bus_params.vmax
    angmin = branch_params.min_ang
    angmax = branch_params.max_ang
    bf = branch_params.bus_from
    bt = branch_params.bus_to
    L = branch_params.number_of_branches

    R_min = Dict{Int, Float64}()
    R_max = Dict{Int, Float64}()
    T_min = Dict{Int, Float64}()
    T_max = Dict{Int, Float64}()

    @inbounds for l in 1:L
        i = bf[l]
        j = bt[l]

        vm_fr_min = vmin[i]
        vm_fr_max = vmax[i]
        vm_to_min = vmin[j]
        vm_to_max = vmax[j]
        amin = angmin[l]
        amax = angmax[l]

        R_min[l] = -Inf
        R_max[l] = Inf
        T_min[l] = -Inf
        T_max[l] = Inf

        if amin >= 0
            R_max[l] = min(R_max[l], vm_fr_max * vm_to_max * cos(amin))
            R_min[l] = max(R_min[l], vm_fr_min * vm_to_min * cos(amax))
            T_max[l] = min(T_max[l], vm_fr_max * vm_to_max * sin(amax))
            T_min[l] = max(T_min[l], vm_fr_min * vm_to_min * sin(amin))
        end
        if amax <= 0
            R_max[l] = min(R_max[l], vm_fr_max * vm_to_max * cos(amax))
            R_min[l] = max(R_min[l], vm_fr_min * vm_to_min * cos(amin))
            T_max[l] = min(T_max[l], vm_fr_min * vm_to_min * sin(amax))
            T_min[l] = max(T_min[l], vm_fr_max * vm_to_max * sin(amin))
        end
        if amin < 0 && amax > 0
            R_max[l] = min(R_max[l], vm_fr_max * vm_to_max * 1.0)
            R_min[l] = max(R_min[l], vm_fr_min * vm_to_min * min(cos(amin), cos(amax)))
            T_max[l] = min(T_max[l], vm_fr_max * vm_to_max * sin(amax))
            T_min[l] = max(T_min[l], vm_fr_max * vm_to_max * sin(amin))
        end
    end

    return R_min, R_max, T_min, T_max
end

# ---------------------------------------------------------------------
# Criação de R/T_var com bounds e starts
# ---------------------------------------------------------------------
@inline function _get_bus_value(V, bus_id::Int, bus_idx_map::Dict{Int, Int})
    if V isa AbstractDict
        return get(V, bus_id, nothing)
    else
        return V[bus_idx_map[bus_id]]
    end
end

function create_rt_variables_with_bounds!(
    model,
    matpower_data,
    branch_params;
    bounded::Bool = true,
)
    T = get(model.ext, :time_periods, 1)
    active =
        findall(i -> branch_params.existing[i] == 1, 1:branch_params.number_of_branches)

    @variable(model, R[l in active, t in 1:T], base_name = "R")
    @variable(model, T_var[l in active, t in 1:T], base_name = "T_var")

    Vm = matpower_data["bus_params"].vstart
    Va = matpower_data["bus_params"].initial_voltage_angle
    map = branch_params.bus_idx_map
    shf = branch_params.tap_angle

    @inbounds for l in active, t in 1:T
        bus_i = branch_params.bus_from[l]
        bus_j = branch_params.bus_to[l]

        Vm_i = _get_bus_value(Vm, bus_i, map)
        Vm_j = _get_bus_value(Vm, bus_j, map)
        Va_i = _get_bus_value(Va, bus_i, map)
        Va_j = _get_bus_value(Va, bus_j, map)

        if !(
            Vm_i === nothing || Vm_j === nothing || Va_i === nothing || Va_j === nothing
        ) &&
           !any(isnan.((Vm_i, Vm_j, Va_i, Va_j, shf[l])))
            set_start_value(R[l, t], Vm_i * Vm_j * cos(Va_i - Va_j - shf[l]))
            set_start_value(T_var[l, t], Vm_i * Vm_j * sin(Va_i - Va_j - shf[l]))
        end
    end

    if bounded
        R_min, R_max, T_min, T_max =
            _rt_product_bounds(matpower_data["bus_params"], branch_params)
        @inbounds for l in active, t in 1:T
            isfinite(R_min[l]) && set_lower_bound(R[l, t], R_min[l])
            isfinite(R_max[l]) && set_upper_bound(R[l, t], R_max[l])
            isfinite(T_min[l]) && set_lower_bound(T_var[l, t], T_min[l])
            isfinite(T_max[l]) && set_upper_bound(T_var[l, t], T_max[l])
        end
    end

    model[:R] = R
    model[:T_var] = T_var
    return nothing
end

# ---------------------------------------------------------------------
# Tangent cuts: T_var <= tan(angmax)*R  e  T_var >= tan(angmin)*R
# ---------------------------------------------------------------------
function add_rt_tangent_angle_cuts!(model, branch_params)
    T = get(model.ext, :time_periods, 1)
    active =
        findall(i -> branch_params.existing[i] == 1, 1:branch_params.number_of_branches)

    R = model[:R]
    T_var = model[:T_var]

    @constraint(model, [l in active, t in 1:T; !isnan(branch_params.max_ang[l])],
        T_var[l, t] <= tan(branch_params.max_ang[l]) * R[l, t]
    )
    @constraint(model, [l in active, t in 1:T; !isnan(branch_params.min_ang[l])],
        T_var[l, t] >= tan(branch_params.min_ang[l]) * R[l, t]
    )
    return nothing
end

function compute_socp_branch_expressions(model, params)
    T = get(model.ext, :time_periods, 1)
    active = findall(i -> params.existing[i] == 1, 1:params.number_of_branches)
    inv_tap = [1.0 / (isnan(t) || t == 0.0 ? 1.0 : t) for t in params.tap_position]

    R = model[:R]
    T_var = model[:T_var]

    tap = params.tap_position
    shift = params.tap_angle
    tr = tap .* cos.(shift)
    ti = tap .* sin.(shift)

    @expression(model, p_flow_branches[l in active, t in 1:T],
        params.g[l] * inv_tap[l]^2 * model[:v_magnitude_sq][params.bus_from[l], t] +
        (-params.g[l] * tr[l] + params.b[l] * ti[l]) * inv_tap[l]^2 * R[l, t] +
        (-params.b[l] * tr[l] - params.g[l] * ti[l]) * inv_tap[l]^2 * T_var[l, t]
    )
    @expression(model, q_flow_branches[l in active, t in 1:T],
        -(params.b[l] + params.b_shunt[l] / 2) * inv_tap[l]^2 *
        model[:v_magnitude_sq][params.bus_from[l], t] -
        (-params.b[l] * tr[l] - params.g[l] * ti[l]) * inv_tap[l]^2 * R[l, t] +
        (-params.g[l] * tr[l] + params.b[l] * ti[l]) * inv_tap[l]^2 * T_var[l, t]
    )

    @expression(model, p_reverse_flow_branches[l in active, t in 1:T],
        params.g[l] * model[:v_magnitude_sq][params.bus_to[l], t] +
        (-params.g[l] * tr[l] - params.b[l] * ti[l]) * inv_tap[l]^2 * R[l, t] +
        (+params.b[l] * tr[l] - params.g[l] * ti[l]) * inv_tap[l]^2 * T_var[l, t]
    )
    @expression(model, q_reverse_flow_branches[l in active, t in 1:T],
        -(params.b[l] + params.b_shunt[l] / 2) *
        model[:v_magnitude_sq][params.bus_to[l], t] +
        (+params.b[l] * tr[l] - params.g[l] * ti[l]) * inv_tap[l]^2 * R[l, t] +
        (+params.g[l] * tr[l] + params.b[l] * ti[l]) * inv_tap[l]^2 * T_var[l, t]   # = (-g*tr-b*ti)*(1/|t|^2)*(-T_var)
    )
    model[:p_flow_branches] = p_flow_branches
    model[:q_flow_branches] = q_flow_branches
    model[:p_reverse_flow_branches] = p_reverse_flow_branches
    model[:q_reverse_flow_branches] = q_reverse_flow_branches
    return nothing
end

# ---------------------------------------------------------------------
# Constraints (inclui RSOC, limites térmicos e cortes tangentes)
# ---------------------------------------------------------------------
function create_socp_branch_constraints(model, params)
    T = get(model.ext, :time_periods, 1)
    active = findall(i -> params.existing[i] == 1, 1:params.number_of_branches)

    # ângulo em theta (mantém seu padrão)
    @constraint(model, branch_angle_lb[i in active, t in 1:T; !isnan(params.min_ang[i])],
        params.min_ang[i] <=
        model[:theta][params.bus_from[i], t] - model[:theta][params.bus_to[i], t]
    )
    @constraint(model, branch_angle_ub[i in active, t in 1:T; !isnan(params.max_ang[i])],
        model[:theta][params.bus_from[i], t] - model[:theta][params.bus_to[i], t] <=
        params.max_ang[i]
    )

    # cone rotacionado: wr^2+wi^2 ≤ w_i*w_j  ⇒  [w_i/√2, w_j/√2, R, T] ∈ RSOC
    @constraint(model, [l in active, t in 1:T],
        [
            model[:v_magnitude_sq][params.bus_from[l], t] / sqrt(2),
            model[:v_magnitude_sq][params.bus_to[l], t] / sqrt(2),
            model[:R][l, t],
            model[:T_var][l, t],
        ] in JuMP.RotatedSecondOrderCone()
    )

    # limites térmicos (SOCP)
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

    # cortes de ângulo em (R,T_var) estilo PM
    add_rt_tangent_angle_cuts!(model, params)

    return nothing
end

# ------------------------------------------------------------------------
# Dispatching Branch Variables, Expressions e Constraints por tipo
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

    # 1) criar R/T_var com bounds e starts
    create_rt_variables_with_bounds!(model, matpower_data, params)

    # 2) expressões (usam R/T_var criados)
    compute_socp_branch_expressions(model, params)

    # 3) constraints (RSOC, térmicos e cortes tangentes)
    create_socp_branch_constraints(model, params)
    return nothing
end
