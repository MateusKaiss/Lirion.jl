# ============================================================
# SOCWRNL — smooth-NLP form of the SOCWR W-space relaxation
# ============================================================
#   SOCWR     <-> SOCWRConicPowerModel  (MOI cones, conic solver)
#   SOCWRNL   <-> SOCWRPowerModel       (quadratic ineqs, NLP solver)

function add_socwrnl_branch_constraints!(
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

    # one Jabr row per bus pair; smooth-NLP form of the RotatedSecondOrderCone
    #   [w_i/√2, w_j/√2, R, T] ∈ RSOC  ⟺  2·(w_i/√2)·(w_j/√2) ≥ R²+T²
    #                                   ⟺  w_i·w_j ≥ R² + T²
    # written as the quadratic inequality  R² + T² − w_i·w_j ≤ 0.
    pair_reps = get(model.ext, :socwr_pair_reps, active)
    @constraint(model, jabr_nl[l in pair_reps, t in 1:T],
        model[:R][l, t]^2 + model[:T_var][l, t]^2 -
        model[:v_magnitude_sq][br.from[l], t] * model[:v_magnitude_sq][br.to[l], t] <= 0
    )

    # Thermal limits, smooth-NLP form of  [rate_a, p, q] ∈ SecondOrderCone
    #   rate_a ≥ ‖(p,q)‖  ⟺  p² + q² ≤ rate_a²
    # written as  p² + q² − rate_a² ≤ 0 (rate_a is a data constant).
    @constraint(model, thermal_fwd_nl[l in active, t in 1:T],
        model[:p_flow_branches][l, t]^2 + model[:q_flow_branches][l, t]^2 -
        br.rate_a[l]^2 <= 0
    )

    @constraint(model, thermal_rev_nl[l in active, t in 1:T],
        model[:p_reverse_flow_branches][l, t]^2 +
        model[:q_reverse_flow_branches][l, t]^2 -
        br.rate_a[l]^2 <= 0
    )

    add_rt_tangent_angle_cuts!(model, br)
    add_socwr_lnc_cuts!(model, bus, br)

    return nothing
end

# ============================================================
# SOCWRNL dispatch — reuse the SOCWR W-space machinery verbatim
# ============================================================

add_bus_variables!(model::JuMP.Model, bus::BusData, ::Type{SOCWRNL}) =
    add_bus_variables!(model, bus, SOCWR)

add_branch_variables!(model::JuMP.Model, sys::SystemData, ::Type{SOCWRNL}) =
    add_branch_variables!(model, sys, SOCWR)

add_branch_expressions!(model::JuMP.Model, sys::SystemData, ::Type{SOCWRNL}) =
    add_branch_expressions!(model, sys, SOCWR)

# The ONLY stage that differs: quadratic inequalities instead of MOI cones.
add_branch_constraints!(model::JuMP.Model, sys::SystemData, ::Type{SOCWRNL}) =
    add_socwrnl_branch_constraints!(model, sys.bus, sys.branch)

add_generator_variables!(model::JuMP.Model, gen::GenData, ::Type{SOCWRNL}) =
    add_generator_variables!(model, gen, SOCWR)

add_generator_constraints!(model::JuMP.Model, gen::GenData, ::Type{SOCWRNL}) =
    add_generator_constraints!(model, gen, SOCWR)

add_storage_constraints!(model::JuMP.Model, stor::StorageData, ::Type{SOCWRNL}) =
    add_storage_constraints!(model, stor, SOCWR)

add_active_power_balance_constraints!(
    model::JuMP.Model, sys::SystemData, ::Type{SOCWRNL},
) = add_active_power_balance_constraints!(model, sys, SOCWR)

add_reactive_power_balance_constraints!(
    model::JuMP.Model, sys::SystemData, ::Type{SOCWRNL},
) = add_reactive_power_balance_constraints!(model, sys, SOCWR)

add_objective_expressions!(model::JuMP.Model, obj::ObjectiveData, ::Type{SOCWRNL}) =
    add_objective_expressions!(model, obj, SOCWR)
