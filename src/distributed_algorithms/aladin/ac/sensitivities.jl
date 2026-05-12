# ============================================================
# ALADIN-AC — Per-iteration sensitivities
# ============================================================

# ------------------------------------------------------------
# Hessian + Jacobian refresh via MOI.Nonlinear.Evaluator
# ------------------------------------------------------------

function refresh_sensitivities_ac!(
    zone::AladinZoneDataAC;
    rho::Float64 = 0.0,
    sigma_theta::Float64 = 0.0,
    sigma_vm::Float64 = 0.0,
    sigma_pg::Float64 = 0.0,
    sigma_qg::Float64 = 0.0,
    sigma_other::Float64 = 0.0,

    dual_tol::Float64 = 1e-3,
    dep_tol::Float64 = 1e-10,
)
    model = zone.model
    n_vars = zone.n_vars
    has_duals = JuMP.has_duals(model)

    nlp = MOI.Nonlinear.Model()
    rows = JuMP.ConstraintRef[]

    for (F, S) in JuMP.list_of_constraint_types(model)
        F <: JuMP.VariableRef && continue

        for ci in JuMP.all_constraints(model, F, S)
            keep = if S === MOI.EqualTo{Float64}
                true
            elseif has_duals
                dval = try
                    JuMP.dual(ci)
                catch
                    0.0
                end
                abs(dval) > dual_tol
            else
                false
            end
            keep || continue

            cobj = JuMP.constraint_object(ci)
            MOI.Nonlinear.add_constraint(nlp, cobj.func, cobj.set)
            push!(rows, ci)
        end
    end

    obj_expr = _build_ac_cost_expr(zone)
    if obj_expr !== nothing
        MOI.Nonlinear.set_objective(nlp, obj_expr)
    end

    ordered_vars = [JuMP.index(v) for v in zone.all_vars]
    evaluator = MOI.Nonlinear.Evaluator(
        nlp,
        MOI.Nonlinear.SparseReverseMode(),
        ordered_vars,
    )
    MOI.initialize(evaluator, [:Jac, :Hess])

    y_primal = zone.y

    J_struct = MOI.jacobian_structure(evaluator)
    V_J = zeros(length(J_struct))
    MOI.eval_constraint_jacobian(evaluator, V_J, y_primal)
    n_rows = length(rows)
    if n_rows == 0
        zone.active_jac = SparseArrays.spzeros(Float64, 0, n_vars)
    else
        I_J = [t[1] for t in J_struct]
        C_J = [t[2] for t in J_struct]
        C_full = SparseArrays.sparse(I_J, C_J, V_J, n_rows, n_vars)
        zone.active_jac = _delete_dependent_rows(C_full; tol = dep_tol)
    end

    H_struct = MOI.hessian_lagrangian_structure(evaluator)
    V_H = zeros(length(H_struct))
    y_duals = Float64[-JuMP.dual(ci) for ci in rows]
    if !isempty(H_struct)
        MOI.eval_hessian_lagrangian(evaluator, V_H, y_primal, 1.0, y_duals)
    end

    I_H = [t[1] for t in H_struct]
    C_H = [t[2] for t in H_struct]
    H_lower = SparseArrays.sparse(I_H, C_H, V_H, n_vars, n_vars)

    H_full =
        H_lower + SparseArrays.transpose(H_lower) -
        SparseArrays.spdiagm(0 => LinearAlgebra.diag(H_lower))

    if rho > 0.0
        add_rho_sigma! = (col_set, sigma) -> begin
            s = rho * sigma
            s > 0.0 || return
            for c in col_set
                H_full[c, c] += s
            end
        end
        add_rho_sigma!(zone.theta_cols, sigma_theta)
        add_rho_sigma!(zone.vm_cols, sigma_vm)
        add_rho_sigma!(zone.pg_cols, sigma_pg)
        add_rho_sigma!(zone.qg_cols, sigma_qg)
        add_rho_sigma!(zone.curt_cols, sigma_other)
    end

    if n_vars > 0
        delta = 1e-6
        @inbounds for c in 1:n_vars
            H_full[c, c] += delta
        end
        @assert all(isfinite, H_full.nzval) "AC Hessian has non-finite entries"
    end

    zone.hess = H_full
    return nothing
end

# ------------------------------------------------------------
# Sensitivity refresh hook
# ------------------------------------------------------------

@inline function refresh_lagrangian_hessian!(
    zone::AladinZoneDataAC,
    ws::AbstractALADINWorkspace,
)
    refresh_sensitivities_ac!(
        zone;
        rho = ws.rho[],
        sigma_theta = ws.sigma_theta,
        sigma_vm = ws.sigma_vm,
        sigma_pg = ws.sigma_pg,
        sigma_qg = ws.sigma_qg,
        sigma_other = ws.sigma_other,
    )
    return nothing
end

# ------------------------------------------------------------
# AC cost expression
# ------------------------------------------------------------

function _build_ac_cost_expr(zone::AladinZoneDataAC)
    parts = Any[]
    @inbounds for k in eachindex(zone.pg_cols)
        col = zone.pg_cols[k]
        vi = JuMP.index(zone.all_vars[col])
        c2 = zone.pg_c2_sbase2[k] / 2.0
        c1 = zone.pg_c1_sbase[k]
        if c2 != 0.0
            push!(parts, :($c2 * $vi * $vi))
        end
        if c1 != 0.0
            push!(parts, :($c1 * $vi))
        end
    end
    isempty(parts) && return nothing
    length(parts) == 1 && return parts[1]
    return Expr(:call, :+, parts...)
end
