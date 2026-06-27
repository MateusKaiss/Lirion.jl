# ============================================================
# ALADIN-AC — Per-iteration sensitivities (MOI.Nonlinear.Evaluator)
# ============================================================

# ============================================================
# AC cost expression
# ============================================================

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

# ============================================================
# Active-constraint NLP + evaluator
# ============================================================

@inline function _ac_keep_constraint(
    ci,
    ::Type{S},
    has_duals::Bool,
    dual_tol::Float64,
) where {S}
    S === MOI.EqualTo{Float64} && return true
    has_duals || return false
    dval = try
        JuMP.dual(ci)
    catch
        0.0
    end
    return abs(dval) > dual_tol
end

function _build_ac_evaluator(zone::AladinZoneDataAC, dual_tol::Float64)
    model = zone.local_model
    has_duals = JuMP.has_duals(model)

    nlp = MOI.Nonlinear.Model()
    rows = JuMP.ConstraintRef[]

    for (F, S) in JuMP.list_of_constraint_types(model)
        F <: JuMP.VariableRef && continue
        for ci in JuMP.all_constraints(model, F, S)
            _ac_keep_constraint(ci, S, has_duals, dual_tol) || continue
            cobj = JuMP.constraint_object(ci)
            MOI.Nonlinear.add_constraint(nlp, cobj.func, cobj.set)
            push!(rows, ci)
        end
    end

    obj_expr = _build_ac_cost_expr(zone)
    obj_expr === nothing || MOI.Nonlinear.set_objective(nlp, obj_expr)

    ordered_vars = [JuMP.index(v) for v in zone.all_vars]
    evaluator = MOI.Nonlinear.Evaluator(
        nlp,
        MOI.Nonlinear.SparseReverseMode(),
        ordered_vars,
    )
    MOI.initialize(evaluator, [:Jac, :Hess])
    return evaluator, rows
end

# ============================================================
# Active-constraint Jacobian
# ============================================================

function _extract_active_jacobian_ac(
    evaluator,
    rows::Vector{JuMP.ConstraintRef},
    y::Vector{Float64},
    n_vars::Int,
    dep_tol::Float64,
)::SparseArrays.SparseMatrixCSC{Float64, Int}
    n_rows = length(rows)
    n_rows == 0 && return SparseArrays.spzeros(Float64, 0, n_vars)

    J_struct = MOI.jacobian_structure(evaluator)
    V_J = zeros(length(J_struct))
    MOI.eval_constraint_jacobian(evaluator, V_J, y)

    I_J = [t[1] for t in J_struct]
    C_J = [t[2] for t in J_struct]
    C_full = SparseArrays.sparse(I_J, C_J, V_J, n_rows, n_vars)
    return _delete_dependent_rows(C_full; tol = dep_tol)
end

# ============================================================
# Lagrangian Hessian (symmetrised from lower triangle)
# ============================================================

function _extract_lagrangian_hessian_ac(
    evaluator,
    rows::Vector{JuMP.ConstraintRef},
    y::Vector{Float64},
    n_vars::Int,
)::SparseArrays.SparseMatrixCSC{Float64, Int}
    H_struct = MOI.hessian_lagrangian_structure(evaluator)
    V_H = zeros(length(H_struct))
    y_duals = Float64[-JuMP.dual(ci) for ci in rows]
    if !isempty(H_struct)
        MOI.eval_hessian_lagrangian(evaluator, V_H, y, 1.0, y_duals)
    end

    I_H = [t[1] for t in H_struct]
    C_H = [t[2] for t in H_struct]
    H_lower = SparseArrays.sparse(I_H, C_H, V_H, n_vars, n_vars)

    return H_lower + SparseArrays.transpose(H_lower) -
           SparseArrays.spdiagm(0 => LinearAlgebra.diag(H_lower))
end

# ============================================================
# Diagonal regularisation: proximal ρ·σ + Tikhonov
# ============================================================

function _add_prox_curvature_ac!(
    H::SparseArrays.SparseMatrixCSC{Float64, Int},
    zone::AladinZoneDataAC,
    rho::Float64;
    sigma_theta::Float64,
    sigma_vm::Float64,
    sigma_pg::Float64,
    sigma_qg::Float64,
    sigma_other::Float64,
)
    _accumulate_prox_diagonal!(
        rho,
        (
            (zone.theta_cols, sigma_theta),
            (zone.vm_cols, sigma_vm),
            (zone.pg_cols, sigma_pg),
            (zone.qg_cols, sigma_qg),
            (zone.curt_cols, sigma_other),
        ),
    ) do c, s
        return H[c, c] += s
    end
    return nothing
end

@inline function _add_tikhonov_ac!(
    H::SparseArrays.SparseMatrixCSC{Float64, Int},
    n_vars::Int;
    delta::Float64 = 1e-6,
)
    @inbounds for c in 1:n_vars
        H[c, c] += delta
    end
    @assert all(isfinite, H.nzval) "AC Hessian has non-finite entries"
    return nothing
end

# ============================================================
# Orchestrator
# ============================================================

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
    y = zone.y
    n_vars = zone.n_vars

    evaluator, rows = _build_ac_evaluator(zone, dual_tol)

    zone.active_jac = _extract_active_jacobian_ac(evaluator, rows, y, n_vars, dep_tol)
    H = _extract_lagrangian_hessian_ac(evaluator, rows, y, n_vars)
    
    _add_prox_curvature_ac!(
        H, zone, rho;
        sigma_theta = sigma_theta,
        sigma_vm = sigma_vm,
        sigma_pg = sigma_pg,
        sigma_qg = sigma_qg,
        sigma_other = sigma_other,
    )
    n_vars > 0 && _add_tikhonov_ac!(H, n_vars)
    zone.hess = H
    return nothing
end

# ============================================================
# Sensitivity refresh hook
# ============================================================

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
