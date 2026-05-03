# ============================================================
# ALADIN-AC — Per-iteration sensitivities
# ============================================================
#
# Hot-path kernels refreshed every ALADIN iteration:
#   - `extract_solve_values!` snapshots Phase-1 NLP solution into
#     `zone.y`;
#   - `compute_cost_gradient!` writes ∇f(y) into `zone.g`;
#   - `refresh_sensitivities_ac!` rebuilds `zone.hess` and
#     `zone.jac_cache` using a `MOI.Nonlinear.Evaluator` (AC
#     equalities/inequalities are non-linear so the Jacobian is not
#     constant; the Hessian needs a Levenberg-Marquardt floor to
#     enforce H_i + ρΣ_i ≻ 0 of HFD Theorem 4);
#   - `refresh_lagrangian_hessian!` is the public hook used by the
#     coordinator pipeline.
#
# Struct definition and static initialisation live in `zone_data.jl`;
# active-set logic lives in `active_set.jl`.

# ============================================================
# Phase-1 → buffer
# ============================================================

function extract_solve_values!(zone::AladinZoneDataAC, model::JuMP.Model)
    y = zone.y
    @inbounds for i in 1:zone.n_vars
        y[i] = JuMP.value(zone.all_vars[i])
    end
    return nothing
end

# ============================================================
# Cost gradient ∇f(y)
# ============================================================

function compute_cost_gradient!(zone::AladinZoneDataAC)
    g = zone.g
    fill!(g, 0.0)
    y = zone.y
    @inbounds for k in eachindex(zone.pg_cols)
        col = zone.pg_cols[k]
        g[col] = zone.pg_c2_sbase2[k] * y[col] + zone.pg_c1_sbase[k]
    end
    return nothing
end

# ============================================================
# Hessian + Jacobian refresh via MOI.Nonlinear.Evaluator
# ============================================================

function refresh_sensitivities_ac!(
    zone::AladinZoneDataAC;
    rho::Float64 = 0.0,
    sigma_theta::Float64 = 0.0,
    sigma_vm::Float64 = 0.0,
    sigma_pg::Float64 = 0.0,
    sigma_qg::Float64 = 0.0,
    sigma_other::Float64 = 0.0,
    # AC active-set tolerance is intentionally looser (1e-3) than DC/SOCWR
    # (1e-6). AC subproblems are non-convex NLPs solved by Ipopt, whose
    # multipliers are typically less precise than those returned by QP/SOCP
    # solvers used in DC/SOCWR. A tighter threshold flags too many
    # constraints as active and breaks the coordinator QP. Documented as
    # canonical-equivalent — model-specific tolerances are standard practice.
    dual_tol::Float64 = 1e-3,
    dep_tol::Float64 = 1e-10,
)
    model = zone.model
    n_vars = zone.n_vars
    has_duals = JuMP.has_duals(model)

    nlp = MOI.Nonlinear.Model()
    rows = JuMP.ConstraintRef[]   # real constraints corresponding to aux NLP rows

    # ---- Add active constraints to the aux NLP ----
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

    # ---- Set objective = generation cost ----
    obj_expr = _build_ac_cost_expr(zone)
    if obj_expr !== nothing
        MOI.Nonlinear.set_objective(nlp, obj_expr)
    end

    # ---- Build the evaluator once, initialise for Jac + Hess ----
    ordered_vars = [JuMP.index(v) for v in zone.all_vars]
    evaluator = MOI.Nonlinear.Evaluator(
        nlp,
        MOI.Nonlinear.SparseReverseMode(),
        ordered_vars,
    )
    MOI.initialize(evaluator, [:Jac, :Hess])

    y_primal = zone.y

    # ---- Jacobian ----
    J_struct = MOI.jacobian_structure(evaluator)
    V_J = zeros(length(J_struct))
    MOI.eval_constraint_jacobian(evaluator, V_J, y_primal)
    n_rows = length(rows)
    if n_rows == 0
        zone.jac_cache = SparseArrays.spzeros(Float64, 0, n_vars)
    else
        I_J = [t[1] for t in J_struct]
        C_J = [t[2] for t in J_struct]
        C_full = SparseArrays.sparse(I_J, C_J, V_J, n_rows, n_vars)
        zone.jac_cache = _delete_dependent_rows(C_full; tol = dep_tol)
    end

    # ---- Hessian of Lagrangian ----
    H_struct = MOI.hessian_lagrangian_structure(evaluator)
    V_H = zeros(length(H_struct))
    y_duals = Float64[-JuMP.dual(ci) for ci in rows]
    if !isempty(H_struct)
        MOI.eval_hessian_lagrangian(evaluator, V_H, y_primal, 1.0, y_duals)
    end

    I_H = [t[1] for t in H_struct]
    C_H = [t[2] for t in H_struct]
    H_lower = SparseArrays.sparse(I_H, C_H, V_H, n_vars, n_vars)
    # MOI returns only the lower triangle — symmetrise explicitly.
    H_full =
        H_lower + SparseArrays.transpose(H_lower) -
        SparseArrays.spdiagm(0 => LinearAlgebra.diag(H_lower))

    if rho > 0.0
        add! = (col_set, σ) -> begin
            s = rho * σ
            s > 0.0 || return
            for c in col_set
                H_full[c, c] += s
            end
        end
        add!(zone.theta_cols, sigma_theta)
        add!(zone.vm_cols, sigma_vm)
        add!(zone.pg_cols, sigma_pg)
        add!(zone.qg_cols, sigma_qg)
        add!(zone.curt_cols, sigma_other)
    end

    # Levenberg-Marquardt regularisation: add δ on the diagonal so the
    # Hessian satisfies the H_i + ρΣ_i ≻ 0 hypothesis of HFD Theorem 4.
    # This replaces the previous dense eigendecomposition (O(n³)/iter) with
    # an O(n) shift; the perturbation is bounded by δ in operator norm and
    # matches SOCWR's convention (`_refresh_hess_diag!`, δ=1e-6).
    if n_vars > 0
        δ = 1e-6
        @inbounds for c in 1:n_vars
            H_full[c, c] += δ
        end
        # Sanity check: NaN/Inf in the Hessian indicates a numerical issue
        # upstream (Ipopt failure, tap_angle wrap, etc.) — fail loudly here
        # rather than propagate garbage into the coordinator QP.
        @assert all(isfinite, H_full.nzval) "AC Hessian has non-finite entries"
    end

    zone.hess = H_full
    return nothing
end

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
