# ============================================================
# ALADIN-SOCWR — Active-set detection & conic Jacobian assembly
# ============================================================
#
# Implements the canonical HFD 2016 Algorithm 1 active-set treatment
# for conic problems:
#
#   - bound detection on box constraints (`detect_active_bounds!`);
#   - scalar affine inequalities walked via dual threshold;
#   - vector cone constraints (RSOC, SOC) walked via dual norm
#     threshold; at each active cone we compute the outward normal
#     in the cone's image space and pull it back to the zone's
#     variable layout via the affine function's Jacobian.
#
# The full active Jacobian `[C_eq ; C_ineq_active ; C_cone_active]`
# is rank-reduced via `_delete_dependent_rows` (defined in
# `aladin/dc/active_set.jl` and reused here).

# ============================================================
# Active bounds at current iterate
# ============================================================

function detect_active_bounds!(zone::AladinZoneDataSOCWR; tol::Float64 = 1e-6)
    active = zone.active_bound_cols
    empty!(active)
    y = zone.y
    @inbounds for i in 1:zone.n_vars
        lo = zone.lb[i]
        hi = zone.ub[i]
        v = y[i]
        if !isnan(lo) && v - lo < tol
            push!(active, i)
        elseif !isnan(hi) && hi - v < tol
            push!(active, i)
        end
    end
    return nothing
end

# ============================================================
# Cone active-face normals (canonical HFD Alg. 1 for conic problems)
# ============================================================
# At an active (primal on the cone boundary) conic constraint
# `func(y) ∈ K`, HFD 2016 Algorithm 1 requires preserving the
# active-face hyperplane in the Phase-2 coordinator QP via a row
# `nᵀ · ∇func · Δy = 0`, where n is the outward normal of K at
# `func(y*)` and ∇func is the Jacobian of the affine constraint
# function w.r.t. the zone variables.
#
# We derive n from the squared-slack representation of each cone
# (so both n and its scaling are unambiguous):
#
#   Rotated SOC  K_rot = {(u, v, t₁, t₂, ...) : 2uv ≥ Σ tₖ², u,v ≥ 0}
#     slack g = Σ tₖ² - 2uv  ≤ 0
#     ∇g = (-2v, -2u, 2t₁, 2t₂, ...)
#
#   Standard SOC K_soc = {(t, x₁, x₂, ...) : t² ≥ Σ xₖ², t ≥ 0}
#     slack g = Σ xₖ² - t²   ≤ 0
#     ∇g = (-2t, 2x₁, 2x₂, ...)
#
# For affine functions in cones the row is constant-direction in
# variable space but its location in the cone depends on y*, so it
# must be refreshed every iteration (unlike the equality Jacobian
# C_eq which is structurally constant).

@inline function _cone_active_normal(
    ::MOI.RotatedSecondOrderCone,
    v::AbstractVector{Float64},
)
    k = length(v)
    n = zeros(Float64, k)
    n[1] = -2.0 * v[2]
    n[2] = -2.0 * v[1]
    @inbounds for i in 3:k
        n[i] = 2.0 * v[i]
    end
    return n
end

@inline function _cone_active_normal(
    ::MOI.SecondOrderCone,
    v::AbstractVector{Float64},
)
    k = length(v)
    n = zeros(Float64, k)
    n[1] = -2.0 * v[1]
    @inbounds for i in 2:k
        n[i] = 2.0 * v[i]
    end
    return n
end

@inline function _eval_aff_at_y(
    f::JuMP.AffExpr,
    y::Vector{Float64},
    var_to_col::Dict{JuMP.VariableRef, Int},
)
    acc = f.constant
    for (v, coef) in f.terms
        col = get(var_to_col, v, 0)
        col == 0 && continue
        acc += coef * y[col]
    end
    return acc
end

# ============================================================
# Active-constraint Jacobian
# ============================================================
# C_active = [ C_eq ; C_ineq_active ; C_cone_active ]
#
# SOCWR constraint categories:
#   * affine equalities (power balance, curtailment links) → C_eq
#     (precomputed at build time, structurally constant)
#   * affine scalar inequalities (branch angle bounds)     → walked
#     via scalar dual threshold
#   * vector-valued cone constraints (Jabr RSOC, thermal SOC) → walked
#     via cone dual norm threshold; at each active cone we compute the
#     outward normal in the cone's image space and pull it back to the
#     zone's variable layout via the affine function's Jacobian.
#
# Inactive cones contribute zero to the Lagrangian Hessian (the affine
# functions have zero Hessian) and zero rows here; active cones
# contribute one row per constraint (the tangent hyperplane).

function build_active_jacobian(
    zone::AladinZoneDataSOCWR;
    dual_tol::Float64 = 1e-6,
    cone_tol::Float64 = 1e-6,
    dep_tol::Float64 = 1e-10,
)::SparseArrays.SparseMatrixCSC{Float64, Int}
    model = zone.model
    var_to_col = zone.var_to_col
    n_vars = zone.n_vars
    y = zone.y

    C_eq = zone.C_eq

    ineq_rows = Int[]
    ineq_cols = Int[]
    ineq_vals = Float64[]
    row_offset = 0

    has_duals = JuMP.has_duals(model)

    for (F, S) in JuMP.list_of_constraint_types(model)
        F <: JuMP.VariableRef && continue

        # --- Scalar affine inequality (angle bounds, etc.) ---
        if F <: JuMP.AffExpr
            S === MOI.EqualTo{Float64} && continue
            for ci in JuMP.all_constraints(model, F, S)
                has_duals || continue
                dval = try
                    JuMP.dual(ci)
                catch
                    0.0
                end
                (dval isa Real) || continue
                abs(dval) > dual_tol || continue

                cobj = JuMP.constraint_object(ci)
                row_offset += 1
                _push_affine_row_offset!(
                    ineq_rows, ineq_cols, ineq_vals,
                    row_offset, cobj.func, var_to_col,
                )
            end
            continue
        end

        # --- Vector cone constraints (RSOC, SOC) ---
        if F <: AbstractVector
            is_rsoc = S <: MOI.RotatedSecondOrderCone
            is_soc = S <: MOI.SecondOrderCone
            (is_rsoc || is_soc) || continue

            for ci in JuMP.all_constraints(model, F, S)
                has_duals || continue
                dval_vec = try
                    JuMP.dual(ci)
                catch
                    nothing
                end
                (dval_vec isa AbstractVector) || continue
                isempty(dval_vec) && continue
                LinearAlgebra.norm(dval_vec) > cone_tol || continue

                cobj = JuMP.constraint_object(ci)
                func = cobj.func
                cone_set = cobj.set

                m = length(func)
                v_star = Vector{Float64}(undef, m)
                ok = true
                @inbounds for k in 1:m
                    fk = func[k]
                    if fk isa JuMP.AffExpr
                        v_star[k] = _eval_aff_at_y(fk, y, var_to_col)
                    else
                        ok = false
                        break
                    end
                end
                ok || continue

                nrm = _cone_active_normal(cone_set, v_star)

                row_offset += 1
                @inbounds for k in 1:m
                    nk = nrm[k]
                    nk == 0.0 && continue
                    fk = func[k]
                    fk isa JuMP.AffExpr || continue
                    for (v, coef) in fk.terms
                        col = get(var_to_col, v, 0)
                        col == 0 && continue
                        push!(ineq_rows, row_offset)
                        push!(ineq_cols, col)
                        push!(ineq_vals, nk * coef)
                    end
                end
            end
        end
    end

    C_stacked = if row_offset == 0
        C_eq
    else
        C_ineq =
            SparseArrays.sparse(ineq_rows, ineq_cols, ineq_vals, row_offset, n_vars)
        vcat(C_eq, C_ineq)
    end

    return _delete_dependent_rows(C_stacked; tol = dep_tol)
end
