# ============================================================
# ALADIN-SOCWR — Conic Jacobian assembly
# ============================================================

# ============================================================
# Cone active-face normals
# ============================================================

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

    eq_jac = zone.eq_jac

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
        eq_jac
    else
        C_ineq =
            SparseArrays.sparse(ineq_rows, ineq_cols, ineq_vals, row_offset, n_vars)
        vcat(eq_jac, C_ineq)
    end

    return _delete_dependent_rows(C_stacked; tol = dep_tol)
end
