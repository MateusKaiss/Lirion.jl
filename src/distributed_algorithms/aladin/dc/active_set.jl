# ============================================================
# ALADIN-DC — Active-constraint Jacobian
# ============================================================

function build_active_jacobian(
    zone::AladinZoneData;
    dual_tol::Float64 = 1e-6,
    dep_tol::Float64 = 1e-10,
)::SparseArrays.SparseMatrixCSC{Float64, Int}
    model = zone.model
    var_to_col = zone.var_to_col
    n_vars = zone.n_vars

    eq_jac = zone.eq_jac
    eq_nrows = size(eq_jac, 1)

    ineq_rows = Int[]
    ineq_cols = Int[]
    ineq_vals = Vector{Float64}()
    row_offset = 0

    has_duals = JuMP.has_duals(model)

    for (F, S) in JuMP.list_of_constraint_types(model)
        F <: JuMP.VariableRef && continue
        S === MOI.EqualTo{Float64} && continue

        for ci in JuMP.all_constraints(model, F, S)
            if has_duals
                dval = try
                    JuMP.dual(ci)
                catch
                    0.0
                end
                abs(dval) > dual_tol || continue
            else
                continue
            end

            cobj = JuMP.constraint_object(ci)
            row_offset += 1
            _push_affine_row_offset!(
                ineq_rows, ineq_cols, ineq_vals,
                row_offset, cobj.func, var_to_col,
            )
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

@inline function _push_affine_row_offset!(
    rows::Vector{Int},
    cols::Vector{Int},
    vals::Vector{Float64},
    row::Int,
    func::JuMP.AffExpr,
    var_to_col::Dict{JuMP.VariableRef, Int},
)
    for (v, coef) in func.terms
        col = get(var_to_col, v, 0)
        col == 0 && continue
        push!(rows, row)
        push!(cols, col)
        push!(vals, coef)
    end
    return nothing
end

@inline _push_affine_row_offset!(
    ::Vector{Int},
    ::Vector{Int},
    ::Vector{Float64},
    ::Int,
    ::Any,
    ::Dict,
) = nothing
