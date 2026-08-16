# ============================================================
# ALADIN — Zone-data helpers shared across DC / AC / SOCWR
# ============================================================

# ============================================================
# Variable layout
# ============================================================

function _build_var_layout(model::JuMP.Model)
    all_vars = JuMP.all_variables(model)
    n = length(all_vars)
    var_to_col = Dict{JuMP.VariableRef, Int}()
    sizehint!(var_to_col, n)
    @inbounds for i in 1:n
        var_to_col[all_vars[i]] = i
    end
    return all_vars, n, var_to_col
end

function _init_from_start_values(all_vars::Vector{JuMP.VariableRef})
    n = length(all_vars)
    v = Vector{Float64}(undef, n)
    @inbounds for k in 1:n
        sv = JuMP.start_value(all_vars[k])
        v[k] = sv === nothing ? 0.0 : Float64(sv)
    end
    return v
end

# ============================================================
# Box bounds
# ============================================================

function _collect_bounds(all_vars::Vector{JuMP.VariableRef})
    n = length(all_vars)
    lb = fill(NaN, n)
    ub = fill(NaN, n)
    @inbounds for i in 1:n
        v = all_vars[i]
        if JuMP.has_lower_bound(v)
            lb[i] = JuMP.lower_bound(v)
        end
        if JuMP.has_upper_bound(v)
            ub[i] = JuMP.upper_bound(v)
        end
    end
    return lb, ub
end


@inline function _push_affine_row!(
    rows::Vector{Int},
    cols::Vector{Int},
    vals::Vector{Float64},
    row::Int,
    func::JuMP.AffExpr,
    var_to_col::Dict{JuMP.VariableRef, Int};
    scale::Float64 = 1.0,
)
    for (v, coef) in func.terms
        col = get(var_to_col, v, 0)
        col == 0 && continue
        push!(rows, row)
        push!(cols, col)
        push!(vals, scale * coef)
    end
    return nothing
end

@inline function _push_affine_row!(
    ::Vector{Int},
    ::Vector{Int},
    ::Vector{Float64},
    ::Int,
    func,
    ::Dict,
)
    return error(
        "ALADIN sensitivity extraction: expected affine equality constraint, " *
        "got $(typeof(func)). Canonical ALADIN DC/SOCWR only supports linear " *
        "local equalities.",
    )
end

# ============================================================
# Equality Jacobian (sparse, from JuMP affine constraints)
# ============================================================

function _extract_equality_jacobian(
    model::JuMP.Model,
    var_to_col::Dict{JuMP.VariableRef, Int},
    n_vars::Int,
)
    rows = Int[]
    cols = Int[]
    vals = Float64[]
    nrows = 0

    for (F, S) in JuMP.list_of_constraint_types(model)
        F <: JuMP.VariableRef && continue
        S === MOI.EqualTo{Float64} || continue

        for ci in JuMP.all_constraints(model, F, S)
            cobj = JuMP.constraint_object(ci)
            func = cobj.func
            nrows += 1
            _push_affine_row!(rows, cols, vals, nrows, func, var_to_col)
        end
    end

    return SparseArrays.sparse(rows, cols, vals, nrows, n_vars)
end

# ============================================================
# Affine inequality rows (active-set Jacobian assembly)
# ============================================================

function _collect_affine_inequality_rows!(
    rows::Vector{Int},
    cols::Vector{Int},
    vals::Vector{Float64},
    row_offset::Int,
    model::JuMP.Model,
    var_to_col::Dict{JuMP.VariableRef, Int},
    has_duals::Bool;
    dual_tol::Float64,
)::Int
    for (F, S) in JuMP.list_of_constraint_types(model)
        F <: JuMP.VariableRef && continue
        S === MOI.EqualTo{Float64} && continue
        F <: JuMP.AffExpr || continue

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
            _push_affine_row!(rows, cols, vals, row_offset, cobj.func, var_to_col)
        end
    end
    return row_offset
end

# ============================================================
# LICQ enforcement (rank-revealing row reduction)
# ============================================================

function _delete_dependent_rows(
    C::SparseArrays.SparseMatrixCSC{Float64, Int};
    tol::Float64 = 1e-10,
)::SparseArrays.SparseMatrixCSC{Float64, Int}
    n_rows = size(C, 1)
    n_rows <= 1 && return C

    A = Matrix(C')
    if size(A, 1) == 0 || size(A, 2) == 0
        return C
    end

    F = LinearAlgebra.qr(A, LinearAlgebra.ColumnNorm())
    R = F.R
    r_diag_len = min(size(R, 1), size(R, 2))
    pivot_scale = r_diag_len >= 1 ? max(abs(R[1, 1]), 1.0) : 1.0
    rnk = 0
    @inbounds for i in 1:r_diag_len
        if abs(R[i, i]) > tol * pivot_scale
            rnk += 1
        end
    end

    rnk == n_rows && return C
    keep = sort(F.p[1:rnk])
    return C[keep, :]
end
