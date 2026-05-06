# ============================================================
# ALADIN-DC — Active-set detection & Jacobian assembly
# ============================================================
#
# At each iteration the coordinator QP needs the Jacobian of currently
# *active* constraints (LICQ-rank rows). For the DC formulation:
#   - equalities are constant (cached in `zone.eq_jac` at build time);
#   - inequalities at box bounds are detected via `detect_active_bounds!`;
#   - active inequality rows are appended to `eq_jac` and the stacked
#     matrix is run through `_delete_dependent_rows` for rank revealing.
#
# The QR with column-norm pivoting in `_delete_dependent_rows` is what
# enforces LICQ at the algorithm level — see header note in that
# function for the densification trade-off.

# ============================================================
# Active bounds at current iterate
# ============================================================

function detect_active_bounds!(zone::AladinZoneData; tol::Float64 = 1e-6)
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
# Active-constraint Jacobian (equalities + active inequalities)
# ============================================================

function build_active_jacobian(
    zone::AladinZoneData;
    dual_tol::Float64 = 1e-6,
    dep_tol::Float64 = 1e-10,
)::SparseArrays.SparseMatrixCSC{Float64, Int}
    model = zone.model
    var_to_col = zone.var_to_col
    n_vars = zone.n_vars

    # Start with the precomputed equality rows.
    eq_jac = zone.eq_jac
    eq_nrows = size(eq_jac, 1)

    # Collect active inequality rows by walking constraint types.
    ineq_rows = Int[]
    ineq_cols = Int[]
    ineq_vals = Vector{Float64}()
    row_offset = 0

    has_duals = JuMP.has_duals(model)

    for (F, S) in JuMP.list_of_constraint_types(model)
        F <: JuMP.VariableRef && continue
        S === MOI.EqualTo{Float64} && continue  # already in eq_jac

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

# ============================================================
# Rank-revealing row reduction (LICQ enforcement)
# ============================================================

function _delete_dependent_rows(
    C::SparseArrays.SparseMatrixCSC{Float64, Int};
    tol::Float64 = 1e-10,
)::SparseArrays.SparseMatrixCSC{Float64, Int}
    n_rows = size(C, 1)
    n_rows <= 1 && return C

    # Dense transpose: columns of A correspond to rows of C.
    # Sparse QR with column pivoting is not available in Julia stdlib
    # (`LinearAlgebra.qr(::SparseMatrixCSC, ::ColumnNorm)` is unsupported);
    # for robust LICQ rank detection we therefore densify here. This is
    # O(n_rows · n_vars) memory and O(n_rows²) time per call. Acceptable
    # for case ≤ 200 buses; for larger cases this is a known performance
    # ceiling and could be replaced by a sparse rank-revealing factorisation
    # (SuiteSparse SPQR or a streaming QR) as future work.
    A = Matrix(C')
    if size(A, 1) == 0 || size(A, 2) == 0
        return C
    end

    F = LinearAlgebra.qr(A, LinearAlgebra.ColumnNorm())

    # Count diagonal entries of R above tolerance → rank.
    r_diag_len = min(size(F.R, 1), size(F.R, 2))
    rnk = 0
    @inbounds for i in 1:r_diag_len
        if abs(F.R[i, i]) > tol
            rnk += 1
        end
    end

    rnk == n_rows && return C
    keep = sort(F.p[1:rnk])
    return C[keep, :]
end

# ============================================================
# Helper: push affine row of an inequality into COO buffers
# (separate from the equality version in `zone_data.jl` because
# inequalities use a runtime row offset rather than nrows counter)
# ============================================================

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
