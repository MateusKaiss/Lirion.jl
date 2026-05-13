# ============================================================
# ALADIN-DC — Active-constraint Jacobian
# ============================================================

function build_active_jacobian(
    zone::AladinZoneData;
    dual_tol::Float64 = 1e-6,
    dep_tol::Float64 = 1e-10,
)::SparseArrays.SparseMatrixCSC{Float64, Int}
    ineq_rows = Int[]
    ineq_cols = Int[]
    ineq_vals = Float64[]

    row_offset = _collect_affine_inequality_rows!(
        ineq_rows, ineq_cols, ineq_vals, 0,
        zone.model, zone.var_to_col, JuMP.has_duals(zone.model);
        dual_tol = dual_tol,
    )

    C_stacked = if row_offset == 0
        zone.eq_jac
    else
        C_ineq = SparseArrays.sparse(
            ineq_rows, ineq_cols, ineq_vals, row_offset, zone.n_vars,
        )
        vcat(zone.eq_jac, C_ineq)
    end

    return _delete_dependent_rows(C_stacked; tol = dep_tol)
end
