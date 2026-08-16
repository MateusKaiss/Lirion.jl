# ============================================================
# ALADIN-AC — Active-constraint Jacobian
# ============================================================

@inline build_active_jacobian(
    zone::AladinZoneDataAC;
    kwargs...,
)::SparseArrays.SparseMatrixCSC{Float64, Int} = zone.active_jac
