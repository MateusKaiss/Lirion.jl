# ============================================================
# ALADIN-DC — Per-iteration sensitivities
# ============================================================

# ------------------------------------------------------------
# Hessian refresh
# ------------------------------------------------------------

@inline refresh_lagrangian_hessian!(::AladinZoneData, ::AbstractALADINWorkspace) =
    nothing
