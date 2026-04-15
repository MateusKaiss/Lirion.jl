# ======================================================================
# ADMM core support utilities
# ======================================================================

@inline _wrap_angle_diff(Δ::Float64) = atan(sin(Δ), cos(Δ))

@inline stop_if_converged(ws::AbstractADMMWorkspace) =
    (ws.r_primal < ws.eps_primal) && (ws.r_dual < ws.eps_dual)
