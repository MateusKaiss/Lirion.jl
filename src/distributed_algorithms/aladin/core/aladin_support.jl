# =======================================================================
# ALADIN — Core support utilities
# =======================================================================
# Stopping criterion follows Houska-Frasch-Diehl 2016, Theorem 4.1:
#
#   r_primal   = ‖∑_i A_i y_i‖₂    → 0    (consensus feasibility)
#   r_progress = ‖y_i - x_iᵏ‖₂     → 0    (Phase-1 optimality / progress)
#
# Both components must fall below their Boyd-style tolerances.
# =======================================================================

abstract type AbstractALADINWorkspace <: AbstractDistributedWorkspace end

@inline stop_if_converged(ws::AbstractALADINWorkspace) =
    (ws.r_primal < ws.eps_primal) && (ws.r_progress < ws.eps_progress)

@inline function _aladin_totals(partitions)
    parts = partitions.parts
    total = 0.0
    @inbounds @simd for k in eachindex(parts)
        total += parts[k].generation_cost
    end
    return total
end
