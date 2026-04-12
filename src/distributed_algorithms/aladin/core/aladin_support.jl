# =======================================================================
# ALADIN — Core support utilities
# =======================================================================

# ---------------------------------------------------------------
# Abstract workspace type (parent of all per-model ALADIN workspaces)
# ---------------------------------------------------------------

abstract type AbstractALADINWorkspace end

# ---------------------------------------------------------------
# Generic convergence / helper predicates
# ---------------------------------------------------------------

@inline function stop_if_converged_aladin(ws::AbstractALADINWorkspace)
    return ws.r_primal < ws.eps_primal && ws.r_dual < ws.eps_dual
end

@inline function _aladin_totals(partitions)
    parts = partitions.parts
    total = 0.0
    @inbounds @simd for k in eachindex(parts)
        total += parts[k].generation_cost
    end
    return total
end

# ---------------------------------------------------------------
# Logging tags
# ---------------------------------------------------------------

_aladin_print_tag(ws::AbstractALADINWorkspace) = "[ALADIN]"
