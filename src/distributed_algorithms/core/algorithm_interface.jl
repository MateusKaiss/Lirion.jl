# =======================================================================
# Abstract Distributed Algorithm Interface
# =======================================================================

# ---------------------------------------------------------------
# Tagging — short label used in logging and diagnostics.
# ---------------------------------------------------------------

algorithm_tag(::ADMMAlgorithm) = :admm
# algorithm_tag(::ALADINAlgorithm) = :aladin  # added in aladin/

# ---------------------------------------------------------------
# Coupling helpers — shared between ADMM and ALADIN
# ---------------------------------------------------------------

const VALID_COUPLINGS = (:state, :power, :hybrid)

@inline function _check_coupling(coupling::Symbol)
    coupling in VALID_COUPLINGS ||
        throw(ArgumentError("Invalid coupling=$coupling. Allowed: $(VALID_COUPLINGS)"))
    return coupling
end

# ---------------------------------------------------------------
# Generic timing / pretty-print utilities
# ---------------------------------------------------------------

@inline _elapsed_s(t0_ns::Integer) = (time_ns() - t0_ns) / 1e9

@inline function _pretty_time(sec::Real)
    if sec < 60
        return "$(round(sec; digits = 3)) s"
    elseif sec < 3600
        return "$(round(sec / 60; digits = 2)) min"
    else
        return "$(round(sec / 3600; digits = 2)) h"
    end
end
