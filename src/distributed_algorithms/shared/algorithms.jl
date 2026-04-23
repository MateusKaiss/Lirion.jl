# ============================================================
# Type accessors
# ============================================================

@inline formulation(::ADMM{C, F}) where {C, F} = F
@inline formulation(::ALADIN{C, F}) where {C, F} = F
@inline coupling(::ADMM{C, F}) where {C, F} = C
@inline coupling(::ALADIN{C, F}) where {C, F} = C

@inline algorithm_family(::ADMM) = :admm
@inline algorithm_family(::ALADIN) = :aladin

# ============================================================
# Algorithm tags
# ============================================================

algorithm_tag(::ADMMAlgorithm) = :admm
algorithm_tag(::ALADINAlgorithm) = :aladin

# ============================================================
# Time utilities
# ============================================================

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
