# ======================================================================
# Distributed Algorithm Accessors
# ======================================================================
#
# The abstract type `DistributedAlgorithm` and the tag structs
# `ADMM{C,F}` / `ALADIN{C,F}` are declared in src/types.jl (the
# canonical location for all top-level types). This file adds the
# accessor + utility layer used by the rest of the
# distributed_algorithms/ subtree so no call site needs to remember the
# (V, M) parameter order or re-implement trivial timing helpers.
#
# Canonical reading:
#
#   F <: OPFModel      — formulation (DC, AC, SOCWR, …)
#   C <: ADMMVariant   — coupling strategy (State, Power, Hybrid)
#
# ======================================================================

# ----------------------------------------------------------------------
# Type-level accessors — extract (F, C) from any distributed tag
# ----------------------------------------------------------------------

@inline formulation(::ADMM{C, F}) where {C, F} = F
@inline formulation(::ALADIN{C, F}) where {C, F} = F
@inline coupling(::ADMM{C, F}) where {C, F} = C
@inline coupling(::ALADIN{C, F}) where {C, F} = C

@inline algorithm_family(::ADMM) = :admm
@inline algorithm_family(::ALADIN) = :aladin

# ----------------------------------------------------------------------
# Short diagnostic tags used by logging and history modules
# ----------------------------------------------------------------------
# `algorithm_tag(::ALADINAlgorithm) = :aladin` is added in
# aladin/run_aladin.jl so this file doesn't need to know about the
# ALADIN subtree types.

algorithm_tag(::ADMMAlgorithm) = :admm

# ----------------------------------------------------------------------
# Timing utilities shared by both run_admm and run_aladin
# ----------------------------------------------------------------------

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
