# ======================================================================
# Distributed Algorithm Accessors
# ======================================================================
#
# The abstract type `DistributedAlgorithm` and the tag structs
# `ADMM{C,F}` / `ALADIN{C,F}` are declared in src/types.jl (canonical
# location for all top-level types). This file adds the accessor layer
# used by the rest of the distributed_algorithms/ subtree so no call
# site needs to remember the (V, M) parameter order.
#
# Canonical reading:
#
#   F <: OPFModel      — formulation (DC, AC, SOCWR, …)
#   C <: ADMMVariant   — coupling strategy (State, Power, Hybrid)
#
# ======================================================================

@inline formulation(::ADMM{C, F}) where {C, F} = F
@inline formulation(::ALADIN{C, F}) where {C, F} = F
@inline coupling(::ADMM{C, F}) where {C, F} = C
@inline coupling(::ALADIN{C, F}) where {C, F} = C

@inline algorithm_family(::ADMM) = :admm
@inline algorithm_family(::ALADIN) = :aladin
