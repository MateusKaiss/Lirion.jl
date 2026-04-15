# ============================================================
# DC WORKSPACES
# ------------------------------------------------------------
# Layout conventions:
#   * Shared-bus quantities (theta) live in `Vector{Float64}` /
#     `Matrix{Float64}` of shape `(n_shared, n_zones)`, column-major,
#     indexed by `part.shared_bus_consensus_idx[i]` × `part.position`.
#   * Interzonal-branch quantities live in `Vector{Float64}` of length
#     `n_inter` and `Matrix{Float64}` of shape `(n_inter, n_zones)`,
#     indexed by position in `interzonal` × `part.position`.
#   * Scratch integer counters are preallocated once.
#
# Phase 4b: per-consensus-variable data (consensus vector, duals
# matrix, rho refs) lives inside a typed `NamedTuple` of `VarState`
# objects, keyed by the canonical variable symbols declared in
# `consensus_vars(F, C)`. Legacy flat names (`ws.theta_consensus`,
# `ws.duals_theta`, `ws.rho_theta`, …) keep working via the
# `getproperty` forwarding layer in `admm/core/components.jl`.
# ============================================================

Base.@kwdef mutable struct DCStateWorkspace <: AbstractADMMWorkspace
    partitions::Partitions{PartitionDC}
    id::Vector{Int}

    vars::NamedTuple{(:theta,), Tuple{VarState}}

    bus_counts::Vector{Int}

    residuals::ResidualMetrics = ResidualMetrics()
    state::IterationState = IterationState()
end

Base.@kwdef mutable struct DCPowerWorkspace <: AbstractADMMWorkspace
    partitions::Partitions{PartitionDC}
    id::Vector{Int}

    interzonal::Vector{Int}

    vars::NamedTuple{(:p_flow,), Tuple{VarState}}

    branch_counts::Vector{Int}

    residuals::ResidualMetrics = ResidualMetrics()
    state::IterationState = IterationState()
end

Base.@kwdef mutable struct DCHybridWorkspace <: AbstractADMMWorkspace
    partitions::Partitions{PartitionDC}
    id::Vector{Int}

    interzonal::Vector{Int}

    vars::NamedTuple{(:theta, :p_flow), Tuple{VarState, VarState}}

    bus_counts::Vector{Int}
    branch_counts::Vector{Int}

    # Per-component residuals (cached for rho_update! without recomputation)
    r_primal_bus::Float64 = 0.0
    r_dual_bus::Float64 = 0.0
    r_primal_branch::Float64 = 0.0
    r_dual_branch::Float64 = 0.0

    residuals::ResidualMetrics = ResidualMetrics()
    state::IterationState = IterationState()
end

# ============================================================
# RHO HELPERS
# ============================================================

_get_rhos(ws::DCStateWorkspace) = (rho_theta = ws.rho_theta[],)
_get_rhos(ws::DCPowerWorkspace) = (rho_p = ws.rho_p[],)
_get_rhos(ws::DCHybridWorkspace) = (
    rho_theta = ws.rho_theta[],
    rho_p = ws.rho_p[],
)

# ============================================================
# BUILDERS
# ============================================================

function build_admm_workspace(sys, config, optimizer_factory, ::ADMM{State, DC})
    partitions = build_partitions(PartitionDC, sys)
    id = partitions.zones
    n_zones = length(partitions)
    n_shared = length(sys.full_shared_bus_ids)

    init_state_partition!(partitions, sys, optimizer_factory, config.rho_theta)

    vars = (
        theta = VarState(n_shared, n_zones, config.rho_theta),
    )

    return DCStateWorkspace(;
        partitions = partitions,
        id = id,
        vars = vars,
        bus_counts = zeros(Int, n_shared),
    )
end

function build_admm_workspace(sys, config, optimizer_factory, ::ADMM{Power, DC})
    partitions = build_partitions(PartitionDC, sys)
    id = partitions.zones
    inter = sys.interzonal_branches
    n_zones = length(partitions)
    n_inter = length(inter)

    init_power_partition!(partitions, sys, optimizer_factory, config.rho_p, inter)

    vars = (
        p_flow = VarState(n_inter, n_zones, config.rho_p),
    )

    return DCPowerWorkspace(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        vars = vars,
        branch_counts = zeros(Int, n_inter),
    )
end

function build_admm_workspace(sys, config, optimizer_factory, ::ADMM{Hybrid, DC})
    partitions = build_partitions(PartitionDC, sys)
    id = partitions.zones
    inter = sys.interzonal_branches
    n_zones = length(partitions)
    n_shared = length(sys.full_shared_bus_ids)
    n_inter = length(inter)

    init_hybrid_partition!(
        partitions,
        sys,
        optimizer_factory,
        config.rho_theta,
        config.rho_p,
        inter,
    )

    vars = (
        theta  = VarState(n_shared, n_zones, config.rho_theta),
        p_flow = VarState(n_inter,  n_zones, config.rho_p),
    )

    return DCHybridWorkspace(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        vars = vars,
        bus_counts = zeros(Int, n_shared),
        branch_counts = zeros(Int, n_inter),
    )
end
