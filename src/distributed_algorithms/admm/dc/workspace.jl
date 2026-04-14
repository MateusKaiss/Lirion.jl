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
# ============================================================

Base.@kwdef mutable struct DCStateWorkspace <: AbstractADMMWorkspace
    partitions::Partitions{PartitionDC}
    id::Vector{Int}

    theta_consensus::Vector{Float64}
    theta_consensus_prev::Vector{Float64}

    duals_theta::Matrix{Float64}

    rho_theta::Base.RefValue{Float64}
    prev_rho_theta::Base.RefValue{Float64}

    bus_counts::Vector{Int}

    residuals::ResidualMetrics = ResidualMetrics()
    state::IterationState = IterationState()
end

Base.@kwdef mutable struct DCPowerWorkspace <: AbstractADMMWorkspace
    partitions::Partitions{PartitionDC}
    id::Vector{Int}

    interzonal::Vector{Int}

    p_flow_consensus::Vector{Float64}
    p_flow_consensus_prev::Vector{Float64}

    duals_p_flow::Matrix{Float64}

    rho_p::Base.RefValue{Float64}
    prev_rho_p::Base.RefValue{Float64}

    branch_counts::Vector{Int}

    residuals::ResidualMetrics = ResidualMetrics()
    state::IterationState = IterationState()
end

Base.@kwdef mutable struct DCHybridWorkspace <: AbstractADMMWorkspace
    partitions::Partitions{PartitionDC}
    id::Vector{Int}

    interzonal::Vector{Int}

    theta_consensus::Vector{Float64}
    theta_consensus_prev::Vector{Float64}

    p_flow_consensus::Vector{Float64}
    p_flow_consensus_prev::Vector{Float64}

    duals_theta::Matrix{Float64}
    duals_p_flow::Matrix{Float64}

    rho_theta::Base.RefValue{Float64}
    rho_p::Base.RefValue{Float64}

    prev_rho_theta::Base.RefValue{Float64}
    prev_rho_p::Base.RefValue{Float64}

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

    rho_theta = Ref(config.rho_theta)

    return DCStateWorkspace(;
        partitions = partitions,
        id = id,
        theta_consensus = zeros(n_shared),
        theta_consensus_prev = zeros(n_shared),
        duals_theta = zeros(n_shared, n_zones),
        rho_theta = rho_theta,
        prev_rho_theta = Ref(rho_theta[]),
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

    rho_p = Ref(config.rho_p)

    return DCPowerWorkspace(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        p_flow_consensus = zeros(n_inter),
        p_flow_consensus_prev = zeros(n_inter),
        duals_p_flow = zeros(n_inter, n_zones),
        rho_p = rho_p,
        prev_rho_p = Ref(rho_p[]),
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

    rho_theta = Ref(config.rho_theta)
    rho_p = Ref(config.rho_p)

    return DCHybridWorkspace(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        theta_consensus = zeros(n_shared),
        theta_consensus_prev = zeros(n_shared),
        p_flow_consensus = zeros(n_inter),
        p_flow_consensus_prev = zeros(n_inter),
        duals_theta = zeros(n_shared, n_zones),
        duals_p_flow = zeros(n_inter, n_zones),
        rho_theta = rho_theta,
        rho_p = rho_p,
        prev_rho_theta = Ref(rho_theta[]),
        prev_rho_p = Ref(rho_p[]),
        bus_counts = zeros(Int, n_shared),
        branch_counts = zeros(Int, n_inter),
    )
end
