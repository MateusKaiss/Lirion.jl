# ============================================================
# DC WORKSPACES
# ============================================================

Base.@kwdef mutable struct DCStateWorkspace <: AbstractADMMWorkspace
    partitions::Partitions{PartitionDC}
    id::Vector{Int}

    consensus_bus_idx_map::Dict{Int, Int}

    theta_consensus::Vector{Float64}
    theta_consensus_prev::Vector{Float64}

    duals_theta::Dict{Int, Vector{Float64}}

    rho_theta::Base.RefValue{Float64}
    prev_rho_theta::Base.RefValue{Float64}

    r_primal::Float64 = 0.0
    r_dual::Float64 = 0.0
    eps_primal::Float64 = 0.0
    eps_dual::Float64 = 0.0

    hist::ADMMHistory = ADMMHistory()
    iter::Int = 0
end

Base.@kwdef mutable struct DCPowerWorkspace <: AbstractADMMWorkspace
    partitions::Partitions{PartitionDC}
    id::Vector{Int}

    interzonal::Vector{Int}

    p_flow_consensus::Dict{Int, Float64}
    p_flow_consensus_prev::Dict{Int, Float64}

    duals_p_flow::Dict{Int, Dict{Int, Float64}}

    rho_p::Base.RefValue{Float64}
    prev_rho_p::Base.RefValue{Float64}

    r_primal::Float64 = 0.0
    r_dual::Float64 = 0.0
    eps_primal::Float64 = 0.0
    eps_dual::Float64 = 0.0

    hist::ADMMHistory = ADMMHistory()
    iter::Int = 0
end

Base.@kwdef mutable struct DCHybridWorkspace <: AbstractADMMWorkspace
    partitions::Partitions{PartitionDC}
    id::Vector{Int}

    consensus_bus_idx_map::Dict{Int, Int}
    interzonal::Vector{Int}

    theta_consensus::Vector{Float64}
    theta_consensus_prev::Vector{Float64}

    p_flow_consensus::Dict{Int, Float64}
    p_flow_consensus_prev::Dict{Int, Float64}

    duals_theta::Dict{Int, Vector{Float64}}

    duals_p_flow::Dict{Int, Dict{Int, Float64}}

    rho_theta::Base.RefValue{Float64}
    rho_p::Base.RefValue{Float64}

    prev_rho_theta::Base.RefValue{Float64}
    prev_rho_p::Base.RefValue{Float64}

    r_primal::Float64 = 0.0
    r_dual::Float64 = 0.0
    eps_primal::Float64 = 0.0
    eps_dual::Float64 = 0.0

    # Per-component residuals (cached for rho_update! without recomputation)
    r_primal_bus::Float64 = 0.0
    r_dual_bus::Float64 = 0.0
    r_primal_branch::Float64 = 0.0
    r_dual_branch::Float64 = 0.0

    hist::ADMMHistory = ADMMHistory()
    iter::Int = 0
end

# ============================================================
# RHO HELPERS
# ============================================================

_get_rhos(ws::DCStateWorkspace) = Dict(:rho_theta => ws.rho_theta[])
_get_rhos(ws::DCPowerWorkspace) = Dict(:rho_p => ws.rho_p[])
_get_rhos(ws::DCHybridWorkspace) = Dict(
    :rho_theta => ws.rho_theta[],
    :rho_p => ws.rho_p[],
)

# ============================================================
# BUILDERS
# ============================================================

function build_admm_workspace(sys, config, optimizer_factory, ::ADMM{State, DC})
    partitions = build_partitions(PartitionDC, sys)
    id = partitions.zones

    consensus_bus_idx_map = build_consensus_bus_idx_map(sys.full_shared_bus_ids)
    n_shared = length(sys.full_shared_bus_ids)

    init_state_partition!(partitions, sys, optimizer_factory, config.rho_theta)

    rho_theta = Ref(config.rho_theta)

    return DCStateWorkspace(;
        partitions = partitions,
        id = id,
        consensus_bus_idx_map = consensus_bus_idx_map,
        theta_consensus = zeros(n_shared),
        theta_consensus_prev = zeros(n_shared),
        duals_theta = Dict(z => zeros(n_shared) for z in id),
        rho_theta = rho_theta,
        prev_rho_theta = Ref(rho_theta[]),
    )
end

function build_admm_workspace(sys, config, optimizer_factory, ::ADMM{Power, DC})
    partitions = build_partitions(PartitionDC, sys)
    id = partitions.zones
    inter = sys.interzonal_branches

    init_power_partition!(partitions, sys, optimizer_factory, config.rho_p, inter)

    rho_p = Ref(config.rho_p)

    return DCPowerWorkspace(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        p_flow_consensus = Dict(e => 0.0 for e in inter),
        p_flow_consensus_prev = Dict(e => 0.0 for e in inter),
        duals_p_flow = Dict(z => Dict{Int, Float64}() for z in id),
        rho_p = rho_p,
        prev_rho_p = Ref(rho_p[]),
    )
end

function build_admm_workspace(sys, config, optimizer_factory, ::ADMM{Hybrid, DC})
    partitions = build_partitions(PartitionDC, sys)
    id = partitions.zones
    inter = sys.interzonal_branches

    consensus_bus_idx_map = build_consensus_bus_idx_map(sys.full_shared_bus_ids)
    n_shared = length(sys.full_shared_bus_ids)

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
        consensus_bus_idx_map = consensus_bus_idx_map,
        interzonal = inter,
        theta_consensus = zeros(n_shared),
        theta_consensus_prev = zeros(n_shared),
        p_flow_consensus = Dict(e => 0.0 for e in inter),
        p_flow_consensus_prev = Dict(e => 0.0 for e in inter),
        duals_theta = Dict(z => zeros(n_shared) for z in id),
        duals_p_flow = Dict(z => Dict{Int, Float64}() for z in id),
        rho_theta = rho_theta,
        rho_p = rho_p,
        prev_rho_theta = Ref(rho_theta[]),
        prev_rho_p = Ref(rho_p[]),
    )
end
