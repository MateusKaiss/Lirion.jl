# ============================================================
# AC WORKSPACES
# ============================================================

# ---------------- STATE ----------------

Base.@kwdef mutable struct ACStateWorkspace <: AbstractADMMWorkspace
    partitions::Partitions{PartitionAC}
    id::Vector{Int}

    consensus_bus_idx_map::Dict{Int, Int}

    theta_consensus::Vector{Float64}
    v_mag_consensus::Vector{Float64}

    theta_consensus_prev::Vector{Float64}
    v_mag_consensus_prev::Vector{Float64}

    duals_theta::Dict{Int, Vector{Float64}}
    duals_v_mag::Dict{Int, Vector{Float64}}

    rho_theta::Base.RefValue{Float64}
    rho_v::Base.RefValue{Float64}

    prev_rho_theta::Base.RefValue{Float64}
    prev_rho_v::Base.RefValue{Float64}

    r_primal::Float64 = 0.0
    r_dual::Float64 = 0.0
    eps_primal::Float64 = 0.0
    eps_dual::Float64 = 0.0

    # Per-component residuals (cached for rho_update! without recomputation)
    r_primal_theta::Float64 = 0.0
    r_primal_v_mag::Float64 = 0.0
    r_dual_theta::Float64 = 0.0
    r_dual_v_mag::Float64 = 0.0

    hist::ADMMHistory = ADMMHistory()
    iter::Int = 0
end

# ---------------- POWER ----------------

Base.@kwdef mutable struct ACPowerWorkspace <: AbstractADMMWorkspace
    partitions::Partitions{PartitionAC}
    id::Vector{Int}

    interzonal::Vector{Int}

    p_flow_consensus::Dict{Int, Float64}
    q_flow_consensus::Dict{Int, Float64}

    p_flow_consensus_prev::Dict{Int, Float64}
    q_flow_consensus_prev::Dict{Int, Float64}

    duals_p_flow::Dict{Int, Dict{Int, Float64}}
    duals_q_flow::Dict{Int, Dict{Int, Float64}}

    rho_p::Base.RefValue{Float64}
    rho_q::Base.RefValue{Float64}

    prev_rho_p::Base.RefValue{Float64}
    prev_rho_q::Base.RefValue{Float64}

    r_primal::Float64 = 0.0
    r_dual::Float64 = 0.0
    eps_primal::Float64 = 0.0
    eps_dual::Float64 = 0.0

    # Per-component residuals (cached for rho_update! without recomputation)
    r_primal_p::Float64 = 0.0
    r_primal_q::Float64 = 0.0
    r_dual_p::Float64 = 0.0
    r_dual_q::Float64 = 0.0

    hist::ADMMHistory = ADMMHistory()
    iter::Int = 0
end

# ---------------- HYBRID ----------------

Base.@kwdef mutable struct ACHybridWorkspace <: AbstractADMMWorkspace
    partitions::Partitions{PartitionAC}
    id::Vector{Int}

    consensus_bus_idx_map::Dict{Int, Int}
    interzonal::Vector{Int}

    theta_consensus::Vector{Float64}
    v_mag_consensus::Vector{Float64}

    theta_consensus_prev::Vector{Float64}
    v_mag_consensus_prev::Vector{Float64}

    p_flow_consensus::Dict{Int, Float64}
    q_flow_consensus::Dict{Int, Float64}

    p_flow_consensus_prev::Dict{Int, Float64}
    q_flow_consensus_prev::Dict{Int, Float64}

    duals_theta::Dict{Int, Vector{Float64}}
    duals_v_mag::Dict{Int, Vector{Float64}}

    duals_p_flow::Dict{Int, Dict{Int, Float64}}
    duals_q_flow::Dict{Int, Dict{Int, Float64}}

    rho_theta::Base.RefValue{Float64}
    rho_v::Base.RefValue{Float64}
    rho_p::Base.RefValue{Float64}
    rho_q::Base.RefValue{Float64}

    prev_rho_theta::Base.RefValue{Float64}
    prev_rho_v::Base.RefValue{Float64}
    prev_rho_p::Base.RefValue{Float64}
    prev_rho_q::Base.RefValue{Float64}

    r_primal::Float64 = 0.0
    r_dual::Float64 = 0.0
    eps_primal::Float64 = 0.0
    eps_dual::Float64 = 0.0

    # Per-component residuals (cached for rho_update! without recomputation)
    r_primal_theta::Float64 = 0.0
    r_primal_v_mag::Float64 = 0.0
    r_dual_theta::Float64 = 0.0
    r_dual_v_mag::Float64 = 0.0
    r_primal_p::Float64 = 0.0
    r_primal_q::Float64 = 0.0
    r_dual_p::Float64 = 0.0
    r_dual_q::Float64 = 0.0

    hist::ADMMHistory = ADMMHistory()
    iter::Int = 0
end

# ============================================================
# RHO HELPERS
# ============================================================

_get_rhos(ws::ACStateWorkspace) = Dict(
    :rho_theta => ws.rho_theta[],
    :rho_v => ws.rho_v[],
)

_get_rhos(ws::ACPowerWorkspace) = Dict(
    :rho_p => ws.rho_p[],
    :rho_q => ws.rho_q[],
)

_get_rhos(ws::ACHybridWorkspace) = Dict(
    :rho_theta => ws.rho_theta[],
    :rho_v => ws.rho_v[],
    :rho_p => ws.rho_p[],
    :rho_q => ws.rho_q[],
)

# ============================================================
# BUILDERS
# ============================================================

# ---------------- STATE ----------------

function build_admm_workspace(sys, config, optimizer_factory, ::ADMM{State, AC})
    partitions = build_partitions(PartitionAC, sys)
    id = partitions.zones

    consensus_bus_idx_map = build_consensus_bus_idx_map(sys.full_shared_bus_ids)
    n_shared = length(sys.full_shared_bus_ids)

    init_state_partition!(
        partitions,
        sys,
        optimizer_factory,
        config.rho_theta,
        config.rho_v,
    )

    rho_theta = Ref(config.rho_theta)
    rho_v = Ref(config.rho_v)

    return ACStateWorkspace(;
        partitions = partitions,
        id = id,
        consensus_bus_idx_map = consensus_bus_idx_map,
        theta_consensus = zeros(n_shared),
        v_mag_consensus = ones(n_shared),
        theta_consensus_prev = zeros(n_shared),
        v_mag_consensus_prev = ones(n_shared),
        duals_theta = Dict(z => zeros(n_shared) for z in id),
        duals_v_mag = Dict(z => zeros(n_shared) for z in id),
        rho_theta = rho_theta,
        rho_v = rho_v,
        prev_rho_theta = Ref(rho_theta[]),
        prev_rho_v = Ref(rho_v[]),
    )
end

# ---------------- POWER ----------------

function build_admm_workspace(sys, config, optimizer_factory, ::ADMM{Power, AC})
    partitions = build_partitions(PartitionAC, sys)
    id = partitions.zones
    inter = sys.interzonal_branches

    init_power_partition!(
        partitions,
        sys,
        optimizer_factory,
        config.rho_p,
        config.rho_q,
        inter,
    )

    rho_p = Ref(config.rho_p)
    rho_q = Ref(config.rho_q)

    return ACPowerWorkspace(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        p_flow_consensus = Dict(e => 0.0 for e in inter),
        q_flow_consensus = Dict(e => 0.0 for e in inter),
        p_flow_consensus_prev = Dict(e => 0.0 for e in inter),
        q_flow_consensus_prev = Dict(e => 0.0 for e in inter),
        duals_p_flow = Dict(z => Dict{Int, Float64}() for z in id),
        duals_q_flow = Dict(z => Dict{Int, Float64}() for z in id),
        rho_p = rho_p,
        rho_q = rho_q,
        prev_rho_p = Ref(rho_p[]),
        prev_rho_q = Ref(rho_q[]),
    )
end

# ---------------- HYBRID ----------------

function build_admm_workspace(sys, config, optimizer_factory, ::ADMM{Hybrid, AC})
    partitions = build_partitions(PartitionAC, sys)
    id = partitions.zones
    inter = sys.interzonal_branches

    consensus_bus_idx_map = build_consensus_bus_idx_map(sys.full_shared_bus_ids)
    n_shared = length(sys.full_shared_bus_ids)

    init_hybrid_partition!(
        partitions,
        sys,
        optimizer_factory,
        config.rho_theta,
        config.rho_v,
        config.rho_p,
        config.rho_q,
        inter,
    )

    rho_theta = Ref(config.rho_theta)
    rho_v = Ref(config.rho_v)
    rho_p = Ref(config.rho_p)
    rho_q = Ref(config.rho_q)

    return ACHybridWorkspace(;
        partitions = partitions,
        id = id,
        consensus_bus_idx_map = consensus_bus_idx_map,
        interzonal = inter,
        theta_consensus = zeros(n_shared),
        v_mag_consensus = ones(n_shared),
        theta_consensus_prev = zeros(n_shared),
        v_mag_consensus_prev = ones(n_shared),
        p_flow_consensus = Dict(e => 0.0 for e in inter),
        q_flow_consensus = Dict(e => 0.0 for e in inter),
        p_flow_consensus_prev = Dict(e => 0.0 for e in inter),
        q_flow_consensus_prev = Dict(e => 0.0 for e in inter),
        duals_theta = Dict(z => zeros(n_shared) for z in id),
        duals_v_mag = Dict(z => zeros(n_shared) for z in id),
        duals_p_flow = Dict(z => Dict{Int, Float64}() for z in id),
        duals_q_flow = Dict(z => Dict{Int, Float64}() for z in id),
        rho_theta = rho_theta,
        rho_v = rho_v,
        rho_p = rho_p,
        rho_q = rho_q,
        prev_rho_theta = Ref(rho_theta[]),
        prev_rho_v = Ref(rho_v[]),
        prev_rho_p = Ref(rho_p[]),
        prev_rho_q = Ref(rho_q[]),
    )
end
