# ============================================================
# AC WORKSPACES
# ============================================================
# ---------------- STATE ----------------

Base.@kwdef mutable struct ACStateWorkspace <: AbstractADMMWorkspace
    partitions::Partitions{PartitionAC}
    id::Vector{Int}

    vars::NamedTuple{(:theta, :v_mag), Tuple{VarState, VarState}}

    # Scratch vectors for circular mean (reused per iter)
    sin_sum::Vector{Float64}
    cos_sum::Vector{Float64}
    bus_counts::Vector{Int}

    # Per-component residuals (cached for rho_update! without recomputation)
    r_primal_theta::Float64 = 0.0
    r_primal_v_mag::Float64 = 0.0
    r_dual_theta::Float64 = 0.0
    r_dual_v_mag::Float64 = 0.0

    r_primal::Float64 = 0.0
    r_dual::Float64 = 0.0
    eps_primal::Float64 = 0.0
    eps_dual::Float64 = 0.0
    iter::Int = 0
    hist::ADMMHistory = ADMMHistory()
end

# ---------------- POWER ----------------

Base.@kwdef mutable struct ACPowerWorkspace <: AbstractADMMWorkspace
    partitions::Partitions{PartitionAC}
    id::Vector{Int}

    interzonal::Vector{Int}

    vars::NamedTuple{(:p_flow, :q_flow), Tuple{VarState, VarState}}

    # Scratch counters for NaN-safe AC branch consensus (per-iter).
    p_branch_counts::Vector{Int}
    q_branch_counts::Vector{Int}

    # Per-component residuals (cached for rho_update! without recomputation)
    r_primal_p::Float64 = 0.0
    r_primal_q::Float64 = 0.0
    r_dual_p::Float64 = 0.0
    r_dual_q::Float64 = 0.0

    r_primal::Float64 = 0.0
    r_dual::Float64 = 0.0
    eps_primal::Float64 = 0.0
    eps_dual::Float64 = 0.0
    iter::Int = 0
    hist::ADMMHistory = ADMMHistory()
end

# ---------------- HYBRID ----------------

Base.@kwdef mutable struct ACHybridWorkspace <: AbstractADMMWorkspace
    partitions::Partitions{PartitionAC}
    id::Vector{Int}

    interzonal::Vector{Int}

    vars::NamedTuple{
        (:theta, :v_mag, :p_flow, :q_flow),
        Tuple{VarState, VarState, VarState, VarState},
    }

    # Scratch buffers
    sin_sum::Vector{Float64}
    cos_sum::Vector{Float64}
    bus_counts::Vector{Int}
    p_branch_counts::Vector{Int}
    q_branch_counts::Vector{Int}

    # Per-component residuals (cached for rho_update! without recomputation)
    r_primal_theta::Float64 = 0.0
    r_primal_v_mag::Float64 = 0.0
    r_dual_theta::Float64 = 0.0
    r_dual_v_mag::Float64 = 0.0
    r_primal_p::Float64 = 0.0
    r_primal_q::Float64 = 0.0
    r_dual_p::Float64 = 0.0
    r_dual_q::Float64 = 0.0

    r_primal::Float64 = 0.0
    r_dual::Float64 = 0.0
    eps_primal::Float64 = 0.0
    eps_dual::Float64 = 0.0
    iter::Int = 0
    hist::ADMMHistory = ADMMHistory()
end

# ============================================================
# RHO HELPERS
# ============================================================

_get_rhos(ws::ACStateWorkspace) = (
    rho_theta = ws.vars.theta.rho[],
    rho_v = ws.vars.v_mag.rho[],
)

_get_rhos(ws::ACPowerWorkspace) = (
    rho_p = ws.vars.p_flow.rho[],
    rho_q = ws.vars.q_flow.rho[],
)

_get_rhos(ws::ACHybridWorkspace) = (
    rho_theta = ws.vars.theta.rho[],
    rho_v = ws.vars.v_mag.rho[],
    rho_p = ws.vars.p_flow.rho[],
    rho_q = ws.vars.q_flow.rho[],
)

# ============================================================
# BUILDERS
# ============================================================

# ---------------- STATE ----------------

function build_admm_workspace(sys, config, optimizer_factory, ::ADMM{State, AC})
    partitions = build_partitions(PartitionAC, sys)
    id = partitions.zones
    n_zones = length(partitions)
    n_shared = length(sys.full_shared_bus_ids)

    init_state_partition!(
        partitions,
        sys,
        optimizer_factory,
        config.rho_theta,
        config.rho_v,
    )

    vars = (
        theta = VarState(n_shared, n_zones, config.rho_theta),
        v_mag = VarState(n_shared, n_zones, config.rho_v; init_value = 1.0),
    )

    return ACStateWorkspace(;
        partitions = partitions,
        id = id,
        vars = vars,
        sin_sum = zeros(n_shared),
        cos_sum = zeros(n_shared),
        bus_counts = zeros(Int, n_shared),
    )
end

# ---------------- POWER ----------------

function build_admm_workspace(sys, config, optimizer_factory, ::ADMM{Power, AC})
    partitions = build_partitions(PartitionAC, sys)
    id = partitions.zones
    inter = sys.interzonal_branches
    n_zones = length(partitions)
    n_inter = length(inter)

    init_power_partition!(
        partitions,
        sys,
        optimizer_factory,
        config.rho_p,
        config.rho_q,
        inter,
    )

    vars = (
        p_flow = VarState(n_inter, n_zones, config.rho_p),
        q_flow = VarState(n_inter, n_zones, config.rho_q),
    )

    return ACPowerWorkspace(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        vars = vars,
        p_branch_counts = zeros(Int, n_inter),
        q_branch_counts = zeros(Int, n_inter),
    )
end

# ---------------- HYBRID ----------------

function build_admm_workspace(sys, config, optimizer_factory, ::ADMM{Hybrid, AC})
    partitions = build_partitions(PartitionAC, sys)
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
        config.rho_v,
        config.rho_p,
        config.rho_q,
        inter,
    )

    vars = (
        theta = VarState(n_shared, n_zones, config.rho_theta),
        v_mag = VarState(n_shared, n_zones, config.rho_v; init_value = 1.0),
        p_flow = VarState(n_inter, n_zones, config.rho_p),
        q_flow = VarState(n_inter, n_zones, config.rho_q),
    )

    return ACHybridWorkspace(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        vars = vars,
        sin_sum = zeros(n_shared),
        cos_sum = zeros(n_shared),
        bus_counts = zeros(Int, n_shared),
        p_branch_counts = zeros(Int, n_inter),
        q_branch_counts = zeros(Int, n_inter),
    )
end
