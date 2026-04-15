# ============================================================
# DC WORKSPACES
# ------------------------------------------------------------

Base.@kwdef mutable struct DCStateWorkspace <: AbstractADMMWorkspace
    partitions::Partitions{PartitionDC}
    id::Vector{Int}

    vars::NamedTuple{(:theta,), Tuple{VarState}}

    bus_counts::Vector{Int}

    r_primal::Float64 = 0.0
    r_dual::Float64 = 0.0
    eps_primal::Float64 = 0.0
    eps_dual::Float64 = 0.0
    iter::Int = 0
    hist::ADMMHistory = ADMMHistory()
end

Base.@kwdef mutable struct DCPowerWorkspace <: AbstractADMMWorkspace
    partitions::Partitions{PartitionDC}
    id::Vector{Int}

    interzonal::Vector{Int}

    vars::NamedTuple{(:p_flow,), Tuple{VarState}}

    branch_counts::Vector{Int}

    r_primal::Float64 = 0.0
    r_dual::Float64 = 0.0
    eps_primal::Float64 = 0.0
    eps_dual::Float64 = 0.0
    iter::Int = 0
    hist::ADMMHistory = ADMMHistory()
end

Base.@kwdef mutable struct DCHybridWorkspace <: AbstractADMMWorkspace
    partitions::Partitions{PartitionDC}
    id::Vector{Int}

    interzonal::Vector{Int}

    vars::NamedTuple{(:theta, :p_flow), Tuple{VarState, VarState}}

    bus_counts::Vector{Int}
    branch_counts::Vector{Int}

    r_primal_bus::Float64 = 0.0
    r_dual_bus::Float64 = 0.0
    r_primal_branch::Float64 = 0.0
    r_dual_branch::Float64 = 0.0

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

_get_rhos(ws::DCStateWorkspace) = (rho_theta = ws.vars.theta.rho[],)
_get_rhos(ws::DCPowerWorkspace) = (rho_p = ws.vars.p_flow.rho[],)
_get_rhos(ws::DCHybridWorkspace) = (
    rho_theta = ws.vars.theta.rho[],
    rho_p = ws.vars.p_flow.rho[],
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
        theta = VarState(n_shared, n_zones, config.rho_theta),
        p_flow = VarState(n_inter, n_zones, config.rho_p),
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
