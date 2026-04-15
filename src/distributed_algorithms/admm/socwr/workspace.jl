# ============================================================
# SOCWR WORKSPACES
# ============================================================

# ============================================================
# STATE (V² + R/T)
# ============================================================

Base.@kwdef mutable struct SOCWRStateWorkspace <: AbstractADMMWorkspace
    partitions::Partitions{PartitionSOCWR}
    id::Vector{Int}

    interzonal::Vector{Int}

    vars::NamedTuple{(:v2, :r, :t), Tuple{VarState, VarState, VarState}}

    bus_counts::Vector{Int}
    branch_counts::Vector{Int}

    # Per-component residuals (cached for rho_update! without recomputation)
    r_primal_bus::Float64 = 0.0
    r_dual_bus::Float64 = 0.0
    r_primal_r::Float64 = 0.0
    r_primal_t::Float64 = 0.0
    r_dual_r::Float64 = 0.0
    r_dual_t::Float64 = 0.0

    r_primal::Float64 = 0.0
    r_dual::Float64 = 0.0
    eps_primal::Float64 = 0.0
    eps_dual::Float64 = 0.0
    iter::Int = 0
    hist::ADMMHistory = ADMMHistory()
end

# ============================================================
# POWER (P/Q)
# ============================================================

Base.@kwdef mutable struct SOCWRPowerWorkspace <: AbstractADMMWorkspace
    partitions::Partitions{PartitionSOCWR}
    id::Vector{Int}

    interzonal::Vector{Int}

    vars::NamedTuple{(:p_flow, :q_flow), Tuple{VarState, VarState}}

    branch_counts::Vector{Int}

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

# ============================================================
# HYBRID (V² + P/Q + R/T)
# ============================================================

Base.@kwdef mutable struct SOCWRHybridWorkspace <: AbstractADMMWorkspace
    partitions::Partitions{PartitionSOCWR}
    id::Vector{Int}

    interzonal::Vector{Int}

    vars::NamedTuple{
        (:v2, :p_flow, :q_flow, :r, :t),
        Tuple{VarState, VarState, VarState, VarState, VarState},
    }

    bus_counts::Vector{Int}
    branch_counts::Vector{Int}

    # Per-component residuals (cached for rho_update! without recomputation)
    r_primal_bus::Float64 = 0.0
    r_dual_bus::Float64 = 0.0
    r_primal_p::Float64 = 0.0
    r_primal_q::Float64 = 0.0
    r_dual_p::Float64 = 0.0
    r_dual_q::Float64 = 0.0
    r_primal_r::Float64 = 0.0
    r_primal_t::Float64 = 0.0
    r_dual_r::Float64 = 0.0
    r_dual_t::Float64 = 0.0

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

_get_rhos(ws::SOCWRStateWorkspace) = (
    rho_v2 = ws.vars.v2.rho[],
    rho_r = ws.vars.r.rho[],
    rho_t = ws.vars.t.rho[],
)

_get_rhos(ws::SOCWRPowerWorkspace) = (
    rho_p = ws.vars.p_flow.rho[],
    rho_q = ws.vars.q_flow.rho[],
)

_get_rhos(ws::SOCWRHybridWorkspace) = (
    rho_v2 = ws.vars.v2.rho[],
    rho_p = ws.vars.p_flow.rho[],
    rho_q = ws.vars.q_flow.rho[],
    rho_r = ws.vars.r.rho[],
    rho_t = ws.vars.t.rho[],
)

# ============================================================
# BUILDERS — STATE
# ============================================================

function build_admm_workspace(sys, config, optimizer_factory, ::ADMM{State, SOCWR})
    partitions = build_partitions(PartitionSOCWR, sys)
    id = partitions.zones
    inter = sys.interzonal_branches
    n_zones = length(partitions)
    n_shared = length(sys.full_shared_bus_ids)
    n_inter = length(inter)

    init_state_partition!(
        partitions,
        sys,
        optimizer_factory,
        config.rho_v2,
        config.rho_r,
        config.rho_t,
        inter,
    )

    vars = (
        v2 = VarState(n_shared, n_zones, config.rho_v2; init_value = 1.0),
        r = VarState(n_inter, n_zones, config.rho_r),
        t = VarState(n_inter, n_zones, config.rho_t),
    )

    return SOCWRStateWorkspace(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        vars = vars,
        bus_counts = zeros(Int, n_shared),
        branch_counts = zeros(Int, n_inter),
    )
end

# ============================================================
# BUILDERS — POWER
# ============================================================

function build_admm_workspace(sys, config, optimizer_factory, ::ADMM{Power, SOCWR})
    partitions = build_partitions(PartitionSOCWR, sys)
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

    return SOCWRPowerWorkspace(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        vars = vars,
        branch_counts = zeros(Int, n_inter),
    )
end

# ============================================================
# BUILDERS — HYBRID
# ============================================================

function build_admm_workspace(sys, config, optimizer_factory, ::ADMM{Hybrid, SOCWR})
    partitions = build_partitions(PartitionSOCWR, sys)
    id = partitions.zones
    inter = sys.interzonal_branches
    n_zones = length(partitions)
    n_shared = length(sys.full_shared_bus_ids)
    n_inter = length(inter)

    init_hybrid_partition!(
        partitions,
        sys,
        optimizer_factory,
        config.rho_v2,
        config.rho_p,
        config.rho_q,
        config.rho_r,
        config.rho_t,
        inter,
    )

    vars = (
        v2 = VarState(n_shared, n_zones, config.rho_v2; init_value = 1.0),
        p_flow = VarState(n_inter, n_zones, config.rho_p),
        q_flow = VarState(n_inter, n_zones, config.rho_q),
        r = VarState(n_inter, n_zones, config.rho_r),
        t = VarState(n_inter, n_zones, config.rho_t),
    )

    return SOCWRHybridWorkspace(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        vars = vars,
        bus_counts = zeros(Int, n_shared),
        branch_counts = zeros(Int, n_inter),
    )
end
