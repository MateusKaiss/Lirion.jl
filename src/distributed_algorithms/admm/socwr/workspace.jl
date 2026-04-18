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

    # SplitBus (Method 3) bookkeeping: sign-flipped p_inj / q_inj consensus
    # per split bus. Empty for TieLine.
    bus_pairs::Vector{NTuple{3, Int}} = NTuple{3, Int}[]
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}} =
        Vector{Dict{Int, JuMP.VariableRef}}()
    q_inj_vars::Vector{Dict{Int, JuMP.VariableRef}} =
        Vector{Dict{Int, JuMP.VariableRef}}()

    vars::NamedTuple{
        (:v2, :p_flow, :q_flow, :r, :t, :p_inj, :q_inj),
        Tuple{
            VarState, VarState, VarState, VarState, VarState, VarState, VarState,
        },
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
    r_primal_p_inj::Float64 = 0.0
    r_dual_p_inj::Float64 = 0.0
    r_primal_q_inj::Float64 = 0.0
    r_dual_q_inj::Float64 = 0.0

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
    rho_p_inj = ws.vars.p_inj.rho[],
    rho_q_inj = ws.vars.q_inj.rho[],
)

# ============================================================
# BUILDERS — STATE
# ============================================================

function build_workspace(sys, config, optimizer_factory, ::ADMM{State, SOCWR})
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

function build_workspace(sys, config, optimizer_factory, ::ADMM{Power, SOCWR})
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

function build_workspace(sys, config, optimizer_factory, ::ADMM{Hybrid, SOCWR})
    partitions = build_partitions(PartitionSOCWR, sys)
    id = partitions.zones
    inter = sys.interzonal_branches
    bpairs = sys.bus_pairs
    n_zones = length(partitions)
    n_shared = length(sys.full_shared_bus_ids)
    n_inter = length(inter)
    n_bus_pairs = length(bpairs)

    # Per-zone p_inj / q_inj containers: populated by init_hybrid_partition!
    # only when bus_pairs is non-empty (SplitBus sectioning).
    p_inj_vars = [Dict{Int, JuMP.VariableRef}() for _ in 1:n_zones]
    q_inj_vars = [Dict{Int, JuMP.VariableRef}() for _ in 1:n_zones]

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
        bpairs,
        p_inj_vars,
        q_inj_vars,
    )

    vars = (
        v2 = VarState(n_shared, n_zones, config.rho_v2; init_value = 1.0),
        p_flow = VarState(n_inter, n_zones, config.rho_p),
        q_flow = VarState(n_inter, n_zones, config.rho_q),
        r = VarState(n_inter, n_zones, config.rho_r),
        t = VarState(n_inter, n_zones, config.rho_t),
        # Virtual injections at split buses (SplitBus Method 3).
        p_inj = VarState(n_bus_pairs, n_zones, config.rho_p),
        q_inj = VarState(n_bus_pairs, n_zones, config.rho_q),
    )

    return SOCWRHybridWorkspace(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        bus_pairs = bpairs,
        p_inj_vars = p_inj_vars,
        q_inj_vars = q_inj_vars,
        vars = vars,
        bus_counts = zeros(Int, n_shared),
        branch_counts = zeros(Int, n_inter),
    )
end
