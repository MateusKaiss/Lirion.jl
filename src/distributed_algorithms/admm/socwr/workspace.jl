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

    # BUS
    v2_consensus::Vector{Float64}
    v2_consensus_prev::Vector{Float64}

    # RT
    r_consensus::Vector{Float64}
    t_consensus::Vector{Float64}

    r_consensus_prev::Vector{Float64}
    t_consensus_prev::Vector{Float64}

    # DUALS
    duals_v2::Matrix{Float64}

    duals_r::Matrix{Float64}
    duals_t::Matrix{Float64}

    # RHOS
    rho_v2::Base.RefValue{Float64}
    rho_r::Base.RefValue{Float64}
    rho_t::Base.RefValue{Float64}

    prev_rho_v2::Base.RefValue{Float64}
    prev_rho_r::Base.RefValue{Float64}
    prev_rho_t::Base.RefValue{Float64}

    bus_counts::Vector{Int}
    branch_counts::Vector{Int}

    r_primal::Float64 = 0.0
    r_dual::Float64 = 0.0
    eps_primal::Float64 = 0.0
    eps_dual::Float64 = 0.0

    # Per-component residuals (cached for rho_update! without recomputation)
    r_primal_bus::Float64 = 0.0
    r_dual_bus::Float64 = 0.0
    r_primal_r::Float64 = 0.0
    r_primal_t::Float64 = 0.0
    r_dual_r::Float64 = 0.0
    r_dual_t::Float64 = 0.0

    hist::ADMMHistory = ADMMHistory()
    iter::Int = 0
end

# ============================================================
# POWER (P/Q)
# ============================================================

Base.@kwdef mutable struct SOCWRPowerWorkspace <: AbstractADMMWorkspace
    partitions::Partitions{PartitionSOCWR}
    id::Vector{Int}

    interzonal::Vector{Int}

    p_flow_consensus::Vector{Float64}
    q_flow_consensus::Vector{Float64}

    p_flow_consensus_prev::Vector{Float64}
    q_flow_consensus_prev::Vector{Float64}

    duals_p::Matrix{Float64}
    duals_q::Matrix{Float64}

    rho_p::Base.RefValue{Float64}
    rho_q::Base.RefValue{Float64}

    prev_rho_p::Base.RefValue{Float64}
    prev_rho_q::Base.RefValue{Float64}

    branch_counts::Vector{Int}

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

# ============================================================
# HYBRID (V² + P/Q + R/T)
# ============================================================

Base.@kwdef mutable struct SOCWRHybridWorkspace <: AbstractADMMWorkspace
    partitions::Partitions{PartitionSOCWR}
    id::Vector{Int}

    interzonal::Vector{Int}

    # BUS
    v2_consensus::Vector{Float64}
    v2_consensus_prev::Vector{Float64}

    # PQ
    p_flow_consensus::Vector{Float64}
    q_flow_consensus::Vector{Float64}

    p_flow_consensus_prev::Vector{Float64}
    q_flow_consensus_prev::Vector{Float64}

    # RT
    r_consensus::Vector{Float64}
    t_consensus::Vector{Float64}

    r_consensus_prev::Vector{Float64}
    t_consensus_prev::Vector{Float64}

    # DUALS
    duals_v2::Matrix{Float64}

    duals_p::Matrix{Float64}
    duals_q::Matrix{Float64}

    duals_r::Matrix{Float64}
    duals_t::Matrix{Float64}

    # RHOS
    rho_v2::Base.RefValue{Float64}
    rho_p::Base.RefValue{Float64}
    rho_q::Base.RefValue{Float64}
    rho_r::Base.RefValue{Float64}
    rho_t::Base.RefValue{Float64}

    prev_rho_v2::Base.RefValue{Float64}
    prev_rho_p::Base.RefValue{Float64}
    prev_rho_q::Base.RefValue{Float64}
    prev_rho_r::Base.RefValue{Float64}
    prev_rho_t::Base.RefValue{Float64}

    bus_counts::Vector{Int}
    branch_counts::Vector{Int}

    r_primal::Float64 = 0.0
    r_dual::Float64 = 0.0
    eps_primal::Float64 = 0.0
    eps_dual::Float64 = 0.0

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

    hist::ADMMHistory = ADMMHistory()
    iter::Int = 0
end

# ============================================================
# RHO HELPERS
# ============================================================

_get_rhos(ws::SOCWRStateWorkspace) = (
    rho_v2 = ws.rho_v2[],
    rho_r = ws.rho_r[],
    rho_t = ws.rho_t[],
)

_get_rhos(ws::SOCWRPowerWorkspace) = (
    rho_p = ws.rho_p[],
    rho_q = ws.rho_q[],
)

_get_rhos(ws::SOCWRHybridWorkspace) = (
    rho_v2 = ws.rho_v2[],
    rho_p = ws.rho_p[],
    rho_q = ws.rho_q[],
    rho_r = ws.rho_r[],
    rho_t = ws.rho_t[],
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

    rho_v2 = Ref(config.rho_v2)
    rho_r = Ref(config.rho_r)
    rho_t = Ref(config.rho_t)

    return SOCWRStateWorkspace(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        v2_consensus = ones(n_shared),
        v2_consensus_prev = ones(n_shared),
        r_consensus = zeros(n_inter),
        t_consensus = zeros(n_inter),
        r_consensus_prev = zeros(n_inter),
        t_consensus_prev = zeros(n_inter),
        duals_v2 = zeros(n_shared, n_zones),
        duals_r = zeros(n_inter, n_zones),
        duals_t = zeros(n_inter, n_zones),
        rho_v2 = rho_v2,
        rho_r = rho_r,
        rho_t = rho_t,
        prev_rho_v2 = Ref(rho_v2[]),
        prev_rho_r = Ref(rho_r[]),
        prev_rho_t = Ref(rho_t[]),
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

    rho_p = Ref(config.rho_p)
    rho_q = Ref(config.rho_q)

    return SOCWRPowerWorkspace(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        p_flow_consensus = zeros(n_inter),
        q_flow_consensus = zeros(n_inter),
        p_flow_consensus_prev = zeros(n_inter),
        q_flow_consensus_prev = zeros(n_inter),
        duals_p = zeros(n_inter, n_zones),
        duals_q = zeros(n_inter, n_zones),
        rho_p = rho_p,
        rho_q = rho_q,
        prev_rho_p = Ref(rho_p[]),
        prev_rho_q = Ref(rho_q[]),
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

    rho_v2 = Ref(config.rho_v2)
    rho_p = Ref(config.rho_p)
    rho_q = Ref(config.rho_q)
    rho_r = Ref(config.rho_r)
    rho_t = Ref(config.rho_t)

    return SOCWRHybridWorkspace(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        v2_consensus = ones(n_shared),
        v2_consensus_prev = ones(n_shared),
        p_flow_consensus = zeros(n_inter),
        q_flow_consensus = zeros(n_inter),
        p_flow_consensus_prev = zeros(n_inter),
        q_flow_consensus_prev = zeros(n_inter),
        r_consensus = zeros(n_inter),
        t_consensus = zeros(n_inter),
        r_consensus_prev = zeros(n_inter),
        t_consensus_prev = zeros(n_inter),
        duals_v2 = zeros(n_shared, n_zones),
        duals_p = zeros(n_inter, n_zones),
        duals_q = zeros(n_inter, n_zones),
        duals_r = zeros(n_inter, n_zones),
        duals_t = zeros(n_inter, n_zones),
        rho_v2 = rho_v2,
        rho_p = rho_p,
        rho_q = rho_q,
        rho_r = rho_r,
        rho_t = rho_t,
        prev_rho_v2 = Ref(rho_v2[]),
        prev_rho_p = Ref(rho_p[]),
        prev_rho_q = Ref(rho_q[]),
        prev_rho_r = Ref(rho_r[]),
        prev_rho_t = Ref(rho_t[]),
        bus_counts = zeros(Int, n_shared),
        branch_counts = zeros(Int, n_inter),
    )
end
