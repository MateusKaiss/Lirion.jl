# ============================================================
# SOCP WORKSPACES
# ============================================================

# ============================================================
# STATE (V² + R/T)
# ============================================================

Base.@kwdef mutable struct SOCPStateWorkspace <: AbstractADMMWorkspace
    partitions::Partitions{PartitionSOCP}
    id::Vector{Int}

    consensus_bus_idx_map::Dict{Int, Int}
    interzonal::Vector{Int}

    # BUS
    v2_consensus::Vector{Float64}
    v2_consensus_prev::Vector{Float64}

    # RT
    r_consensus::Dict{Int, Float64}
    t_consensus::Dict{Int, Float64}

    r_consensus_prev::Dict{Int, Float64}
    t_consensus_prev::Dict{Int, Float64}

    # DUALS
    duals_v2::Dict{Int, Vector{Float64}}

    duals_r::Dict{Int, Dict{Int, Float64}}
    duals_t::Dict{Int, Dict{Int, Float64}}

    # RHOS
    rho_v2::Base.RefValue{Float64}
    rho_r::Base.RefValue{Float64}
    rho_t::Base.RefValue{Float64}

    prev_rho_v2::Base.RefValue{Float64}
    prev_rho_r::Base.RefValue{Float64}
    prev_rho_t::Base.RefValue{Float64}

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

Base.@kwdef mutable struct SOCPPowerWorkspace <: AbstractADMMWorkspace
    partitions::Partitions{PartitionSOCP}
    id::Vector{Int}

    interzonal::Vector{Int}

    p_flow_consensus::Dict{Int, Float64}
    q_flow_consensus::Dict{Int, Float64}

    p_flow_consensus_prev::Dict{Int, Float64}
    q_flow_consensus_prev::Dict{Int, Float64}

    duals_p::Dict{Int, Dict{Int, Float64}}
    duals_q::Dict{Int, Dict{Int, Float64}}

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

# ============================================================
# HYBRID (V² + P/Q + R/T)
# ============================================================

Base.@kwdef mutable struct SOCPHybridWorkspace <: AbstractADMMWorkspace
    partitions::Partitions{PartitionSOCP}
    id::Vector{Int}

    consensus_bus_idx_map::Dict{Int, Int}
    interzonal::Vector{Int}

    # BUS
    v2_consensus::Vector{Float64}
    v2_consensus_prev::Vector{Float64}

    # PQ
    p_flow_consensus::Dict{Int, Float64}
    q_flow_consensus::Dict{Int, Float64}

    p_flow_consensus_prev::Dict{Int, Float64}
    q_flow_consensus_prev::Dict{Int, Float64}

    # RT
    r_consensus::Dict{Int, Float64}
    t_consensus::Dict{Int, Float64}

    r_consensus_prev::Dict{Int, Float64}
    t_consensus_prev::Dict{Int, Float64}

    # DUALS
    duals_v2::Dict{Int, Vector{Float64}}

    duals_p::Dict{Int, Dict{Int, Float64}}
    duals_q::Dict{Int, Dict{Int, Float64}}

    duals_r::Dict{Int, Dict{Int, Float64}}
    duals_t::Dict{Int, Dict{Int, Float64}}

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

_get_rhos(ws::SOCPStateWorkspace) = Dict(
    :rho_v2 => ws.rho_v2[],
    :rho_r => ws.rho_r[],
    :rho_t => ws.rho_t[],
)

_get_rhos(ws::SOCPPowerWorkspace) = Dict(
    :rho_p => ws.rho_p[],
    :rho_q => ws.rho_q[],
)

_get_rhos(ws::SOCPHybridWorkspace) = Dict(
    :rho_v2 => ws.rho_v2[],
    :rho_p => ws.rho_p[],
    :rho_q => ws.rho_q[],
    :rho_r => ws.rho_r[],
    :rho_t => ws.rho_t[],
)

# ============================================================
# BUILDERS — STATE
# ============================================================

function build_admm_workspace(sys, config, optimizer_factory, ::ADMM{State, SOCP})
    partitions = build_partitions(PartitionSOCP, sys)
    id = partitions.zones
    inter = sys.interzonal_branches

    consensus_bus_idx_map = build_consensus_bus_idx_map(sys.full_shared_bus_ids)
    n_shared = length(sys.full_shared_bus_ids)

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

    return SOCPStateWorkspace(;
        partitions = partitions,
        id = id,
        consensus_bus_idx_map = consensus_bus_idx_map,
        interzonal = inter,
        v2_consensus = ones(n_shared),
        v2_consensus_prev = zeros(n_shared),
        r_consensus = Dict(e => 0.0 for e in inter),
        t_consensus = Dict(e => 0.0 for e in inter),
        r_consensus_prev = Dict(e => 0.0 for e in inter),
        t_consensus_prev = Dict(e => 0.0 for e in inter),
        duals_v2 = Dict(z => zeros(n_shared) for z in id),
        duals_r = Dict(z => Dict{Int, Float64}() for z in id),
        duals_t = Dict(z => Dict{Int, Float64}() for z in id),
        rho_v2 = rho_v2,
        rho_r = rho_r,
        rho_t = rho_t,
        prev_rho_v2 = Ref(rho_v2[]),
        prev_rho_r = Ref(rho_r[]),
        prev_rho_t = Ref(rho_t[]),
    )
end

# ============================================================
# BUILDERS — POWER
# ============================================================

function build_admm_workspace(sys, config, optimizer_factory, ::ADMM{Power, SOCP})
    partitions = build_partitions(PartitionSOCP, sys)
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

    return SOCPPowerWorkspace(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        p_flow_consensus = Dict(e => 0.0 for e in inter),
        q_flow_consensus = Dict(e => 0.0 for e in inter),
        p_flow_consensus_prev = Dict(e => 0.0 for e in inter),
        q_flow_consensus_prev = Dict(e => 0.0 for e in inter),
        duals_p = Dict(z => Dict{Int, Float64}() for z in id),
        duals_q = Dict(z => Dict{Int, Float64}() for z in id),
        rho_p = rho_p,
        rho_q = rho_q,
        prev_rho_p = Ref(rho_p[]),
        prev_rho_q = Ref(rho_q[]),
    )
end

# ============================================================
# BUILDERS — HYBRID
# ============================================================

function build_admm_workspace(sys, config, optimizer_factory, ::ADMM{Hybrid, SOCP})
    partitions = build_partitions(PartitionSOCP, sys)
    id = partitions.zones
    inter = sys.interzonal_branches

    consensus_bus_idx_map = build_consensus_bus_idx_map(sys.full_shared_bus_ids)
    n_shared = length(sys.full_shared_bus_ids)

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

    return SOCPHybridWorkspace(;
        partitions = partitions,
        id = id,
        consensus_bus_idx_map = consensus_bus_idx_map,
        interzonal = inter,
        v2_consensus = ones(n_shared),
        v2_consensus_prev = zeros(n_shared),
        p_flow_consensus = Dict(e => 0.0 for e in inter),
        q_flow_consensus = Dict(e => 0.0 for e in inter),
        p_flow_consensus_prev = Dict(e => 0.0 for e in inter),
        q_flow_consensus_prev = Dict(e => 0.0 for e in inter),
        r_consensus = Dict(e => 0.0 for e in inter),
        t_consensus = Dict(e => 0.0 for e in inter),
        r_consensus_prev = Dict(e => 0.0 for e in inter),
        t_consensus_prev = Dict(e => 0.0 for e in inter),
        duals_v2 = Dict(z => zeros(n_shared) for z in id),
        duals_p = Dict(z => Dict{Int, Float64}() for z in id),
        duals_q = Dict(z => Dict{Int, Float64}() for z in id),
        duals_r = Dict(z => Dict{Int, Float64}() for z in id),
        duals_t = Dict(z => Dict{Int, Float64}() for z in id),
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
    )
end
