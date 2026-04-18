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
    # DummyBus (Method 2) bookkeeping: symmetric P/Q consensus on the two
    # halves of each split tie-line. Empty for TieLine.
    branch_pairs::Vector{NTuple{4, Int}} = NTuple{4, Int}[]
    pair_counts::Vector{Int} = Int[]

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
    # DummyBus (Method 2) bookkeeping: symmetric P/Q consensus on the two
    # halves of each split tie-line. Empty for TieLine/SplitBus.
    branch_pairs::Vector{NTuple{4, Int}} = NTuple{4, Int}[]
    pair_counts::Vector{Int} = Int[]

    # SplitBus (Method 3) bookkeeping: sign-flipped p_inj / q_inj consensus
    # per split bus. Empty for TieLine/DummyBus.
    bus_pairs::Vector{NTuple{3, Int}} = NTuple{3, Int}[]
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}} =
        Vector{Dict{Int, JuMP.VariableRef}}()
    q_inj_vars::Vector{Dict{Int, JuMP.VariableRef}} =
        Vector{Dict{Int, JuMP.VariableRef}}()

    vars::NamedTuple{
        (:theta, :v_mag, :p_flow, :q_flow, :p_inj, :q_inj),
        Tuple{VarState, VarState, VarState, VarState, VarState, VarState},
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
    rho_p_inj = ws.vars.p_inj.rho[],
    rho_q_inj = ws.vars.q_inj.rho[],
)

# ============================================================
# BUILDERS
# ============================================================

# ---------------- STATE ----------------

function build_workspace(sys, config, optimizer_factory, ::ADMM{State, AC})
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

function build_workspace(sys, config, optimizer_factory, ::ADMM{Power, AC})
    partitions = build_partitions(PartitionAC, sys)
    id = partitions.zones
    inter = sys.interzonal_branches
    brpairs = sys.branch_pairs
    n_zones = length(partitions)
    n_inter = length(inter)
    n_pairs = length(brpairs)
    n_total = n_inter + n_pairs

    init_power_partition!(
        partitions,
        sys,
        optimizer_factory,
        config.rho_p,
        config.rho_q,
        inter,
        brpairs,
    )

    vars = (
        p_flow = VarState(n_total, n_zones, config.rho_p),
        q_flow = VarState(n_total, n_zones, config.rho_q),
    )

    return ACPowerWorkspace(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        branch_pairs = brpairs,
        pair_counts = fill(2, n_pairs),
        vars = vars,
        p_branch_counts = zeros(Int, n_inter),
        q_branch_counts = zeros(Int, n_inter),
    )
end

# ---------------- HYBRID ----------------

function build_workspace(sys, config, optimizer_factory, ::ADMM{Hybrid, AC})
    partitions = build_partitions(PartitionAC, sys)
    id = partitions.zones
    inter = sys.interzonal_branches
    brpairs = sys.branch_pairs
    bpairs = sys.bus_pairs
    n_zones = length(partitions)
    n_shared = length(sys.full_shared_bus_ids)
    n_inter = length(inter)
    n_pairs = length(brpairs)
    n_bus_pairs = length(bpairs)
    n_total = n_inter + n_pairs

    # Per-zone p_inj / q_inj containers: populated by init_hybrid_partition! only
    # when bus_pairs is non-empty (i.e. SplitBus sectioning).
    p_inj_vars = [Dict{Int, JuMP.VariableRef}() for _ in 1:n_zones]
    q_inj_vars = [Dict{Int, JuMP.VariableRef}() for _ in 1:n_zones]

    init_hybrid_partition!(
        partitions,
        sys,
        optimizer_factory,
        config.rho_theta,
        config.rho_v,
        config.rho_p,
        config.rho_q,
        inter,
        brpairs,
        bpairs,
        p_inj_vars,
        q_inj_vars,
    )

    vars = (
        theta = VarState(n_shared, n_zones, config.rho_theta),
        v_mag = VarState(n_shared, n_zones, config.rho_v; init_value = 1.0),
        # p_flow / q_flow cover interzonal (TieLine) + branch_pairs (DummyBus)
        # in a single VarState: [1..n_inter] | [n_inter+1..n_inter+n_pairs].
        p_flow = VarState(n_total, n_zones, config.rho_p),
        q_flow = VarState(n_total, n_zones, config.rho_q),
        # Virtual injections at split buses (SplitBus Method 3).
        p_inj = VarState(n_bus_pairs, n_zones, config.rho_p),
        q_inj = VarState(n_bus_pairs, n_zones, config.rho_q),
    )

    return ACHybridWorkspace(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        branch_pairs = brpairs,
        pair_counts = fill(2, n_pairs),
        bus_pairs = bpairs,
        p_inj_vars = p_inj_vars,
        q_inj_vars = q_inj_vars,
        vars = vars,
        sin_sum = zeros(n_shared),
        cos_sum = zeros(n_shared),
        bus_counts = zeros(Int, n_shared),
        p_branch_counts = zeros(Int, n_inter),
        q_branch_counts = zeros(Int, n_inter),
    )
end
