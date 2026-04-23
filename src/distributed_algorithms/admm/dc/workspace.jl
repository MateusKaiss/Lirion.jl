# ============================================================
# Workspace structs
# ============================================================

# ------------------------------------------------------------
# State
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

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

Base.@kwdef mutable struct DCPowerWorkspace <: AbstractADMMWorkspace
    partitions::Partitions{PartitionDC}
    id::Vector{Int}

    interzonal::Vector{Int}
    branch_pairs::Vector{NTuple{4, Int}} = NTuple{4, Int}[]
    pair_counts::Vector{Int} = Int[]

    vars::NamedTuple{(:p_flow,), Tuple{VarState}}

    branch_counts::Vector{Int}

    r_primal::Float64 = 0.0
    r_dual::Float64 = 0.0
    eps_primal::Float64 = 0.0
    eps_dual::Float64 = 0.0
    iter::Int = 0
    hist::ADMMHistory = ADMMHistory()
end

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

Base.@kwdef mutable struct DCHybridWorkspace <: AbstractADMMWorkspace
    partitions::Partitions{PartitionDC}
    id::Vector{Int}

    interzonal::Vector{Int}
    branch_pairs::Vector{NTuple{4, Int}} = NTuple{4, Int}[]
    pair_counts::Vector{Int} = Int[]

    bus_pairs::Vector{NTuple{3, Int}} = NTuple{3, Int}[]
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}} =
        Vector{Dict{Int, JuMP.VariableRef}}()

    vars::NamedTuple{(:theta, :p_flow, :p_inj),
        Tuple{VarState, VarState, VarState}}

    bus_counts::Vector{Int}
    branch_counts::Vector{Int}

    r_primal_bus::Float64 = 0.0
    r_dual_bus::Float64 = 0.0
    r_primal_branch::Float64 = 0.0
    r_dual_branch::Float64 = 0.0
    r_primal_bus_pair::Float64 = 0.0
    r_dual_bus_pair::Float64 = 0.0

    r_primal::Float64 = 0.0
    r_dual::Float64 = 0.0
    eps_primal::Float64 = 0.0
    eps_dual::Float64 = 0.0
    iter::Int = 0
    hist::ADMMHistory = ADMMHistory()
end

# ============================================================
# Rho accessors
# ============================================================

_get_rhos(ws::DCStateWorkspace) = (rho_theta = ws.vars.theta.rho[],)
_get_rhos(ws::DCPowerWorkspace) = (rho_p = ws.vars.p_flow.rho[],)
_get_rhos(ws::DCHybridWorkspace) = (
    rho_theta = ws.vars.theta.rho[],
    rho_p = ws.vars.p_flow.rho[],
    rho_p_inj = ws.vars.p_inj.rho[],
)

# ============================================================
# Builders
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function build_workspace(sys, config, optimizer_factory, ::ADMM{State, DC})
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

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function build_workspace(sys, config, optimizer_factory, ::ADMM{Power, DC})
    partitions = build_partitions(PartitionDC, sys)
    id = partitions.zones
    inter = sys.interzonal_branches
    pairs = sys.branch_pairs
    n_zones = length(partitions)
    n_inter = length(inter)
    n_pairs = length(pairs)
    n_total = n_inter + n_pairs

    init_power_partition!(partitions, sys, optimizer_factory, config.rho_p, inter, pairs)

    vars = (
        p_flow = VarState(n_total, n_zones, config.rho_p),
    )

    return DCPowerWorkspace(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        branch_pairs = pairs,
        pair_counts = fill(2, n_pairs),
        vars = vars,
        branch_counts = zeros(Int, n_inter),
    )
end

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

function build_workspace(sys, config, optimizer_factory, ::ADMM{Hybrid, DC})
    partitions = build_partitions(PartitionDC, sys)
    id = partitions.zones
    inter = sys.interzonal_branches
    pairs = sys.branch_pairs
    bpairs = sys.bus_pairs
    n_zones = length(partitions)
    n_shared = length(sys.full_shared_bus_ids)
    n_inter = length(inter)
    n_pairs = length(pairs)
    n_total = n_inter + n_pairs
    n_bus_pairs = length(bpairs)

    p_inj_vars = [Dict{Int, JuMP.VariableRef}() for _ in 1:n_zones]

    sectioning = config.sectioning
    if sectioning isa TieLine
        init_hybrid_partition!(
            partitions, sys, optimizer_factory,
            config.rho_theta, config.rho_p, sectioning, inter,
        )
    elseif sectioning isa DummyBus
        init_hybrid_partition!(
            partitions, sys, optimizer_factory,
            config.rho_theta, config.rho_p, sectioning, inter, pairs,
        )
    elseif sectioning isa SplitBus
        init_hybrid_partition!(
            partitions, sys, optimizer_factory,
            config.rho_theta, config.rho_p, sectioning, inter,
            bpairs, p_inj_vars,
        )
    else
        error("ADMM{Hybrid, DC}: unsupported sectioning $(sectioning)")
    end

    vars = (
        theta = VarState(n_shared, n_zones, config.rho_theta),
        p_flow = VarState(n_total, n_zones, config.rho_p),
        p_inj = VarState(n_bus_pairs, n_zones, config.rho_p),
    )

    return DCHybridWorkspace(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        branch_pairs = pairs,
        pair_counts = fill(2, n_pairs),
        bus_pairs = bpairs,
        p_inj_vars = p_inj_vars,
        vars = vars,
        bus_counts = zeros(Int, n_shared),
        branch_counts = zeros(Int, n_inter),
    )
end
