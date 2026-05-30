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

Base.@kwdef mutable struct DCPowerWorkspace{S <: Sectioning} <: AbstractADMMWorkspace
    partitions::Partitions{PartitionDC}
    id::Vector{Int}

    sectioning::S

    interzonal::Vector{Int}
    branch_pairs::Vector{NTuple{4, Int}} = NTuple{4, Int}[]
    pair_counts::Vector{Int} = Int[]

    # `p_flow` carries tie-line consensus; `p_dummy` carries dummy-line
    # consensus (Method 2). Exactly one is non-empty depending on `S`.
    vars::NamedTuple{(:p_flow, :p_dummy), Tuple{VarState, VarState}}

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

Base.@kwdef mutable struct DCHybridWorkspace{S <: Sectioning} <: AbstractADMMWorkspace
    partitions::Partitions{PartitionDC}
    id::Vector{Int}

    sectioning::S

    interzonal::Vector{Int}
    branch_pairs::Vector{NTuple{4, Int}} = NTuple{4, Int}[]
    pair_counts::Vector{Int} = Int[]

    bus_pairs::Vector{NTuple{3, Int}} = NTuple{3, Int}[]
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}} =
        Vector{Dict{Int, JuMP.VariableRef}}()

    # SplitBus N-way injection: one copy per (split bus, zone). `inj_groups[g] =
    # (bus_id, [copy_idx...])`; `inj_copy_zone[c]` is copy c's zone position.
    inj_groups::Vector{Tuple{Int, Vector{Int}}} = Tuple{Int, Vector{Int}}[]
    inj_copy_zone::Vector{Int} = Int[]

    vars::NamedTuple{(:theta, :p_flow, :p_dummy, :p_inj),
        Tuple{VarState, VarState, VarState, VarState}}

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

    init_state_partition!(partitions, sys, optimizer_factory)

    vars = (
        theta = VarState(n_shared, n_zones, config.rho_theta),
    )

    ws = DCStateWorkspace(;
        partitions = partitions,
        id = id,
        vars = vars,
        bus_counts = zeros(Int, n_shared),
    )
    refresh_rho!(ws)
    return ws
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

    sectioning = config.sectioning
    if sectioning isa TieLine
        n_p_flow = n_inter
        n_p_dummy = 0
    elseif sectioning isa DummyBus
        n_p_flow = 0
        n_p_dummy = n_pairs
    else
        error("ADMM{Power, DC}: unsupported sectioning $(sectioning)")
    end

    init_power_partition!(partitions, sys, optimizer_factory, inter)

    p_flow = VarState(n_p_flow, n_zones, config.rho_p)
    vars = (
        p_flow = p_flow,
        # p_dummy shares p_flow's rho refs so adaptive rho/refresh kernels
        # update both consistently with a single schedule.
        p_dummy = VarState(n_p_dummy, n_zones, p_flow.rho, p_flow.prev_rho),
    )

    ws = DCPowerWorkspace(;
        partitions = partitions,
        id = id,
        sectioning = sectioning,
        interzonal = inter,
        branch_pairs = pairs,
        pair_counts = fill(2, n_pairs),
        vars = vars,
        branch_counts = zeros(Int, n_inter),
    )
    refresh_rho!(ws)
    return ws
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
    n_bus_pairs = length(bpairs)

    p_inj_vars = [Dict{Int, JuMP.VariableRef}() for _ in 1:n_zones]

    inj_groups = Tuple{Int, Vector{Int}}[]
    inj_copy_zone = Int[]
    n_p_inj = 0

    sectioning = config.sectioning
    if sectioning isa TieLine
        n_p_flow = n_inter
        n_p_dummy = 0
        init_hybrid_partition!(
            partitions, sys, optimizer_factory, sectioning, inter,
        )
    elseif sectioning isa DummyBus
        n_p_flow = 0
        n_p_dummy = n_pairs
        init_hybrid_partition!(
            partitions, sys, optimizer_factory, sectioning, inter,
        )
    elseif sectioning isa SplitBus
        n_p_flow = 0
        n_p_dummy = 0
        for (b_id, zones) in sys.bus_groups
            cidxs = Int[]
            for z in zones
                n_p_inj += 1
                push!(cidxs, n_p_inj)
                push!(inj_copy_zone, z)
            end
            push!(inj_groups, (b_id, cidxs))
        end
        init_hybrid_partition!(
            partitions, sys, optimizer_factory, sectioning, inter,
            bpairs, p_inj_vars,
        )
    else
        error("ADMM{Hybrid, DC}: unsupported sectioning $(sectioning)")
    end

    p_flow = VarState(n_p_flow, n_zones, config.rho_p)
    vars = (
        theta = VarState(n_shared, n_zones, config.rho_theta),
        p_flow = p_flow,
        p_dummy = VarState(n_p_dummy, n_zones, p_flow.rho, p_flow.prev_rho),
        p_inj = VarState(n_p_inj, n_zones, config.rho_p),
    )

    ws = DCHybridWorkspace(;
        partitions = partitions,
        id = id,
        sectioning = sectioning,
        interzonal = inter,
        branch_pairs = pairs,
        pair_counts = fill(2, n_pairs),
        bus_pairs = bpairs,
        p_inj_vars = p_inj_vars,
        inj_groups = inj_groups,
        inj_copy_zone = inj_copy_zone,
        vars = vars,
        bus_counts = zeros(Int, n_shared),
        branch_counts = zeros(Int, n_inter),
    )
    refresh_rho!(ws)
    return ws
end
