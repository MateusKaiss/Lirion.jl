# ============================================================
# Workspace structs
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

Base.@kwdef mutable struct SOCWRStateWorkspace <: AbstractADMMWorkspace
    partitions::Partitions{PartitionSOCWR}
    id::Vector{Int}

    interzonal::Vector{Int}

    vars::NamedTuple{(:v2, :r, :t), Tuple{VarState, VarState, VarState}}

    bus_counts::Vector{Int}
    branch_counts::Vector{Int}

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

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

Base.@kwdef mutable struct SOCWRPowerWorkspace{S <: Sectioning} <: AbstractADMMWorkspace
    partitions::Partitions{PartitionSOCWR}
    id::Vector{Int}

    sectioning::S

    interzonal::Vector{Int}
    branch_pairs::Vector{NTuple{4, Int}} = NTuple{4, Int}[]
    pair_counts::Vector{Int} = Int[]

    # `p_flow`/`q_flow` carry tie-line consensus; `p_dummy`/`q_dummy`
    # carry dummy-line consensus (Method 2). Exactly one pair is non-empty
    # depending on `S`.
    vars::NamedTuple{
        (:p_flow, :q_flow, :p_dummy, :q_dummy),
        Tuple{VarState, VarState, VarState, VarState},
    }

    branch_counts::Vector{Int}

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

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

Base.@kwdef mutable struct SOCWRHybridWorkspace{S <: Sectioning} <: AbstractADMMWorkspace
    partitions::Partitions{PartitionSOCWR}
    id::Vector{Int}

    sectioning::S

    interzonal::Vector{Int}
    branch_pairs::Vector{NTuple{4, Int}} = NTuple{4, Int}[]
    pair_counts::Vector{Int} = Int[]

    bus_pairs::Vector{NTuple{3, Int}} = NTuple{3, Int}[]
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}} =
        Vector{Dict{Int, JuMP.VariableRef}}()
    q_inj_vars::Vector{Dict{Int, JuMP.VariableRef}} =
        Vector{Dict{Int, JuMP.VariableRef}}()

    vars::NamedTuple{
        (:v2, :p_flow, :q_flow, :p_dummy, :q_dummy, :r, :t, :p_inj, :q_inj),
        Tuple{
            VarState, VarState, VarState, VarState, VarState,
            VarState, VarState, VarState, VarState,
        },
    }

    bus_counts::Vector{Int}
    branch_counts::Vector{Int}

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
# Rho accessors
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
# Builders
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function build_workspace(sys, config, optimizer_factory, ::ADMM{Power, SOCWR})
    partitions = build_partitions(PartitionSOCWR, sys)
    id = partitions.zones
    inter = sys.interzonal_branches
    brpairs = sys.branch_pairs
    n_zones = length(partitions)
    n_inter = length(inter)
    n_pairs = length(brpairs)

    sectioning = config.sectioning
    if sectioning isa TieLine
        n_flow = n_inter
        n_dummy = 0
    elseif sectioning isa DummyBus
        n_flow = 0
        n_dummy = n_pairs
    else
        error("ADMM{Power, SOCWR}: unsupported sectioning $(sectioning)")
    end

    init_power_partition!(
        partitions,
        sys,
        optimizer_factory,
        config.rho_p,
        config.rho_q,
        inter,
        brpairs,
    )

    p_flow = VarState(n_flow, n_zones, config.rho_p)
    q_flow = VarState(n_flow, n_zones, config.rho_q)
    vars = (
        p_flow = p_flow,
        q_flow = q_flow,
        p_dummy = VarState(n_dummy, n_zones, p_flow.rho, p_flow.prev_rho),
        q_dummy = VarState(n_dummy, n_zones, q_flow.rho, q_flow.prev_rho),
    )

    return SOCWRPowerWorkspace(;
        partitions = partitions,
        id = id,
        sectioning = sectioning,
        interzonal = inter,
        branch_pairs = brpairs,
        pair_counts = fill(2, n_pairs),
        vars = vars,
        branch_counts = zeros(Int, n_inter),
    )
end

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

function build_workspace(sys, config, optimizer_factory, ::ADMM{Hybrid, SOCWR})
    partitions = build_partitions(PartitionSOCWR, sys)
    id = partitions.zones
    inter = sys.interzonal_branches
    brpairs = sys.branch_pairs
    bpairs = sys.bus_pairs
    n_zones = length(partitions)
    n_shared = length(sys.full_shared_bus_ids)
    n_inter = length(inter)
    n_pairs = length(brpairs)
    n_bus_pairs = length(bpairs)

    p_inj_vars = [Dict{Int, JuMP.VariableRef}() for _ in 1:n_zones]
    q_inj_vars = [Dict{Int, JuMP.VariableRef}() for _ in 1:n_zones]

    sectioning = config.sectioning
    if sectioning isa TieLine
        n_flow = n_inter
        n_dummy = 0
        n_rt = n_inter
        init_hybrid_partition!(
            partitions, sys, optimizer_factory,
            config.rho_v2, config.rho_p, config.rho_q, config.rho_r, config.rho_t,
            sectioning, inter,
        )
    elseif sectioning isa DummyBus
        n_flow = 0
        n_dummy = n_pairs
        n_rt = 0
        init_hybrid_partition!(
            partitions, sys, optimizer_factory,
            config.rho_v2, config.rho_p, config.rho_q, config.rho_r, config.rho_t,
            sectioning, inter, brpairs,
        )
    elseif sectioning isa SplitBus
        n_flow = 0
        n_dummy = 0
        n_rt = 0
        init_hybrid_partition!(
            partitions, sys, optimizer_factory,
            config.rho_v2, config.rho_p, config.rho_q, config.rho_r, config.rho_t,
            sectioning, inter, bpairs, p_inj_vars, q_inj_vars,
        )
    else
        error("ADMM{Hybrid, SOCWR}: unsupported sectioning $(sectioning)")
    end

    p_flow = VarState(n_flow, n_zones, config.rho_p)
    q_flow = VarState(n_flow, n_zones, config.rho_q)
    vars = (
        v2 = VarState(n_shared, n_zones, config.rho_v2; init_value = 1.0),
        p_flow = p_flow,
        q_flow = q_flow,
        p_dummy = VarState(n_dummy, n_zones, p_flow.rho, p_flow.prev_rho),
        q_dummy = VarState(n_dummy, n_zones, q_flow.rho, q_flow.prev_rho),
        r = VarState(n_rt, n_zones, config.rho_r),
        t = VarState(n_rt, n_zones, config.rho_t),
        p_inj = VarState(n_bus_pairs, n_zones, config.rho_p),
        q_inj = VarState(n_bus_pairs, n_zones, config.rho_q),
    )

    return SOCWRHybridWorkspace(;
        partitions = partitions,
        id = id,
        sectioning = sectioning,
        interzonal = inter,
        branch_pairs = brpairs,
        pair_counts = fill(2, n_pairs),
        bus_pairs = bpairs,
        p_inj_vars = p_inj_vars,
        q_inj_vars = q_inj_vars,
        vars = vars,
        bus_counts = zeros(Int, n_shared),
        branch_counts = zeros(Int, n_inter),
    )
end
