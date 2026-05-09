# ============================================================
# Workspace structs
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

Base.@kwdef mutable struct SOCWRStateWorkspaceALADIN <: AbstractALADINWorkspace
    partitions::Partitions{PartitionSOCWR}
    id::Vector{Int}
    interzonal::Vector{Int}

    zones::Vector{AladinZoneDataSOCWR}
    bus_pairs::Vector{SharedBusPairSOCWR}
    interzonal_line_pairs::Vector{InterzonalLinePairSOCWR}

    vars::NamedTuple{(:v2, :r, :t), Tuple{LambdaState, LambdaState, LambdaState}}

    qp_optimizer::Any

    rho::Base.RefValue{Float64}
    mu::Base.RefValue{Float64}

    alpha1::Float64 = 1.0
    alpha2::Float64 = 1.0
    alpha3::Float64 = 1.0
    alpha2_init::Float64 = 1.0
    line_search::Bool = true

    sigma_v2::Float64 = 10.0
    sigma_r::Float64 = 5.0
    sigma_t::Float64 = 5.0
    sigma_theta::Float64 = 1.0
    sigma_pg::Float64 = 1.0
    sigma_qg::Float64 = 1.0
    sigma_other::Float64 = 1.0

    r_primal::Float64 = 0.0
    r_progress::Float64 = 0.0
    lambda_step::Float64 = 0.0
    eps_primal::Float64 = 0.0
    eps_progress::Float64 = 0.0
    iter::Int = 0
    hist::ALADINHistory = ALADINHistory()
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

Base.@kwdef mutable struct SOCWRPowerWorkspaceALADIN <: AbstractALADINWorkspace
    partitions::Partitions{PartitionSOCWR}
    id::Vector{Int}
    interzonal::Vector{Int}

    zones::Vector{AladinZoneDataSOCWR}
    interzonal_line_pairs::Vector{InterzonalLinePairSOCWR}

    vars::NamedTuple{(:p_flow, :q_flow), Tuple{LambdaState, LambdaState}}

    qp_optimizer::Any

    rho::Base.RefValue{Float64}
    mu::Base.RefValue{Float64}

    alpha1::Float64 = 1.0
    alpha2::Float64 = 1.0
    alpha3::Float64 = 1.0
    alpha2_init::Float64 = 1.0
    line_search::Bool = true

    sigma_v2::Float64 = 5.0
    sigma_r::Float64 = 5.0
    sigma_t::Float64 = 5.0
    sigma_theta::Float64 = 1.0
    sigma_pg::Float64 = 1.0
    sigma_qg::Float64 = 1.0
    sigma_other::Float64 = 1.0
    sigma_pflow::Float64 = 10.0
    sigma_qflow::Float64 = 10.0

    proximal_mode::Symbol = :primal

    r_primal::Float64 = 0.0
    r_progress::Float64 = 0.0
    lambda_step::Float64 = 0.0
    eps_primal::Float64 = 0.0
    eps_progress::Float64 = 0.0
    iter::Int = 0
    hist::ALADINHistory = ALADINHistory()
end

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

Base.@kwdef mutable struct SOCWRHybridWorkspaceALADIN <: AbstractALADINWorkspace
    partitions::Partitions{PartitionSOCWR}
    id::Vector{Int}
    interzonal::Vector{Int}

    zones::Vector{AladinZoneDataSOCWR}
    bus_pairs::Vector{SharedBusPairSOCWR}
    interzonal_line_pairs::Vector{InterzonalLinePairSOCWR}
    split_bus_pairs::Vector{SplitBusPairSOCWR} = SplitBusPairSOCWR[]
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}} =
        Vector{Dict{Int, JuMP.VariableRef}}()
    q_inj_vars::Vector{Dict{Int, JuMP.VariableRef}} =
        Vector{Dict{Int, JuMP.VariableRef}}()

    vars::NamedTuple{
        (:v2, :p_flow, :q_flow, :r, :t, :p_inj, :q_inj),
        Tuple{
            LambdaState, LambdaState, LambdaState, LambdaState, LambdaState,
            LambdaState, LambdaState,
        },
    }

    qp_optimizer::Any

    rho::Base.RefValue{Float64}
    mu::Base.RefValue{Float64}

    alpha1::Float64 = 1.0
    alpha2::Float64 = 1.0
    alpha3::Float64 = 1.0
    alpha2_init::Float64 = 1.0
    line_search::Bool = true

    sigma_v2::Float64 = 10.0
    sigma_r::Float64 = 5.0
    sigma_t::Float64 = 5.0
    sigma_theta::Float64 = 1.0
    sigma_pg::Float64 = 1.0
    sigma_qg::Float64 = 1.0
    sigma_other::Float64 = 1.0

    r_primal::Float64 = 0.0
    r_progress::Float64 = 0.0
    lambda_step::Float64 = 0.0
    eps_primal::Float64 = 0.0
    eps_progress::Float64 = 0.0
    iter::Int = 0
    hist::ALADINHistory = ALADINHistory()
end

# ============================================================
# Rho accessors
# ============================================================

_get_rhos(ws::SOCWRStateWorkspaceALADIN) = (rho = ws.rho[],)
_get_rhos(ws::SOCWRPowerWorkspaceALADIN) = (rho = ws.rho[],)
_get_rhos(ws::SOCWRHybridWorkspaceALADIN) = (rho = ws.rho[],)

# ============================================================
# Builders
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function build_workspace(
    sys::DistributedSystem,
    cfg::ALADINConfig,
    optimizer_factory,
    ::ALADIN{State, SOCWR},
)
    partitions = build_partitions(PartitionSOCWR, sys)
    id = partitions.zones
    inter = sys.interzonal_branches

    zones, bus_pairs, interzonal_line_pairs =
        init_state_partition!(partitions, sys, optimizer_factory, cfg)
    n_bus = length(bus_pairs)
    n_edge = length(interzonal_line_pairs)

    ws = SOCWRStateWorkspaceALADIN(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        zones = zones,
        bus_pairs = bus_pairs,
        interzonal_line_pairs = interzonal_line_pairs,
        vars = (
            v2 = LambdaState(n_bus),
            r = LambdaState(n_edge),
            t = LambdaState(n_edge),
        ),
        qp_optimizer = _resolve_qp_optimizer(optimizer_factory),
        rho = Ref(cfg.rho),
        mu = Ref(cfg.mu),
        alpha1 = cfg.alpha1,
        alpha2 = cfg.alpha2,
        alpha3 = cfg.alpha3,
        alpha2_init = cfg.alpha2,
    )
    refresh_rho!(ws)
    return ws
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function build_workspace(
    sys::DistributedSystem,
    cfg::ALADINConfig,
    optimizer_factory,
    ::ALADIN{Power, SOCWR},
)
    partitions = build_partitions(PartitionSOCWR, sys)
    id = partitions.zones
    inter = sys.interzonal_branches

    zones, interzonal_line_pairs =
        init_power_partition!(partitions, sys, optimizer_factory, cfg)
    n_edge = length(interzonal_line_pairs)

    ws = SOCWRPowerWorkspaceALADIN(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        zones = zones,
        interzonal_line_pairs = interzonal_line_pairs,
        vars = (
            p_flow = LambdaState(n_edge),
            q_flow = LambdaState(n_edge),
        ),
        qp_optimizer = _resolve_qp_optimizer(optimizer_factory),
        rho = Ref(cfg.rho),
        mu = Ref(cfg.mu),
        alpha1 = cfg.alpha1,
        alpha2 = cfg.alpha2,
        alpha3 = cfg.alpha3,
        alpha2_init = cfg.alpha2,
        sigma_pflow = cfg.sigma_pflow,
        sigma_qflow = cfg.sigma_qflow,
        proximal_mode = cfg.proximal_mode,
    )
    refresh_rho!(ws)
    return ws
end

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

function build_workspace(
    sys::DistributedSystem,
    cfg::ALADINConfig,
    optimizer_factory,
    ::ALADIN{Hybrid, SOCWR},
)
    partitions = build_partitions(PartitionSOCWR, sys)
    id = partitions.zones
    inter = sys.interzonal_branches
    sectioning = cfg.sectioning

    if sectioning isa TieLine
        zones, bus_pairs, interzonal_line_pairs =
            init_hybrid_partition!(partitions, sys, optimizer_factory, cfg, sectioning)
        split_pairs = SplitBusPairSOCWR[]
        p_inj_vars = Vector{Dict{Int, JuMP.VariableRef}}()
        q_inj_vars = Vector{Dict{Int, JuMP.VariableRef}}()
    elseif sectioning isa SplitBus
        zones, bus_pairs, interzonal_line_pairs, split_pairs, p_inj_vars, q_inj_vars =
            init_hybrid_partition!(partitions, sys, optimizer_factory, cfg, sectioning)
    else
        error("ALADIN{Hybrid, SOCWR}: unsupported sectioning $(sectioning)")
    end

    n_bus = length(bus_pairs)
    n_edge = length(interzonal_line_pairs)
    n_split = length(split_pairs)

    ws = SOCWRHybridWorkspaceALADIN(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        zones = zones,
        bus_pairs = bus_pairs,
        interzonal_line_pairs = interzonal_line_pairs,
        split_bus_pairs = split_pairs,
        p_inj_vars = p_inj_vars,
        q_inj_vars = q_inj_vars,
        vars = (
            v2 = LambdaState(n_bus),
            p_flow = LambdaState(n_edge),
            q_flow = LambdaState(n_edge),
            r = LambdaState(n_edge),
            t = LambdaState(n_edge),
            p_inj = LambdaState(n_split),
            q_inj = LambdaState(n_split),
        ),
        qp_optimizer = _resolve_qp_optimizer(optimizer_factory),
        rho = Ref(cfg.rho),
        mu = Ref(cfg.mu),
        alpha1 = cfg.alpha1,
        alpha2 = cfg.alpha2,
        alpha3 = cfg.alpha3,
        alpha2_init = cfg.alpha2,
    )
    refresh_rho!(ws)
    return ws
end
