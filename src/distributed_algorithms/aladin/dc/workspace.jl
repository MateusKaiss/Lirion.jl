# ============================================================
# Workspace structs
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

Base.@kwdef mutable struct DCStateWorkspaceALADIN <: AbstractALADINWorkspace
    partitions::Partitions{PartitionDC}
    id::Vector{Int}

    zones::Vector{AladinZoneData}
    bus_pairs::Vector{SharedBusPair}

    vars::NamedTuple{(:theta,), Tuple{LambdaState}}

    qp_optimizer::Any

    rho::Base.RefValue{Float64}
    mu::Base.RefValue{Float64}

    alpha1::Float64 = 1.0
    alpha2::Float64 = 1.0
    alpha3::Float64 = 1.0
    alpha2_init::Float64 = 1.0
    line_search::Bool = true

    sigma_theta_shared::Float64 = 10.0
    sigma_pg::Float64 = 1.0
    sigma_other::Float64 = 1.0
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
# Power
# ------------------------------------------------------------

Base.@kwdef mutable struct DCPowerWorkspaceALADIN <: AbstractALADINWorkspace
    partitions::Partitions{PartitionDC}
    id::Vector{Int}
    interzonal::Vector{Int}

    zones::Vector{AladinZoneData}
    interzonal_line_pairs::Vector{InterzonalLinePair}
    dummy_line_pairs::Vector{DummyLinePairDC} = DummyLinePairDC[]

    vars::NamedTuple{(:p_flow, :p_dummy), Tuple{LambdaState, LambdaState}}

    qp_optimizer::Any

    rho::Base.RefValue{Float64}
    mu::Base.RefValue{Float64}

    alpha1::Float64 = 1.0
    alpha2::Float64 = 1.0
    alpha3::Float64 = 1.0
    alpha2_init::Float64 = 1.0
    line_search::Bool = true

    sigma_theta_shared::Float64 = 10.0
    sigma_pg::Float64 = 1.0
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

Base.@kwdef mutable struct DCHybridWorkspaceALADIN <: AbstractALADINWorkspace
    partitions::Partitions{PartitionDC}
    id::Vector{Int}
    interzonal::Vector{Int}

    zones::Vector{AladinZoneData}
    bus_pairs::Vector{SharedBusPair}
    interzonal_line_pairs::Vector{InterzonalLinePair}
    split_bus_pairs::Vector{SplitBusPair} = SplitBusPair[]
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}} =
        Vector{Dict{Int, JuMP.VariableRef}}()
    dummy_line_pairs::Vector{DummyLinePairDC} = DummyLinePairDC[]

    vars::NamedTuple{
        (:theta, :p_flow, :p_inj, :p_dummy),
        Tuple{LambdaState, LambdaState, LambdaState, LambdaState},
    }

    qp_optimizer::Any

    rho_bus::Base.RefValue{Float64}
    rho_edge::Base.RefValue{Float64}
    rho_inj::Base.RefValue{Float64}
    mu::Base.RefValue{Float64}

    alpha1::Float64 = 1.0
    alpha2::Float64 = 1.0
    alpha3::Float64 = 1.0
    alpha2_init::Float64 = 1.0
    line_search::Bool = true

    sigma_theta_shared::Float64 = 10.0
    sigma_pg::Float64 = 1.0
    sigma_other::Float64 = 1.0
    proximal_mode::Symbol = :primal

    r_primal_bus::Float64 = 0.0
    r_primal_edge::Float64 = 0.0
    r_progress_bus::Float64 = 0.0
    r_progress_edge::Float64 = 0.0

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

_get_rhos(ws::DCStateWorkspaceALADIN) = (rho = ws.rho[],)
_get_rhos(ws::DCPowerWorkspaceALADIN) = (rho = ws.rho[],)
_get_rhos(ws::DCHybridWorkspaceALADIN) = (
    rho_bus = ws.rho_bus[],
    rho_edge = ws.rho_edge[],
    rho_inj = ws.rho_inj[],
)

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
    ::ALADIN{State, DC},
)
    partitions = build_partitions(PartitionDC, sys)
    id = partitions.zones

    zones, bus_pairs = init_state_partition!(partitions, sys, optimizer_factory, cfg)
    n_pairs = length(bus_pairs)

    return DCStateWorkspaceALADIN(;
        partitions = partitions,
        id = id,
        zones = zones,
        bus_pairs = bus_pairs,
        vars = (theta = LambdaState(n_pairs),),
        qp_optimizer = _resolve_qp_optimizer(optimizer_factory),
        rho = Ref(cfg.rho),
        mu = Ref(cfg.mu),
        alpha1 = cfg.alpha1,
        alpha2 = cfg.alpha2,
        alpha3 = cfg.alpha3,
        alpha2_init = cfg.alpha2,
        sigma_theta_shared = cfg.sigma_theta_shared,
        sigma_pg = cfg.sigma_pg,
        sigma_other = cfg.sigma_other,
        proximal_mode = cfg.proximal_mode,
    )
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function build_workspace(
    sys::DistributedSystem,
    cfg::ALADINConfig,
    optimizer_factory,
    ::ALADIN{Power, DC},
)
    partitions = build_partitions(PartitionDC, sys)
    id = partitions.zones
    inter = sys.interzonal_branches

    zones, interzonal_line_pairs, dummy_pairs =
        init_power_partition!(partitions, sys, optimizer_factory, cfg)
    n_pairs = length(interzonal_line_pairs)
    n_dummy = length(dummy_pairs)

    return DCPowerWorkspaceALADIN(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        zones = zones,
        interzonal_line_pairs = interzonal_line_pairs,
        dummy_line_pairs = dummy_pairs,
        vars = (
            p_flow = LambdaState(n_pairs),
            p_dummy = LambdaState(n_dummy),
        ),
        qp_optimizer = _resolve_qp_optimizer(optimizer_factory),
        rho = Ref(cfg.rho),
        mu = Ref(cfg.mu),
        alpha1 = cfg.alpha1,
        alpha2 = cfg.alpha2,
        alpha3 = cfg.alpha3,
        alpha2_init = cfg.alpha2,
        sigma_theta_shared = cfg.sigma_theta_shared,
        sigma_pg = cfg.sigma_pg,
        sigma_other = cfg.sigma_other,
        sigma_pflow = cfg.sigma_pflow,
        sigma_qflow = cfg.sigma_qflow,
        proximal_mode = cfg.proximal_mode,
    )
end

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

function build_workspace(
    sys::DistributedSystem,
    cfg::ALADINConfig,
    optimizer_factory,
    ::ALADIN{Hybrid, DC},
)
    partitions = build_partitions(PartitionDC, sys)
    id = partitions.zones
    inter = sys.interzonal_branches
    sectioning = cfg.sectioning

    if sectioning isa TieLine
        zones, bus_pairs, interzonal_line_pairs =
            init_hybrid_partition!(partitions, sys, optimizer_factory, cfg, sectioning)
        split_pairs = SplitBusPair[]
        dummy_pairs = DummyLinePairDC[]
        p_inj_vars = Vector{Dict{Int, JuMP.VariableRef}}()
    elseif sectioning isa DummyBus
        zones, bus_pairs, interzonal_line_pairs, dummy_pairs =
            init_hybrid_partition!(partitions, sys, optimizer_factory, cfg, sectioning)
        split_pairs = SplitBusPair[]
        p_inj_vars = Vector{Dict{Int, JuMP.VariableRef}}()
    elseif sectioning isa SplitBus
        zones, bus_pairs, interzonal_line_pairs, split_pairs, p_inj_vars =
            init_hybrid_partition!(partitions, sys, optimizer_factory, cfg, sectioning)
        dummy_pairs = DummyLinePairDC[]
    else
        error("ALADIN{Hybrid, DC}: unsupported sectioning $(sectioning)")
    end

    n_bus = length(bus_pairs)
    n_edge = length(interzonal_line_pairs)
    n_split = length(split_pairs)
    n_dummy = length(dummy_pairs)

    return DCHybridWorkspaceALADIN(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        zones = zones,
        bus_pairs = bus_pairs,
        interzonal_line_pairs = interzonal_line_pairs,
        split_bus_pairs = split_pairs,
        p_inj_vars = p_inj_vars,
        dummy_line_pairs = dummy_pairs,
        vars = (
            theta = LambdaState(n_bus),
            p_flow = LambdaState(n_edge),
            p_inj = LambdaState(n_split),
            p_dummy = LambdaState(n_dummy),
        ),
        qp_optimizer = _resolve_qp_optimizer(optimizer_factory),
        rho_bus = Ref(cfg.rho),
        rho_edge = Ref(cfg.rho),
        rho_inj = Ref(cfg.rho),
        mu = Ref(cfg.mu),
        alpha1 = cfg.alpha1,
        alpha2 = cfg.alpha2,
        alpha3 = cfg.alpha3,
        alpha2_init = cfg.alpha2,
        sigma_theta_shared = cfg.sigma_theta_shared,
        sigma_pg = cfg.sigma_pg,
        sigma_other = cfg.sigma_other,
        proximal_mode = cfg.proximal_mode,
    )
end
