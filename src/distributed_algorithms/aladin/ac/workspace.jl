# ============================================================
# Workspace structs
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

Base.@kwdef mutable struct ACStateWorkspaceALADIN <: AbstractALADINWorkspace
    partitions::Partitions{PartitionAC}
    id::Vector{Int}
    interzonal::Vector{Int}

    zones::Vector{AladinZoneDataAC}
    bus_pairs::Vector{SharedBusPairAC}

    vars::NamedTuple{(:theta, :v_mag), Tuple{LambdaState, LambdaState}}

    qp_optimizer::Any

    rho::Base.RefValue{Float64}
    mu::Base.RefValue{Float64}

    alpha1::Float64 = 1.0
    alpha2::Float64 = 1.0
    alpha3::Float64 = 1.0
    alpha2_init::Float64 = 1.0
    line_search::Bool = true

    sigma_theta::Float64 = 10.0
    sigma_vm::Float64 = 5.0
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

Base.@kwdef mutable struct ACPowerWorkspaceALADIN <: AbstractALADINWorkspace
    partitions::Partitions{PartitionAC}
    id::Vector{Int}
    interzonal::Vector{Int}

    zones::Vector{AladinZoneDataAC}
    edge_pairs::Vector{SharedEdgePairAC}
    dummy_branch_pairs::Vector{DummyBranchPairAC} = DummyBranchPairAC[]

    vars::NamedTuple{
        (:p_flow, :q_flow, :p_dummy, :q_dummy),
        Tuple{LambdaState, LambdaState, LambdaState, LambdaState},
    }

    qp_optimizer::Any

    rho::Base.RefValue{Float64}
    mu::Base.RefValue{Float64}

    alpha1::Float64 = 1.0
    alpha2::Float64 = 1.0
    alpha3::Float64 = 1.0
    alpha2_init::Float64 = 1.0
    line_search::Bool = true

    sigma_theta::Float64 = 10.0
    sigma_vm::Float64 = 5.0
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

Base.@kwdef mutable struct ACHybridWorkspaceALADIN <: AbstractALADINWorkspace
    partitions::Partitions{PartitionAC}
    id::Vector{Int}
    interzonal::Vector{Int}

    zones::Vector{AladinZoneDataAC}
    bus_pairs::Vector{SharedBusPairAC}
    edge_pairs::Vector{SharedEdgePairAC}
    split_bus_pairs::Vector{SplitBusPairAC} = SplitBusPairAC[]
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}} =
        Vector{Dict{Int, JuMP.VariableRef}}()
    q_inj_vars::Vector{Dict{Int, JuMP.VariableRef}} =
        Vector{Dict{Int, JuMP.VariableRef}}()
    dummy_branch_pairs::Vector{DummyBranchPairAC} = DummyBranchPairAC[]

    vars::NamedTuple{
        (:theta, :v_mag, :p_flow, :q_flow, :p_inj, :q_inj, :p_dummy, :q_dummy),
        Tuple{
            LambdaState, LambdaState, LambdaState, LambdaState,
            LambdaState, LambdaState, LambdaState, LambdaState,
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

    sigma_theta::Float64 = 10.0
    sigma_vm::Float64 = 5.0
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

_get_rhos(ws::ACStateWorkspaceALADIN) = (rho = ws.rho[],)
_get_rhos(ws::ACPowerWorkspaceALADIN) = (rho = ws.rho[],)
_get_rhos(ws::ACHybridWorkspaceALADIN) = (rho = ws.rho[],)

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
    ::ALADIN{State, AC},
)
    partitions = build_partitions(PartitionAC, sys)
    id = partitions.zones
    inter = sys.interzonal_branches

    zones, bus_pairs = init_state_partition!(partitions, sys, optimizer_factory, cfg)
    n_bus = length(bus_pairs)

    ws = ACStateWorkspaceALADIN(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        zones = zones,
        bus_pairs = bus_pairs,
        vars = (
            theta = LambdaState(n_bus),
            v_mag = LambdaState(n_bus),
        ),
        qp_optimizer = _resolve_qp_optimizer(optimizer_factory),
        rho = Ref(cfg.rho),
        mu = Ref(cfg.mu),
        alpha1 = cfg.alpha1,
        alpha2 = cfg.alpha2,
        alpha3 = cfg.alpha3,
        alpha2_init = cfg.alpha2,
    )
    _refresh_quad_static_ac_all!(ws.zones, ws.rho[], ws)
    return ws
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function build_workspace(
    sys::DistributedSystem,
    cfg::ALADINConfig,
    optimizer_factory,
    ::ALADIN{Power, AC},
)
    partitions = build_partitions(PartitionAC, sys)
    id = partitions.zones
    inter = sys.interzonal_branches

    zones, edge_pairs, dummy_pairs =
        init_power_partition!(partitions, sys, optimizer_factory, cfg)
    n_edge = length(edge_pairs)
    n_dummy = length(dummy_pairs)

    ws = ACPowerWorkspaceALADIN(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        zones = zones,
        edge_pairs = edge_pairs,
        dummy_branch_pairs = dummy_pairs,
        vars = (
            p_flow = LambdaState(n_edge),
            q_flow = LambdaState(n_edge),
            p_dummy = LambdaState(n_dummy),
            q_dummy = LambdaState(n_dummy),
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
    _refresh_quad_static_ac_all!(ws.zones, ws.rho[], ws)
    return ws
end

# ------------------------------------------------------------
# Hybrid
# ------------------------------------------------------------

function build_workspace(
    sys::DistributedSystem,
    cfg::ALADINConfig,
    optimizer_factory,
    ::ALADIN{Hybrid, AC},
)
    partitions = build_partitions(PartitionAC, sys)
    id = partitions.zones
    inter = sys.interzonal_branches
    sectioning = cfg.sectioning

    if sectioning isa TieLine
        zones, bus_pairs, edge_pairs =
            init_hybrid_partition!(partitions, sys, optimizer_factory, cfg, sectioning)
        split_pairs = SplitBusPairAC[]
        dummy_pairs = DummyBranchPairAC[]
        p_inj_vars = Vector{Dict{Int, JuMP.VariableRef}}()
        q_inj_vars = Vector{Dict{Int, JuMP.VariableRef}}()
    elseif sectioning isa DummyBus
        zones, bus_pairs, edge_pairs, dummy_pairs =
            init_hybrid_partition!(partitions, sys, optimizer_factory, cfg, sectioning)
        split_pairs = SplitBusPairAC[]
        p_inj_vars = Vector{Dict{Int, JuMP.VariableRef}}()
        q_inj_vars = Vector{Dict{Int, JuMP.VariableRef}}()
    elseif sectioning isa SplitBus
        zones, bus_pairs, edge_pairs, split_pairs, p_inj_vars, q_inj_vars =
            init_hybrid_partition!(partitions, sys, optimizer_factory, cfg, sectioning)
        dummy_pairs = DummyBranchPairAC[]
    else
        error("ALADIN{Hybrid, AC}: unsupported sectioning $(sectioning)")
    end

    n_bus = length(bus_pairs)
    n_edge = length(edge_pairs)
    n_split = length(split_pairs)
    n_dummy = length(dummy_pairs)

    ws = ACHybridWorkspaceALADIN(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        zones = zones,
        bus_pairs = bus_pairs,
        edge_pairs = edge_pairs,
        split_bus_pairs = split_pairs,
        p_inj_vars = p_inj_vars,
        q_inj_vars = q_inj_vars,
        dummy_branch_pairs = dummy_pairs,
        vars = (
            theta = LambdaState(n_bus),
            v_mag = LambdaState(n_bus),
            p_flow = LambdaState(n_edge),
            q_flow = LambdaState(n_edge),
            p_inj = LambdaState(n_split),
            q_inj = LambdaState(n_split),
            p_dummy = LambdaState(n_dummy),
            q_dummy = LambdaState(n_dummy),
        ),
        qp_optimizer = _resolve_qp_optimizer(optimizer_factory),
        rho = Ref(cfg.rho),
        mu = Ref(cfg.mu),
        alpha1 = cfg.alpha1,
        alpha2 = cfg.alpha2,
        alpha3 = cfg.alpha3,
        alpha2_init = cfg.alpha2,
    )
    _refresh_quad_static_ac_all!(ws.zones, ws.rho[], ws)
    return ws
end
