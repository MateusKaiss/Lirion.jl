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
    globalization::Symbol = :armijo_l2
    lam_bar::Float64 = 0.0        # frozen L1 exact-penalty weight for :alg3 (0 = unset)
    lam_bar_fixed::Bool = false   # :alg3_original — λ̄ = κ̄ user-fixed (paper-literal)
    alg3_c_status::Symbol = :none # :alg3_original — outcome of the last branch (c)

    sigma_theta::Float64 = 10.0
    sigma_vm::Float64 = 5.0
    sigma_pg::Float64 = 1.0
    sigma_qg::Float64 = 1.0
    sigma_other::Float64 = 1.0

    r_primal::Float64 = 0.0
    r_step_primal::Float64 = 0.0
    lambda_step::Float64 = 0.0
    eps_primal::Float64 = 0.0
    eps_step_primal::Float64 = 0.0
    r_dual::Float64 = 0.0
    eps_dual::Float64 = 0.0
    stopping_criterion::Symbol = :boyd
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
    interzonal_line_pairs::Vector{InterzonalLinePairAC}
    dummy_line_pairs::Vector{DummyLinePairAC} = DummyLinePairAC[]

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
    globalization::Symbol = :armijo_l2
    lam_bar::Float64 = 0.0        # frozen L1 exact-penalty weight for :alg3 (0 = unset)
    lam_bar_fixed::Bool = false   # :alg3_original — λ̄ = κ̄ user-fixed (paper-literal)
    alg3_c_status::Symbol = :none # :alg3_original — outcome of the last branch (c)

    sigma_theta::Float64 = 10.0
    sigma_vm::Float64 = 5.0
    sigma_pg::Float64 = 1.0
    sigma_qg::Float64 = 1.0
    sigma_other::Float64 = 1.0
    sigma_pflow::Float64 = 10.0
    sigma_qflow::Float64 = 10.0

    proximal_mode::Symbol = :primal

    r_primal::Float64 = 0.0
    r_step_primal::Float64 = 0.0
    lambda_step::Float64 = 0.0
    eps_primal::Float64 = 0.0
    eps_step_primal::Float64 = 0.0
    r_dual::Float64 = 0.0
    eps_dual::Float64 = 0.0
    stopping_criterion::Symbol = :boyd
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
    interzonal_line_pairs::Vector{InterzonalLinePairAC}
    split_bus_pairs::Vector{SplitBusPairAC} = SplitBusPairAC[]
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}} =
        Vector{Dict{Int, JuMP.VariableRef}}()
    q_inj_vars::Vector{Dict{Int, JuMP.VariableRef}} =
        Vector{Dict{Int, JuMP.VariableRef}}()
    dummy_line_pairs::Vector{DummyLinePairAC} = DummyLinePairAC[]

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
    globalization::Symbol = :armijo_l2
    lam_bar::Float64 = 0.0        # frozen L1 exact-penalty weight for :alg3 (0 = unset)
    lam_bar_fixed::Bool = false   # :alg3_original — λ̄ = κ̄ user-fixed (paper-literal)
    alg3_c_status::Symbol = :none # :alg3_original — outcome of the last branch (c)

    sigma_theta::Float64 = 10.0
    sigma_vm::Float64 = 5.0
    sigma_pg::Float64 = 1.0
    sigma_qg::Float64 = 1.0
    sigma_other::Float64 = 1.0
    sigma_pflow::Float64 = 10.0
    sigma_qflow::Float64 = 10.0

    proximal_mode::Symbol = :primal

    r_primal::Float64 = 0.0
    r_step_primal::Float64 = 0.0
    lambda_step::Float64 = 0.0
    eps_primal::Float64 = 0.0
    eps_step_primal::Float64 = 0.0
    r_dual::Float64 = 0.0
    eps_dual::Float64 = 0.0
    stopping_criterion::Symbol = :boyd
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
        qp_optimizer = _resolve_qp_optimizer(cfg.qp_optimizer),
        stopping_criterion = cfg.stopping_criterion,
        globalization = cfg.globalization,
        lam_bar = _alg3_original_lam_bar(cfg),
        lam_bar_fixed = _alg3_original_lam_bar(cfg) > 0.0,
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
    ::ALADIN{Power, AC},
)
    partitions = build_partitions(PartitionAC, sys)
    id = partitions.zones
    inter = sys.interzonal_branches

    zones, interzonal_line_pairs, dummy_pairs =
        init_power_partition!(partitions, sys, optimizer_factory, cfg)
    n_edge = length(interzonal_line_pairs)
    n_dummy = length(dummy_pairs)

    ws = ACPowerWorkspaceALADIN(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        zones = zones,
        interzonal_line_pairs = interzonal_line_pairs,
        dummy_line_pairs = dummy_pairs,
        vars = (
            p_flow = LambdaState(n_edge),
            q_flow = LambdaState(n_edge),
            p_dummy = LambdaState(n_dummy),
            q_dummy = LambdaState(n_dummy),
        ),
        qp_optimizer = _resolve_qp_optimizer(cfg.qp_optimizer),
        stopping_criterion = cfg.stopping_criterion,
        globalization = cfg.globalization,
        lam_bar = _alg3_original_lam_bar(cfg),
        lam_bar_fixed = _alg3_original_lam_bar(cfg) > 0.0,
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
    ::ALADIN{Hybrid, AC},
)
    partitions = build_partitions(PartitionAC, sys)
    id = partitions.zones
    inter = sys.interzonal_branches
    sectioning = cfg.sectioning

    if sectioning isa TieLine
        zones, bus_pairs, interzonal_line_pairs =
            init_hybrid_partition!(partitions, sys, optimizer_factory, cfg, sectioning)
        split_pairs = SplitBusPairAC[]
        dummy_pairs = DummyLinePairAC[]
        p_inj_vars = Vector{Dict{Int, JuMP.VariableRef}}()
        q_inj_vars = Vector{Dict{Int, JuMP.VariableRef}}()
    elseif sectioning isa DummyBus
        zones, bus_pairs, interzonal_line_pairs, dummy_pairs =
            init_hybrid_partition!(partitions, sys, optimizer_factory, cfg, sectioning)
        split_pairs = SplitBusPairAC[]
        p_inj_vars = Vector{Dict{Int, JuMP.VariableRef}}()
        q_inj_vars = Vector{Dict{Int, JuMP.VariableRef}}()
    elseif sectioning isa SplitBus
        zones, bus_pairs, interzonal_line_pairs, split_pairs, p_inj_vars, q_inj_vars =
            init_hybrid_partition!(partitions, sys, optimizer_factory, cfg, sectioning)
        dummy_pairs = DummyLinePairAC[]
    else
        error("ALADIN{Hybrid, AC}: unsupported sectioning $(sectioning)")
    end

    n_bus = length(bus_pairs)
    n_edge = length(interzonal_line_pairs)
    n_split = length(split_pairs)
    n_dummy = length(dummy_pairs)

    ws = ACHybridWorkspaceALADIN(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        zones = zones,
        bus_pairs = bus_pairs,
        interzonal_line_pairs = interzonal_line_pairs,
        split_bus_pairs = split_pairs,
        p_inj_vars = p_inj_vars,
        q_inj_vars = q_inj_vars,
        dummy_line_pairs = dummy_pairs,
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
        qp_optimizer = _resolve_qp_optimizer(cfg.qp_optimizer),
        stopping_criterion = cfg.stopping_criterion,
        globalization = cfg.globalization,
        lam_bar = _alg3_original_lam_bar(cfg),
        lam_bar_fixed = _alg3_original_lam_bar(cfg) > 0.0,
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
# :alg3_original — initial λ̄
# ------------------------------------------------------------

@inline _alg3_original_lam_bar(cfg::ALADINConfig) =
    cfg.globalization === :alg3_original ? cfg.alg3_lam_bar : 0.0
