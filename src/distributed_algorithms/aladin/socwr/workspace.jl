# ============================================================
# ALADIN-SOCWR — workspaces (State / Power / Hybrid)

# ------------------------------------------------------------
# Pair descriptors
# ------------------------------------------------------------

struct SharedBusPairSOCWR
    consensus_idx::Int
    bus_id::Int
    zone_i::Int
    v2_col_i::Int
    zone_j::Int
    v2_col_j::Int
end

struct SharedEdgePairSOCWR
    consensus_idx::Int
    edge_global::Int
    zone_i::Int
    edge_local_i::Int   # position in zone_i.edge_r_cols / edge_t_cols
    zone_j::Int
    edge_local_j::Int
end

# ------------------------------------------------------------
# STATE workspace
# ------------------------------------------------------------

Base.@kwdef mutable struct SOCWRStateWorkspaceALADIN <: AbstractALADINWorkspace
    partitions::Partitions{PartitionSOCWR}
    id::Vector{Int}
    interzonal::Vector{Int}

    zones::Vector{AladinZoneDataSOCWR}

    # State consensus on v² at shared buses and on R/T at interzonal
    # branches. One scalar dual per consensus row.
    bus_pairs::Vector{SharedBusPairSOCWR}
    edge_pairs::Vector{SharedEdgePairSOCWR}

    lambda_v2::Vector{Float64}
    lambda_r::Vector{Float64}
    lambda_t::Vector{Float64}

    lambda_v2_prev::Vector{Float64}
    lambda_r_prev::Vector{Float64}
    lambda_t_prev::Vector{Float64}

    qp_optimizer::Any

    rho::Base.RefValue{Float64}
    mu::Base.RefValue{Float64}

    alpha1::Float64 = 1.0
    alpha2::Float64 = 1.0
    alpha3::Float64 = 1.0

    # Initial (cap) value of α₂ for the Armijo line search each iter.
    # Set from cfg.alpha2 at workspace build time; the effective α₂
    # used by primal_update! is stored in `alpha2` and refreshed each
    # iteration by `_line_search_step!`.
    alpha2_init::Float64 = 1.0
    line_search::Bool = true

    # Σ_i proximal weights per variable class. The shared state
    # (v², R, T) gets a stronger pull than free interior variables,
    # matching the structural sparsity of the Phase-2 consensus rows.
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

    hist::ALADINHistory = ALADINHistory()
    iter::Int = 0
end

# ------------------------------------------------------------
# POWER workspace
# ------------------------------------------------------------

Base.@kwdef mutable struct SOCWRPowerWorkspaceALADIN <: AbstractALADINWorkspace
    partitions::Partitions{PartitionSOCWR}
    id::Vector{Int}
    interzonal::Vector{Int}

    zones::Vector{AladinZoneDataSOCWR}

    edge_pairs::Vector{SharedEdgePairSOCWR}

    lambda_p::Vector{Float64}
    lambda_q::Vector{Float64}
    lambda_p_prev::Vector{Float64}
    lambda_q_prev::Vector{Float64}

    qp_optimizer::Any

    rho::Base.RefValue{Float64}
    mu::Base.RefValue{Float64}

    alpha1::Float64 = 1.0
    alpha2::Float64 = 1.0
    alpha3::Float64 = 1.0

    # Initial (cap) value of α₂ for the Armijo line search each iter.
    # Set from cfg.alpha2 at workspace build time; the effective α₂
    # used by primal_update! is stored in `alpha2` and refreshed each
    # iteration by `_line_search_step!`.
    alpha2_init::Float64 = 1.0
    line_search::Bool = true

    sigma_v2::Float64 = 5.0
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

    hist::ALADINHistory = ALADINHistory()
    iter::Int = 0
end

# ------------------------------------------------------------
# HYBRID workspace
# ------------------------------------------------------------

Base.@kwdef mutable struct SOCWRHybridWorkspaceALADIN <: AbstractALADINWorkspace
    partitions::Partitions{PartitionSOCWR}
    id::Vector{Int}
    interzonal::Vector{Int}

    zones::Vector{AladinZoneDataSOCWR}

    bus_pairs::Vector{SharedBusPairSOCWR}
    edge_pairs::Vector{SharedEdgePairSOCWR}

    lambda_v2::Vector{Float64}
    lambda_r::Vector{Float64}
    lambda_t::Vector{Float64}
    lambda_p::Vector{Float64}
    lambda_q::Vector{Float64}

    lambda_v2_prev::Vector{Float64}
    lambda_r_prev::Vector{Float64}
    lambda_t_prev::Vector{Float64}
    lambda_p_prev::Vector{Float64}
    lambda_q_prev::Vector{Float64}

    qp_optimizer::Any

    rho::Base.RefValue{Float64}
    mu::Base.RefValue{Float64}

    alpha1::Float64 = 1.0
    alpha2::Float64 = 1.0
    alpha3::Float64 = 1.0

    # Initial (cap) value of α₂ for the Armijo line search each iter.
    # Set from cfg.alpha2 at workspace build time; the effective α₂
    # used by primal_update! is stored in `alpha2` and refreshed each
    # iteration by `_line_search_step!`.
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

    hist::ALADINHistory = ALADINHistory()
    iter::Int = 0
end

# ------------------------------------------------------------
# Rho helpers
# ------------------------------------------------------------

_get_rhos(ws::SOCWRStateWorkspaceALADIN) = (rho = ws.rho[],)
_get_rhos(ws::SOCWRPowerWorkspaceALADIN) = (rho = ws.rho[],)
_get_rhos(ws::SOCWRHybridWorkspaceALADIN) = (rho = ws.rho[],)

# ------------------------------------------------------------
# Pair builders
# ------------------------------------------------------------

function _build_shared_bus_pairs_socwr(
    partitions::Partitions{PartitionSOCWR},
    zones::Vector{AladinZoneDataSOCWR},
)::Vector{SharedBusPairSOCWR}
    acc = Dict{Int, Vector{Tuple{Int, Int}}}()  # cidx -> [(zpos, v2_col)]
    bus_of = Dict{Int, Int}()

    for (zpos, part) in enumerate(partitions.parts)
        zdata = zones[zpos]
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            cidx = part.shared_bus_consensus_idx[i]
            cidx == 0 && continue
            col = zdata.shared_v2_cols[i]
            v = get!(acc, cidx, Tuple{Int, Int}[])
            push!(v, (zpos, col))
            bus_of[cidx] = part.zone_shared_bus_ids[i]
        end
    end

    pairs = SharedBusPairSOCWR[]
    for cidx in sort!(collect(keys(acc)))
        entries = acc[cidx]
        length(entries) == 2 || continue
        zi, ci = entries[1]
        zj, cj = entries[2]
        push!(pairs, SharedBusPairSOCWR(cidx, bus_of[cidx], zi, ci, zj, cj))
    end
    return pairs
end

function _build_shared_edge_pairs_socwr(
    partitions::Partitions{PartitionSOCWR},
    interzonal::Vector{Int},
)::Vector{SharedEdgePairSOCWR}
    pairs = SharedEdgePairSOCWR[]
    for (cidx, e) in enumerate(interzonal)
        entries = Tuple{Int, Int}[]
        for (zpos, part) in enumerate(partitions.parts)
            part.branch_idx_full_to_zone[e] == 0 && continue
            local_k = findfirst(==(e), part.zone_edges)
            local_k === nothing && continue
            push!(entries, (zpos, local_k))
        end
        length(entries) == 2 || continue
        zi, li = entries[1]
        zj, lj = entries[2]
        push!(pairs, SharedEdgePairSOCWR(cidx, e, zi, li, zj, lj))
    end
    return pairs
end

# ------------------------------------------------------------
# Zone-data vector builder
# ------------------------------------------------------------

function _build_zone_data_vector_socwr(
    partitions::Partitions{PartitionSOCWR};
    include_edges::Bool,
)::Vector{AladinZoneDataSOCWR}
    parts = partitions.parts
    zones = Vector{AladinZoneDataSOCWR}(undef, length(parts))
    @inbounds for i in eachindex(parts)
        zones[i] = build_aladin_zone_data_socwr(parts[i]; include_edges = include_edges)
    end
    return zones
end

# ------------------------------------------------------------
# Workspace builder — STATE
# ------------------------------------------------------------

function build_aladin_workspace(
    sys::DistributedSystem,
    cfg::ALADINConfig,
    optimizer_factory,
    ::ALADIN{State, SOCWR},
)
    partitions = build_partitions(PartitionSOCWR, sys)
    id = partitions.zones
    inter = sys.interzonal_branches

    # Route cfg.rho to every ADMM-SOCWR rho slot. ADMM's init fills
    # in the local model structure (constraints, base cost) that
    # ALADIN then overrides via `update_objective!`.
    init_state_partition!(
        partitions, sys, optimizer_factory,
        cfg.rho, cfg.rho, cfg.rho,
        inter,
    )

    zones = _build_zone_data_vector_socwr(partitions; include_edges = true)
    bus_pairs = _build_shared_bus_pairs_socwr(partitions, zones)
    edge_pairs = _build_shared_edge_pairs_socwr(partitions, inter)

    n_bus = length(bus_pairs)
    n_edge = length(edge_pairs)

    return SOCWRStateWorkspaceALADIN(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        zones = zones,
        bus_pairs = bus_pairs,
        edge_pairs = edge_pairs,
        lambda_v2 = zeros(n_bus),
        lambda_r = zeros(n_edge),
        lambda_t = zeros(n_edge),
        lambda_v2_prev = zeros(n_bus),
        lambda_r_prev = zeros(n_edge),
        lambda_t_prev = zeros(n_edge),
        qp_optimizer = _resolve_qp_optimizer(optimizer_factory),
        rho = Ref(cfg.rho),
        mu = Ref(cfg.mu),
        alpha1 = cfg.alpha1,
        alpha2 = cfg.alpha2,
        alpha3 = cfg.alpha3,
        alpha2_init = cfg.alpha2,
    )
end

# ------------------------------------------------------------
# Workspace builder — POWER
# ------------------------------------------------------------

function build_aladin_workspace(
    sys::DistributedSystem,
    cfg::ALADINConfig,
    optimizer_factory,
    ::ALADIN{Power, SOCWR},
)
    partitions = build_partitions(PartitionSOCWR, sys)
    id = partitions.zones
    inter = sys.interzonal_branches

    init_power_partition!(
        partitions, sys, optimizer_factory,
        cfg.rho, cfg.rho,
        inter,
    )

    zones = _build_zone_data_vector_socwr(partitions; include_edges = true)
    edge_pairs = _build_shared_edge_pairs_socwr(partitions, inter)
    n_edge = length(edge_pairs)

    return SOCWRPowerWorkspaceALADIN(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        zones = zones,
        edge_pairs = edge_pairs,
        lambda_p = zeros(n_edge),
        lambda_q = zeros(n_edge),
        lambda_p_prev = zeros(n_edge),
        lambda_q_prev = zeros(n_edge),
        qp_optimizer = _resolve_qp_optimizer(optimizer_factory),
        rho = Ref(cfg.rho),
        mu = Ref(cfg.mu),
        alpha1 = cfg.alpha1,
        alpha2 = cfg.alpha2,
        alpha3 = cfg.alpha3,
        alpha2_init = cfg.alpha2,
    )
end

# ------------------------------------------------------------
# Workspace builder — HYBRID
# ------------------------------------------------------------

function build_aladin_workspace(
    sys::DistributedSystem,
    cfg::ALADINConfig,
    optimizer_factory,
    ::ALADIN{Hybrid, SOCWR},
)
    partitions = build_partitions(PartitionSOCWR, sys)
    id = partitions.zones
    inter = sys.interzonal_branches

    init_hybrid_partition!(
        partitions, sys, optimizer_factory,
        cfg.rho, cfg.rho, cfg.rho, cfg.rho, cfg.rho,
        inter,
    )

    zones = _build_zone_data_vector_socwr(partitions; include_edges = true)
    bus_pairs = _build_shared_bus_pairs_socwr(partitions, zones)
    edge_pairs = _build_shared_edge_pairs_socwr(partitions, inter)

    n_bus = length(bus_pairs)
    n_edge = length(edge_pairs)

    return SOCWRHybridWorkspaceALADIN(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        zones = zones,
        bus_pairs = bus_pairs,
        edge_pairs = edge_pairs,
        lambda_v2 = zeros(n_bus),
        lambda_r = zeros(n_edge),
        lambda_t = zeros(n_edge),
        lambda_p = zeros(n_edge),
        lambda_q = zeros(n_edge),
        lambda_v2_prev = zeros(n_bus),
        lambda_r_prev = zeros(n_edge),
        lambda_t_prev = zeros(n_edge),
        lambda_p_prev = zeros(n_edge),
        lambda_q_prev = zeros(n_edge),
        qp_optimizer = _resolve_qp_optimizer(optimizer_factory),
        rho = Ref(cfg.rho),
        mu = Ref(cfg.mu),
        alpha1 = cfg.alpha1,
        alpha2 = cfg.alpha2,
        alpha3 = cfg.alpha3,
        alpha2_init = cfg.alpha2,
    )
end
