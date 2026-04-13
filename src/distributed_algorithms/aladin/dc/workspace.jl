# ============================================================
# ALADIN-DC WORKSPACES (canonical Houska-Frasch-Diehl 2016)
# ------------------------------------------------------------
# Houska, B., Frasch, J., Diehl, M. (2016).
#   "An augmented Lagrangian based algorithm for distributed
#    non-convex optimization." SIAM Journal on Optimization.
#
# Algorithm 1:
#   Phase 1 (parallel) per zone i:
#     y_i = argmin_{x_i} f_i(x_i) + λᵀ A_i x_i + (ρ/2) ‖x_i - x_iᵏ‖²
#           s.t. local equality + inequality constraints
#     Extract g_i = ∇f_i(y_i), B_i ≈ ∇²L_i, C_i = active Jacobian
#
#   Phase 2 (centralized QP):
#     min_{Δy,s} ∑_i [g_iᵀ Δy_i + ½ Δy_iᵀ B_i Δy_i] + λᵀ s + (μ/2) ‖s‖²
#     s.t.  ∑_i A_i (y_i + Δy_i) = s        (consensus)
#           C_i Δy_i = 0                    (active-constraint preservation)
#           Δy_i,k = 0 for k ∈ active bounds (box-active preservation)
#     Extract ν = dual multiplier of the consensus row.
#
#   Update:
#     x_iᵏ⁺¹ = x_iᵏ + α₁ (y_i - x_iᵏ) + α₂ Δy_i
#     λᵏ⁺¹   = λᵏ   + α₃ (ν - λᵏ)
# ============================================================

# ------------------------------------------------------------
# Consensus pair: for tie-line sectioning (K=2) each shared
# coupling variable is referenced by exactly two zones.
# ------------------------------------------------------------

struct SharedBusPair
    consensus_idx::Int   # index into sys.full_shared_bus_ids (global)
    bus_id::Int
    zone_i::Int          # position of zone i in partitions.parts
    col_i::Int           # column of θ[bus_id] in zone i's flat variable layout
    zone_j::Int
    col_j::Int
end

struct SharedEdgePair
    consensus_idx::Int   # index into sys.interzonal_branches
    edge_global::Int     # global branch index
    zone_i::Int
    edge_local_i::Int    # position of this edge in zone_i's zone_edges vector
    zone_j::Int
    edge_local_j::Int
end

# ------------------------------------------------------------
# DC State workspace
# ------------------------------------------------------------

Base.@kwdef mutable struct DCStateWorkspaceALADIN <: AbstractALADINWorkspace
    partitions::Partitions{PartitionDC}
    id::Vector{Int}

    # Per-zone canonical ALADIN data (variable layout, sensitivities, etc.)
    zones::Vector{AladinZoneData}

    # Pairwise bus consensus constraints (one per shared bus, for K=2).
    bus_pairs::Vector{SharedBusPair}

    # Global dual: one scalar per consensus row (aligned with bus_pairs).
    lambda_bus::Vector{Float64}

    # Phase-2 QP optimizer factory.
    qp_optimizer::Any

    # Parameters (held as Refs so adaptive schemes can update them).
    rho::Base.RefValue{Float64}
    mu::Base.RefValue{Float64}

    # Line-search parameters (default 1.0 → full step).
    alpha1::Float64 = 1.0
    alpha2::Float64 = 1.0
    alpha3::Float64 = 1.0

    # HFD Σ_i diagonal weights (per variable class).
    sigma_theta_shared::Float64 = 10.0
    sigma_theta_other::Float64 = 10.0
    sigma_pg::Float64 = 1.0
    sigma_other::Float64 = 1.0

    # Convergence state — canonical HFD 2016 Theorem 4.1 criterion.
    #   r_primal   = ‖∑_i A_i y_i‖₂   (consensus violation)
    #   r_progress = ‖y_i - x_iᵏ‖₂    (Phase-1 progress / KKT proxy)
    # `lambda_step` is kept as a diagnostic only (not part of stopping rule).
    r_primal::Float64 = 0.0
    r_progress::Float64 = 0.0
    lambda_step::Float64 = 0.0
    eps_primal::Float64 = 0.0
    eps_progress::Float64 = 0.0

    # Previous lambda (for dual-residual computation).
    lambda_bus_prev::Vector{Float64}

    hist::ALADINHistory = ALADINHistory()
    iter::Int = 0
end

# ------------------------------------------------------------
# DC Power workspace
# ------------------------------------------------------------

Base.@kwdef mutable struct DCPowerWorkspaceALADIN <: AbstractALADINWorkspace
    partitions::Partitions{PartitionDC}
    id::Vector{Int}
    interzonal::Vector{Int}

    zones::Vector{AladinZoneData}
    edge_pairs::Vector{SharedEdgePair}

    lambda_edge::Vector{Float64}
    lambda_edge_prev::Vector{Float64}

    qp_optimizer::Any

    rho::Base.RefValue{Float64}
    mu::Base.RefValue{Float64}

    alpha1::Float64 = 1.0
    alpha2::Float64 = 1.0
    alpha3::Float64 = 1.0

    sigma_theta_shared::Float64 = 10.0
    sigma_theta_other::Float64 = 10.0
    sigma_pg::Float64 = 1.0
    sigma_other::Float64 = 1.0

    # Convergence state — canonical HFD 2016 Theorem 4.1 criterion.
    r_primal::Float64 = 0.0
    r_progress::Float64 = 0.0
    lambda_step::Float64 = 0.0
    eps_primal::Float64 = 0.0
    eps_progress::Float64 = 0.0

    hist::ALADINHistory = ALADINHistory()
    iter::Int = 0
end

# ------------------------------------------------------------
# DC Hybrid workspace — both bus and edge couplings simultaneously
# ------------------------------------------------------------

Base.@kwdef mutable struct DCHybridWorkspaceALADIN <: AbstractALADINWorkspace
    partitions::Partitions{PartitionDC}
    id::Vector{Int}
    interzonal::Vector{Int}

    zones::Vector{AladinZoneData}
    bus_pairs::Vector{SharedBusPair}
    edge_pairs::Vector{SharedEdgePair}

    lambda_bus::Vector{Float64}
    lambda_edge::Vector{Float64}
    lambda_bus_prev::Vector{Float64}
    lambda_edge_prev::Vector{Float64}

    qp_optimizer::Any

    rho_bus::Base.RefValue{Float64}
    rho_edge::Base.RefValue{Float64}
    mu::Base.RefValue{Float64}

    alpha1::Float64 = 1.0
    alpha2::Float64 = 1.0
    alpha3::Float64 = 1.0

    sigma_theta_shared::Float64 = 10.0
    sigma_theta_other::Float64 = 10.0
    sigma_pg::Float64 = 1.0
    sigma_other::Float64 = 1.0

    # Convergence state — canonical HFD 2016 Theorem 4.1 criterion.
    r_primal::Float64 = 0.0
    r_progress::Float64 = 0.0
    lambda_step::Float64 = 0.0
    eps_primal::Float64 = 0.0
    eps_progress::Float64 = 0.0

    # Per-component diagnostics (not part of the stopping rule)
    r_primal_bus::Float64 = 0.0
    r_primal_edge::Float64 = 0.0
    r_progress_bus::Float64 = 0.0
    r_progress_edge::Float64 = 0.0

    hist::ALADINHistory = ALADINHistory()
    iter::Int = 0
end

# ------------------------------------------------------------
# Rho helpers (mirror ADMM's _get_rhos for history/logging)
# ------------------------------------------------------------

_get_rhos(ws::DCStateWorkspaceALADIN) = (rho = ws.rho[],)
_get_rhos(ws::DCPowerWorkspaceALADIN) = (rho = ws.rho[],)
_get_rhos(ws::DCHybridWorkspaceALADIN) =
    (rho_bus = ws.rho_bus[], rho_edge = ws.rho_edge[])

# ------------------------------------------------------------
# Helpers to build the shared-pair lists
# ------------------------------------------------------------

function _build_shared_bus_pairs(
    partitions::Partitions{PartitionDC},
    zones::Vector{AladinZoneData},
)::Vector{SharedBusPair}
    # For each (bus_id, consensus_idx), collect (zone_position, local_col).
    # With tie-line sectioning each consensus_idx has exactly 2 entries.
    acc = Dict{Int, Vector{Tuple{Int, Int}}}()  # consensus_idx -> [(zone_pos, col)]
    bus_of = Dict{Int, Int}()                   # consensus_idx -> bus_id

    for (zpos, part) in enumerate(partitions.parts)
        zdata = zones[zpos]
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            cidx = part.shared_bus_consensus_idx[i]
            cidx == 0 && continue
            col = zdata.shared_theta_cols[i]
            v = get!(acc, cidx, Tuple{Int, Int}[])
            push!(v, (zpos, col))
            bus_of[cidx] = part.zone_shared_bus_ids[i]
        end
    end

    pairs = SharedBusPair[]
    for cidx in sort!(collect(keys(acc)))
        entries = acc[cidx]
        length(entries) == 2 || continue  # skip non-K=2 cases safely
        zi, ci = entries[1]
        zj, cj = entries[2]
        push!(pairs, SharedBusPair(cidx, bus_of[cidx], zi, ci, zj, cj))
    end
    return pairs
end

function _build_shared_edge_pairs(
    partitions::Partitions{PartitionDC},
    interzonal::Vector{Int},
)::Vector{SharedEdgePair}
    # For each interzonal branch, find the two zones that see it and their
    # local index within each zone's zone_edges vector.
    pairs = SharedEdgePair[]

    for (cidx, e) in enumerate(interzonal)
        entries = Tuple{Int, Int}[]  # [(zone_pos, local_edge_idx)]
        for (zpos, part) in enumerate(partitions.parts)
            part.branch_idx_full_to_zone[e] == 0 && continue
            # position of e in part.zone_edges
            local_k = findfirst(==(e), part.zone_edges)
            local_k === nothing && continue
            push!(entries, (zpos, local_k))
        end
        length(entries) == 2 || continue
        zi, li = entries[1]
        zj, lj = entries[2]
        push!(pairs, SharedEdgePair(cidx, e, zi, li, zj, lj))
    end
    return pairs
end

# ------------------------------------------------------------
# Builders
# ------------------------------------------------------------

function _build_zone_data_vector(
    partitions::Partitions{PartitionDC};
    include_edges::Bool,
)::Vector{AladinZoneData}
    parts = partitions.parts
    zones = Vector{AladinZoneData}(undef, length(parts))
    @inbounds for i in eachindex(parts)
        zones[i] = build_aladin_zone_data(parts[i]; include_edges = include_edges)
    end
    return zones
end

function _resolve_qp_optimizer(user_opt)
    # Default to OSQP with quiet settings; fall back to whatever the user gave.
    user_opt === nothing || return user_opt
    return JuMP.optimizer_with_attributes(
        OSQP.Optimizer,
        "eps_abs" => 1e-8,
        "eps_rel" => 1e-8,
        "max_iter" => 200_000,
        "verbose" => false,
        "polish" => true,
        "polish_refine_iter" => 5,
        "linsys_solver" => "qdldl",
    )
end

function build_aladin_workspace(
    sys::DistributedSystem,
    cfg::ALADINConfig,
    optimizer_factory,
    ::ALADIN{State, DC},
)
    partitions = build_partitions(PartitionDC, sys)
    id = partitions.zones

    init_state_partition!(partitions, sys, optimizer_factory, cfg.rho)

    zones = _build_zone_data_vector(partitions; include_edges = false)
    bus_pairs = _build_shared_bus_pairs(partitions, zones)
    n_pairs = length(bus_pairs)

    return DCStateWorkspaceALADIN(;
        partitions = partitions,
        id = id,
        zones = zones,
        bus_pairs = bus_pairs,
        lambda_bus = zeros(n_pairs),
        lambda_bus_prev = zeros(n_pairs),
        qp_optimizer = _resolve_qp_optimizer(optimizer_factory),
        rho = Ref(cfg.rho),
        mu = Ref(cfg.mu),
        alpha1 = cfg.alpha1,
        alpha2 = cfg.alpha2,
        alpha3 = cfg.alpha3,
        sigma_theta_shared = cfg.sigma_theta_shared,
        sigma_theta_other = cfg.sigma_theta_other,
        sigma_pg = cfg.sigma_pg,
        sigma_other = cfg.sigma_other,
    )
end

function build_aladin_workspace(
    sys::DistributedSystem,
    cfg::ALADINConfig,
    optimizer_factory,
    ::ALADIN{Power, DC},
)
    partitions = build_partitions(PartitionDC, sys)
    id = partitions.zones
    inter = sys.interzonal_branches

    init_power_partition!(partitions, sys, optimizer_factory, cfg.rho, inter)

    zones = _build_zone_data_vector(partitions; include_edges = true)
    edge_pairs = _build_shared_edge_pairs(partitions, inter)
    n_pairs = length(edge_pairs)

    return DCPowerWorkspaceALADIN(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        zones = zones,
        edge_pairs = edge_pairs,
        lambda_edge = zeros(n_pairs),
        lambda_edge_prev = zeros(n_pairs),
        qp_optimizer = _resolve_qp_optimizer(optimizer_factory),
        rho = Ref(cfg.rho),
        mu = Ref(cfg.mu),
        alpha1 = cfg.alpha1,
        alpha2 = cfg.alpha2,
        alpha3 = cfg.alpha3,
        sigma_theta_shared = cfg.sigma_theta_shared,
        sigma_theta_other = cfg.sigma_theta_other,
        sigma_pg = cfg.sigma_pg,
        sigma_other = cfg.sigma_other,
    )
end

function build_aladin_workspace(
    sys::DistributedSystem,
    cfg::ALADINConfig,
    optimizer_factory,
    ::ALADIN{Hybrid, DC},
)
    partitions = build_partitions(PartitionDC, sys)
    id = partitions.zones
    inter = sys.interzonal_branches

    init_hybrid_partition!(partitions, sys, optimizer_factory, cfg.rho, cfg.rho, inter)

    zones = _build_zone_data_vector(partitions; include_edges = true)
    bus_pairs = _build_shared_bus_pairs(partitions, zones)
    edge_pairs = _build_shared_edge_pairs(partitions, inter)
    n_bus = length(bus_pairs)
    n_edge = length(edge_pairs)

    return DCHybridWorkspaceALADIN(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        zones = zones,
        bus_pairs = bus_pairs,
        edge_pairs = edge_pairs,
        lambda_bus = zeros(n_bus),
        lambda_edge = zeros(n_edge),
        lambda_bus_prev = zeros(n_bus),
        lambda_edge_prev = zeros(n_edge),
        qp_optimizer = _resolve_qp_optimizer(optimizer_factory),
        rho_bus = Ref(cfg.rho),
        rho_edge = Ref(cfg.rho),
        mu = Ref(cfg.mu),
        alpha1 = cfg.alpha1,
        alpha2 = cfg.alpha2,
        alpha3 = cfg.alpha3,
        sigma_theta_shared = cfg.sigma_theta_shared,
        sigma_theta_other = cfg.sigma_theta_other,
        sigma_pg = cfg.sigma_pg,
        sigma_other = cfg.sigma_other,
    )
end
