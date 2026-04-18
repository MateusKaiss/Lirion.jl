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

# SplitBus (Method 3) sign-flipped p_inj consensus: for each split bus, both
# zones hold a `p_inj` slack variable and enforce p_inj_A + p_inj_B = 0.
struct SplitBusPair
    consensus_idx::Int
    bus_id::Int
    zone_i::Int
    pinj_col_i::Int
    zone_j::Int
    pinj_col_j::Int
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

    # Duals keyed on canonical variable symbol. DC State couples only θ,
    # so `vars = (theta = LambdaState(n_pairs),)`. Legacy aliases
    # `ws.vars.theta.lambda` / `ws.vars.theta.lambda_prev` forward to
    # `vars.theta.lambda` / `vars.theta.lambda_prev` via the getproperty
    # layer in `aladin/core/components.jl`.
    vars::NamedTuple{(:theta,), Tuple{LambdaState}}

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

    r_primal::Float64 = 0.0
    r_progress::Float64 = 0.0
    lambda_step::Float64 = 0.0
    eps_primal::Float64 = 0.0
    eps_progress::Float64 = 0.0
    iter::Int = 0
    hist::ALADINHistory = ALADINHistory()
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

    # Duals: DC Power couples only branch p_flow. Legacy alias
    # `ws.vars.p_flow.lambda` → `vars.p_flow.lambda`.
    vars::NamedTuple{(:p_flow,), Tuple{LambdaState}}

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

    r_primal::Float64 = 0.0
    r_progress::Float64 = 0.0
    lambda_step::Float64 = 0.0
    eps_primal::Float64 = 0.0
    eps_progress::Float64 = 0.0
    iter::Int = 0
    hist::ALADINHistory = ALADINHistory()
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
    # SplitBus (Method 3): sign-flipped p_inj consensus per split bus.
    split_bus_pairs::Vector{SplitBusPair} = SplitBusPair[]
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}} =
        Vector{Dict{Int, JuMP.VariableRef}}()

    # Duals: DC Hybrid couples bus θ, branch p_flow, and (SplitBus) p_inj.
    vars::NamedTuple{(:theta, :p_flow, :p_inj),
        Tuple{LambdaState, LambdaState, LambdaState}}

    qp_optimizer::Any

    rho_bus::Base.RefValue{Float64}
    rho_edge::Base.RefValue{Float64}
    rho_inj::Base.RefValue{Float64}
    mu::Base.RefValue{Float64}

    alpha1::Float64 = 1.0
    alpha2::Float64 = 1.0
    alpha3::Float64 = 1.0

    sigma_theta_shared::Float64 = 10.0
    sigma_theta_other::Float64 = 10.0
    sigma_pg::Float64 = 1.0
    sigma_other::Float64 = 1.0

    # Per-component diagnostics (not part of the stopping rule)
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

# ------------------------------------------------------------
# Rho helpers (mirror ADMM's _get_rhos for history/logging)
# ------------------------------------------------------------

_get_rhos(ws::DCStateWorkspaceALADIN) = (rho = ws.rho[],)
_get_rhos(ws::DCPowerWorkspaceALADIN) = (rho = ws.rho[],)
_get_rhos(ws::DCHybridWorkspaceALADIN) = (
    rho_bus = ws.rho_bus[],
    rho_edge = ws.rho_edge[],
    rho_inj = ws.rho_inj[],
)

function _build_shared_bus_pairs(
    partitions::Partitions{PartitionDC},
    zones::Vector{AladinZoneData},
)::Vector{SharedBusPair}
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
        length(entries) == 2 || continue
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
    pairs = SharedEdgePair[]

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
        push!(pairs, SharedEdgePair(cidx, e, zi, li, zj, lj))
    end
    return pairs
end

function _build_split_bus_pairs(
    partitions::Partitions{PartitionDC},
    zones::Vector{AladinZoneData},
    sys_bus_pairs::Vector{NTuple{3, Int}},
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}},
)::Vector{SplitBusPair}
    pairs = SplitBusPair[]
    isempty(sys_bus_pairs) && return pairs
    sizehint!(pairs, length(sys_bus_pairs))

    for (k, (bid, za, zb)) in enumerate(sys_bus_pairs)
        v_a = get(p_inj_vars[za], bid, nothing)
        v_b = get(p_inj_vars[zb], bid, nothing)
        (v_a === nothing || v_b === nothing) && continue
        col_a = zones[za].var_to_col[v_a]
        col_b = zones[zb].var_to_col[v_b]
        push!(pairs, SplitBusPair(k, bid, za, col_a, zb, col_b))
    end
    return pairs
end

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

function build_workspace(
    sys::DistributedSystem,
    cfg::ALADINConfig,
    optimizer_factory,
    ::ALADIN{State, DC},
)
    partitions = build_partitions(PartitionDC, sys)
    id = partitions.zones

    build_local_model!(partitions, sys, State, optimizer_factory)

    zones = _build_zone_data_vector(partitions; include_edges = false)
    bus_pairs = _build_shared_bus_pairs(partitions, zones)
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
        sigma_theta_shared = cfg.sigma_theta_shared,
        sigma_theta_other = cfg.sigma_theta_other,
        sigma_pg = cfg.sigma_pg,
        sigma_other = cfg.sigma_other,
    )
end

function build_workspace(
    sys::DistributedSystem,
    cfg::ALADINConfig,
    optimizer_factory,
    ::ALADIN{Power, DC},
)
    partitions = build_partitions(PartitionDC, sys)
    id = partitions.zones
    inter = sys.interzonal_branches

    build_local_model!(partitions, sys, Power, optimizer_factory, inter)

    zones = _build_zone_data_vector(partitions; include_edges = true)
    edge_pairs = _build_shared_edge_pairs(partitions, inter)
    n_pairs = length(edge_pairs)

    return DCPowerWorkspaceALADIN(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        zones = zones,
        edge_pairs = edge_pairs,
        vars = (p_flow = LambdaState(n_pairs),),
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

function build_workspace(
    sys::DistributedSystem,
    cfg::ALADINConfig,
    optimizer_factory,
    ::ALADIN{Hybrid, DC},
)
    partitions = build_partitions(PartitionDC, sys)
    id = partitions.zones
    inter = sys.interzonal_branches
    sys_bus_pairs = sys.bus_pairs

    build_local_model!(partitions, sys, Hybrid, optimizer_factory, inter)

    # SplitBus: attach p_inj variables and wire them into active power balance
    # before we snapshot the variable layout. Reuses ADMM's helper — the
    # partition object is identical.
    n_zones = length(partitions.parts)
    p_inj_vars = [Dict{Int, JuMP.VariableRef}() for _ in 1:n_zones]
    if !isempty(sys_bus_pairs)
        @inbounds for i in 1:n_zones
            _attach_p_inj_variables!(partitions.parts[i], sys_bus_pairs, p_inj_vars[i])
        end
    end

    zones = _build_zone_data_vector(partitions; include_edges = true)
    bus_pairs = _build_shared_bus_pairs(partitions, zones)
    edge_pairs = _build_shared_edge_pairs(partitions, inter)
    split_pairs = _build_split_bus_pairs(partitions, zones, sys_bus_pairs, p_inj_vars)
    n_bus = length(bus_pairs)
    n_edge = length(edge_pairs)
    n_split = length(split_pairs)

    return DCHybridWorkspaceALADIN(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        zones = zones,
        bus_pairs = bus_pairs,
        edge_pairs = edge_pairs,
        split_bus_pairs = split_pairs,
        p_inj_vars = p_inj_vars,
        vars = (
            theta = LambdaState(n_bus),
            p_flow = LambdaState(n_edge),
            p_inj = LambdaState(n_split),
        ),
        qp_optimizer = _resolve_qp_optimizer(optimizer_factory),
        rho_bus = Ref(cfg.rho),
        rho_edge = Ref(cfg.rho),
        rho_inj = Ref(cfg.rho),
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
