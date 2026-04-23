# ============================================================
# Pair descriptors
# ============================================================

struct SharedBusPair
    consensus_idx::Int
    bus_id::Int
    zone_i::Int
    col_i::Int
    zone_j::Int
    col_j::Int
end

struct SharedEdgePair
    consensus_idx::Int
    edge_global::Int
    zone_i::Int
    edge_local_i::Int
    zone_j::Int
    edge_local_j::Int
end

struct SplitBusPair
    consensus_idx::Int
    bus_id::Int
    zone_i::Int
    pinj_col_i::Int
    zone_j::Int
    pinj_col_j::Int
end

struct DummyBranchPairDC
    consensus_idx::Int
    zone_i::Int
    from_col_i::Int
    to_col_i::Int
    b_i::Float64
    zone_j::Int
    from_col_j::Int
    to_col_j::Int
    b_j::Float64
end


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
# Power
# ------------------------------------------------------------

Base.@kwdef mutable struct DCPowerWorkspaceALADIN <: AbstractALADINWorkspace
    partitions::Partitions{PartitionDC}
    id::Vector{Int}
    interzonal::Vector{Int}

    zones::Vector{AladinZoneData}
    edge_pairs::Vector{SharedEdgePair}
    dummy_branch_pairs::Vector{DummyBranchPairDC} = DummyBranchPairDC[]

    vars::NamedTuple{(:p_flow, :p_dummy), Tuple{LambdaState, LambdaState}}

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
# Hybrid
# ------------------------------------------------------------

Base.@kwdef mutable struct DCHybridWorkspaceALADIN <: AbstractALADINWorkspace
    partitions::Partitions{PartitionDC}
    id::Vector{Int}
    interzonal::Vector{Int}

    zones::Vector{AladinZoneData}
    bus_pairs::Vector{SharedBusPair}
    edge_pairs::Vector{SharedEdgePair}
    split_bus_pairs::Vector{SplitBusPair} = SplitBusPair[]
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}} =
        Vector{Dict{Int, JuMP.VariableRef}}()
    dummy_branch_pairs::Vector{DummyBranchPairDC} = DummyBranchPairDC[]

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

    sigma_theta_shared::Float64 = 10.0
    sigma_theta_other::Float64 = 10.0
    sigma_pg::Float64 = 1.0
    sigma_other::Float64 = 1.0

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
    rho_bus  = ws.rho_bus[],
    rho_edge = ws.rho_edge[],
    rho_inj  = ws.rho_inj[],
)


# ============================================================
# Helpers
# ============================================================

function _build_shared_bus_pairs(
    partitions::Partitions{PartitionDC},
    zones::Vector{AladinZoneData},
)::Vector{SharedBusPair}
    acc = Dict{Int, Vector{Tuple{Int, Int}}}()
    bus_of = Dict{Int, Int}()

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

function _build_dummy_branch_pairs_dc(
    partitions::Partitions{PartitionDC},
    zones::Vector{AladinZoneData},
    sys_branch_pairs::Vector{NTuple{4, Int}},
)::Vector{DummyBranchPairDC}
    pairs = DummyBranchPairDC[]
    isempty(sys_branch_pairs) && return pairs
    sizehint!(pairs, length(sys_branch_pairs))

    for (k, (ga, za, gb, zb)) in enumerate(sys_branch_pairs)
        pa = partitions.parts[za]
        pb = partitions.parts[zb]
        la = pa.branch_idx_full_to_zone[ga]
        lb = pb.branch_idx_full_to_zone[gb]
        (la == 0 || lb == 0) && continue

        zda = zones[za]
        zdb = zones[zb]

        fa = pa.branch.from[la]
        ta = pa.branch.to[la]
        fb = pb.branch.from[lb]
        tb = pb.branch.to[lb]

        theta_a = pa.model[:theta]
        theta_b = pb.model[:theta]
        from_col_i = zda.var_to_col[theta_a[fa, 1]]
        to_col_i = zda.var_to_col[theta_a[ta, 1]]
        from_col_j = zdb.var_to_col[theta_b[fb, 1]]
        to_col_j = zdb.var_to_col[theta_b[tb, 1]]

        b_i = pa.branch.b[la]
        b_j = pb.branch.b[lb]

        push!(
            pairs,
            DummyBranchPairDC(
                k, za, from_col_i, to_col_i, b_i,
                zb, from_col_j, to_col_j, b_j,
            ),
        )
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

_prepare_p_inj_vars_dc!(partitions, sys_bus_pairs, ::Union{TieLine, DummyBus}) =
    [Dict{Int, JuMP.VariableRef}() for _ in 1:length(partitions.parts)]

function _prepare_p_inj_vars_dc!(partitions, sys_bus_pairs, ::SplitBus)
    n_zones = length(partitions.parts)
    p_inj_vars = [Dict{Int, JuMP.VariableRef}() for _ in 1:n_zones]
    @inbounds for i in 1:n_zones
        _attach_p_inj_variables!(partitions.parts[i], sys_bus_pairs, p_inj_vars[i])
    end
    return p_inj_vars
end


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
    sys_branch_pairs = sys.branch_pairs

    build_local_model!(partitions, sys, Power, optimizer_factory, inter)

    zones = _build_zone_data_vector(partitions; include_edges = true)
    edge_pairs = _build_shared_edge_pairs(partitions, inter)
    dummy_pairs = _build_dummy_branch_pairs_dc(partitions, zones, sys_branch_pairs)
    n_pairs = length(edge_pairs)
    n_dummy = length(dummy_pairs)

    return DCPowerWorkspaceALADIN(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        zones = zones,
        edge_pairs = edge_pairs,
        dummy_branch_pairs = dummy_pairs,
        vars = (
            p_flow  = LambdaState(n_pairs),
            p_dummy = LambdaState(n_dummy),
        ),
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
    sys_bus_pairs = sys.bus_pairs

    build_local_model!(partitions, sys, Hybrid, optimizer_factory, inter)

    p_inj_vars = _prepare_p_inj_vars_dc!(partitions, sys_bus_pairs, cfg.sectioning)

    zones = _build_zone_data_vector(partitions; include_edges = true)
    bus_pairs = _build_shared_bus_pairs(partitions, zones)
    edge_pairs = _build_shared_edge_pairs(partitions, inter)
    split_pairs = _build_split_bus_pairs(partitions, zones, sys_bus_pairs, p_inj_vars)
    dummy_pairs = _build_dummy_branch_pairs_dc(partitions, zones, sys.branch_pairs)
    n_bus = length(bus_pairs)
    n_edge = length(edge_pairs)
    n_split = length(split_pairs)
    n_dummy = length(dummy_pairs)

    return DCHybridWorkspaceALADIN(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        zones = zones,
        bus_pairs = bus_pairs,
        edge_pairs = edge_pairs,
        split_bus_pairs = split_pairs,
        p_inj_vars = p_inj_vars,
        dummy_branch_pairs = dummy_pairs,
        vars = (
            theta   = LambdaState(n_bus),
            p_flow  = LambdaState(n_edge),
            p_inj   = LambdaState(n_split),
            p_dummy = LambdaState(n_dummy),
        ),
        qp_optimizer = _resolve_qp_optimizer(optimizer_factory),
        rho_bus  = Ref(cfg.rho),
        rho_edge = Ref(cfg.rho),
        rho_inj  = Ref(cfg.rho),
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
