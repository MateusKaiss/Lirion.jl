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
    edge_local_i::Int
    zone_j::Int
    edge_local_j::Int
end

# SplitBus (Method 3) pairs for SOCWR: sign-flipped p_inj / q_inj consensus
# at each split bus. Both zones host p_inj and q_inj slack variables folded
# into the AFFINE active/reactive balance; consensus target is p_A + p_B = 0
# (and similarly for q).
struct SplitBusPairSOCWR
    consensus_idx::Int
    bus_id::Int
    zone_i::Int
    pinj_col_i::Int
    qinj_col_i::Int
    zone_j::Int
    pinj_col_j::Int
    qinj_col_j::Int
end

Base.@kwdef mutable struct SOCWRStateWorkspaceALADIN <: AbstractALADINWorkspace
    partitions::Partitions{PartitionSOCWR}
    id::Vector{Int}
    interzonal::Vector{Int}

    zones::Vector{AladinZoneDataSOCWR}
    bus_pairs::Vector{SharedBusPairSOCWR}
    edge_pairs::Vector{SharedEdgePairSOCWR}

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
# POWER workspace
# ------------------------------------------------------------

Base.@kwdef mutable struct SOCWRPowerWorkspaceALADIN <: AbstractALADINWorkspace
    partitions::Partitions{PartitionSOCWR}
    id::Vector{Int}
    interzonal::Vector{Int}

    zones::Vector{AladinZoneDataSOCWR}

    edge_pairs::Vector{SharedEdgePairSOCWR}

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

    r_primal::Float64 = 0.0
    r_progress::Float64 = 0.0
    lambda_step::Float64 = 0.0
    eps_primal::Float64 = 0.0
    eps_progress::Float64 = 0.0
    iter::Int = 0
    hist::ALADINHistory = ALADINHistory()
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
    # SplitBus (Method 3): sign-flipped p_inj / q_inj consensus per split bus.
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

# ------------------------------------------------------------
# helpers
# ------------------------------------------------------------

_get_rhos(ws::SOCWRStateWorkspaceALADIN) = (rho = ws.rho[],)
_get_rhos(ws::SOCWRPowerWorkspaceALADIN) = (rho = ws.rho[],)
_get_rhos(ws::SOCWRHybridWorkspaceALADIN) = (rho = ws.rho[],)

function _build_shared_bus_pairs_socwr(
    partitions::Partitions{PartitionSOCWR},
    zones::Vector{AladinZoneDataSOCWR},
)::Vector{SharedBusPairSOCWR}
    acc = Dict{Int, Vector{Tuple{Int, Int}}}()
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

function _build_split_bus_pairs_socwr(
    zones::Vector{AladinZoneDataSOCWR},
    sys_bus_pairs::Vector{NTuple{3, Int}},
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}},
    q_inj_vars::Vector{Dict{Int, JuMP.VariableRef}},
)::Vector{SplitBusPairSOCWR}
    pairs = SplitBusPairSOCWR[]
    isempty(sys_bus_pairs) && return pairs
    sizehint!(pairs, length(sys_bus_pairs))

    for (k, (bid, za, zb)) in enumerate(sys_bus_pairs)
        vpa = get(p_inj_vars[za], bid, nothing)
        vqa = get(q_inj_vars[za], bid, nothing)
        vpb = get(p_inj_vars[zb], bid, nothing)
        vqb = get(q_inj_vars[zb], bid, nothing)
        (vpa === nothing || vqa === nothing ||
         vpb === nothing || vqb === nothing) && continue
        pcol_a = zones[za].var_to_col[vpa]
        qcol_a = zones[za].var_to_col[vqa]
        pcol_b = zones[zb].var_to_col[vpb]
        qcol_b = zones[zb].var_to_col[vqb]
        push!(pairs, SplitBusPairSOCWR(
            k, bid, za, pcol_a, qcol_a, zb, pcol_b, qcol_b,
        ))
    end
    return pairs
end

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
# Builders
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

    build_local_model!(partitions, sys, State, optimizer_factory, inter)

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
end

function build_workspace(
    sys::DistributedSystem,
    cfg::ALADINConfig,
    optimizer_factory,
    ::ALADIN{Power, SOCWR},
)
    partitions = build_partitions(PartitionSOCWR, sys)
    id = partitions.zones
    inter = sys.interzonal_branches

    build_local_model!(partitions, sys, Power, optimizer_factory, inter)

    zones = _build_zone_data_vector_socwr(partitions; include_edges = true)
    edge_pairs = _build_shared_edge_pairs_socwr(partitions, inter)
    n_edge = length(edge_pairs)

    return SOCWRPowerWorkspaceALADIN(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        zones = zones,
        edge_pairs = edge_pairs,
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
    )
end

function build_workspace(
    sys::DistributedSystem,
    cfg::ALADINConfig,
    optimizer_factory,
    ::ALADIN{Hybrid, SOCWR},
)
    partitions = build_partitions(PartitionSOCWR, sys)
    id = partitions.zones
    inter = sys.interzonal_branches
    sys_bus_pairs = sys.bus_pairs

    build_local_model!(partitions, sys, Hybrid, optimizer_factory, inter)

    # SplitBus: attach p_inj / q_inj slacks into the AFFINE active/reactive
    # balance at each split bus (reuses the ADMM helper) before snapshotting
    # the variable layout.
    n_zones = length(partitions.parts)
    p_inj_vars = [Dict{Int, JuMP.VariableRef}() for _ in 1:n_zones]
    q_inj_vars = [Dict{Int, JuMP.VariableRef}() for _ in 1:n_zones]
    if !isempty(sys_bus_pairs)
        @inbounds for i in 1:n_zones
            _attach_pq_inj_variables_socwr!(
                partitions.parts[i], sys_bus_pairs,
                p_inj_vars[i], q_inj_vars[i],
            )
        end
    end

    zones = _build_zone_data_vector_socwr(partitions; include_edges = true)
    bus_pairs = _build_shared_bus_pairs_socwr(partitions, zones)
    edge_pairs = _build_shared_edge_pairs_socwr(partitions, inter)
    split_pairs = _build_split_bus_pairs_socwr(
        zones, sys_bus_pairs, p_inj_vars, q_inj_vars,
    )

    n_bus = length(bus_pairs)
    n_edge = length(edge_pairs)
    n_split = length(split_pairs)

    return SOCWRHybridWorkspaceALADIN(;
        partitions = partitions,
        id = id,
        interzonal = inter,
        zones = zones,
        bus_pairs = bus_pairs,
        edge_pairs = edge_pairs,
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
end
