struct SharedBusPairAC
    consensus_idx::Int
    bus_id::Int
    zone_i::Int
    theta_col_i::Int
    vm_col_i::Int
    zone_j::Int
    theta_col_j::Int
    vm_col_j::Int
end

struct SharedEdgePairAC
    consensus_idx::Int
    edge_global::Int
    zone_i::Int
    edge_local_i::Int
    zone_j::Int
    edge_local_j::Int
end

struct SplitBusPairAC
    consensus_idx::Int
    bus_id::Int
    zone_i::Int
    pinj_col_i::Int
    qinj_col_i::Int
    zone_j::Int
    pinj_col_j::Int
    qinj_col_j::Int
end

struct DummyBranchPairAC
    consensus_idx::Int
    zone_i::Int
    local_i::Int
    from_theta_col_i::Int
    to_theta_col_i::Int
    from_vm_col_i::Int
    to_vm_col_i::Int
    g_i::Float64
    b_i::Float64
    b_sh_i::Float64
    inv_tap_i::Float64
    tap_angle_i::Float64
    zone_j::Int
    local_j::Int
    from_theta_col_j::Int
    to_theta_col_j::Int
    from_vm_col_j::Int
    to_vm_col_j::Int
    g_j::Float64
    b_j::Float64
    b_sh_j::Float64
    inv_tap_j::Float64
    tap_angle_j::Float64
end

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

_get_rhos(ws::ACStateWorkspaceALADIN) = (rho = ws.rho[],)
_get_rhos(ws::ACPowerWorkspaceALADIN) = (rho = ws.rho[],)
_get_rhos(ws::ACHybridWorkspaceALADIN) = (rho = ws.rho[],)

function _build_shared_bus_pairs_ac(
    partitions::Partitions{PartitionAC},
    zones::Vector{AladinZoneDataAC},
)::Vector{SharedBusPairAC}
    acc = Dict{Int, Vector{Tuple{Int, Int, Int}}}()
    bus_of = Dict{Int, Int}()

    for (zpos, part) in enumerate(partitions.parts)
        zdata = zones[zpos]
        @inbounds for i in eachindex(part.zone_shared_bus_ids)
            cidx = part.shared_bus_consensus_idx[i]
            cidx == 0 && continue
            tcol = zdata.shared_theta_cols[i]
            vcol = zdata.shared_vm_cols[i]
            v = get!(acc, cidx, Tuple{Int, Int, Int}[])
            push!(v, (zpos, tcol, vcol))
            bus_of[cidx] = part.zone_shared_bus_ids[i]
        end
    end

    pairs = SharedBusPairAC[]
    for cidx in sort!(collect(keys(acc)))
        entries = acc[cidx]
        length(entries) == 2 || continue
        zi, ti, vi = entries[1]
        zj, tj, vj = entries[2]
        push!(pairs, SharedBusPairAC(cidx, bus_of[cidx], zi, ti, vi, zj, tj, vj))
    end
    return pairs
end

function _build_shared_edge_pairs_ac(
    partitions::Partitions{PartitionAC},
    interzonal::Vector{Int},
)::Vector{SharedEdgePairAC}
    pairs = SharedEdgePairAC[]
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
        push!(pairs, SharedEdgePairAC(cidx, e, zi, li, zj, lj))
    end
    return pairs
end

function _build_dummy_branch_pairs_ac(
    partitions::Partitions{PartitionAC},
    zones::Vector{AladinZoneDataAC},
    sys_branch_pairs::Vector{NTuple{4, Int}},
)::Vector{DummyBranchPairAC}
    pairs = DummyBranchPairAC[]
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

        theta_a = pa.model[:theta]
        theta_b = pb.model[:theta]
        vm_a = pa.model[:v_magnitude]
        vm_b = pb.model[:v_magnitude]

        fa = pa.branch.from[la]
        ta = pa.branch.to[la]
        fb = pb.branch.from[lb]
        tb = pb.branch.to[lb]

        from_theta_col_i = zda.var_to_col[theta_a[fa, 1]]
        to_theta_col_i = zda.var_to_col[theta_a[ta, 1]]
        from_vm_col_i = zda.var_to_col[vm_a[fa, 1]]
        to_vm_col_i = zda.var_to_col[vm_a[ta, 1]]

        from_theta_col_j = zdb.var_to_col[theta_b[fb, 1]]
        to_theta_col_j = zdb.var_to_col[theta_b[tb, 1]]
        from_vm_col_j = zdb.var_to_col[vm_b[fb, 1]]
        to_vm_col_j = zdb.var_to_col[vm_b[tb, 1]]

        push!(
            pairs,
            DummyBranchPairAC(
                k,
                za, la,
                from_theta_col_i, to_theta_col_i, from_vm_col_i, to_vm_col_i,
                pa.branch.g[la], pa.branch.b[la], pa.branch.b_sh[la],
                pa.branch.inv_tap[la], pa.branch.tap_angle[la],
                zb, lb,
                from_theta_col_j, to_theta_col_j, from_vm_col_j, to_vm_col_j,
                pb.branch.g[lb], pb.branch.b[lb], pb.branch.b_sh[lb],
                pb.branch.inv_tap[lb], pb.branch.tap_angle[lb],
            ),
        )
    end
    return pairs
end

function _build_split_bus_pairs_ac(
    zones::Vector{AladinZoneDataAC},
    sys_bus_pairs::Vector{NTuple{3, Int}},
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}},
    q_inj_vars::Vector{Dict{Int, JuMP.VariableRef}},
)::Vector{SplitBusPairAC}
    pairs = SplitBusPairAC[]
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
        push!(pairs, SplitBusPairAC(
            k, bid, za, pcol_a, qcol_a, zb, pcol_b, qcol_b,
        ))
    end
    return pairs
end

function _build_zone_data_vector_ac(
    partitions::Partitions{PartitionAC};
    include_edges::Bool,
)::Vector{AladinZoneDataAC}
    parts = partitions.parts
    zones = Vector{AladinZoneDataAC}(undef, length(parts))
    @inbounds for i in eachindex(parts)
        zones[i] = build_aladin_zone_data_ac(parts[i]; include_edges = include_edges)
    end
    return zones
end

function _prepare_pq_inj_vars_ac!(
    partitions,
    sys,
    sys_bus_pairs,
    ::Union{TieLine, DummyBus},
)
    n_zones = length(partitions.parts)
    p_inj_vars = [Dict{Int, JuMP.VariableRef}() for _ in 1:n_zones]
    q_inj_vars = [Dict{Int, JuMP.VariableRef}() for _ in 1:n_zones]
    return p_inj_vars, q_inj_vars
end

function _prepare_pq_inj_vars_ac!(partitions, sys, sys_bus_pairs, ::SplitBus)
    n_zones = length(partitions.parts)
    p_inj_vars = [Dict{Int, JuMP.VariableRef}() for _ in 1:n_zones]
    q_inj_vars = [Dict{Int, JuMP.VariableRef}() for _ in 1:n_zones]
    @inbounds for i in 1:n_zones
        _attach_pq_inj_variables_ac!(
            partitions.parts[i], sys.parts[i], sys_bus_pairs,
            p_inj_vars[i], q_inj_vars[i],
        )
    end
    return p_inj_vars, q_inj_vars
end

function build_workspace(
    sys::DistributedSystem,
    cfg::ALADINConfig,
    optimizer_factory,
    ::ALADIN{State, AC},
)
    partitions = build_partitions(PartitionAC, sys)
    id = partitions.zones
    inter = sys.interzonal_branches

    build_local_model!(partitions, sys, State, optimizer_factory)

    zones = _build_zone_data_vector_ac(partitions; include_edges = false)
    bus_pairs = _build_shared_bus_pairs_ac(partitions, zones)
    n_bus = length(bus_pairs)

    return ACStateWorkspaceALADIN(;
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
    )
end

function build_workspace(
    sys::DistributedSystem,
    cfg::ALADINConfig,
    optimizer_factory,
    ::ALADIN{Power, AC},
)
    partitions = build_partitions(PartitionAC, sys)
    id = partitions.zones
    inter = sys.interzonal_branches
    sys_branch_pairs = sys.branch_pairs

    build_local_model!(partitions, sys, Power, optimizer_factory, inter)

    zones = _build_zone_data_vector_ac(partitions; include_edges = true)
    edge_pairs = _build_shared_edge_pairs_ac(partitions, inter)
    dummy_pairs = _build_dummy_branch_pairs_ac(partitions, zones, sys_branch_pairs)
    n_edge = length(edge_pairs)
    n_dummy = length(dummy_pairs)

    return ACPowerWorkspaceALADIN(;
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
    )
end

function build_workspace(
    sys::DistributedSystem,
    cfg::ALADINConfig,
    optimizer_factory,
    ::ALADIN{Hybrid, AC},
)
    partitions = build_partitions(PartitionAC, sys)
    id = partitions.zones
    inter = sys.interzonal_branches
    sys_bus_pairs = sys.bus_pairs

    build_local_model!(partitions, sys, Hybrid, optimizer_factory, inter)

    p_inj_vars, q_inj_vars = _prepare_pq_inj_vars_ac!(
        partitions, sys, sys_bus_pairs, cfg.sectioning,
    )

    zones = _build_zone_data_vector_ac(partitions; include_edges = true)
    bus_pairs = _build_shared_bus_pairs_ac(partitions, zones)
    edge_pairs = _build_shared_edge_pairs_ac(partitions, inter)
    split_pairs = _build_split_bus_pairs_ac(
        zones, sys_bus_pairs, p_inj_vars, q_inj_vars,
    )
    dummy_pairs = _build_dummy_branch_pairs_ac(partitions, zones, sys.branch_pairs)
    n_bus = length(bus_pairs)
    n_edge = length(edge_pairs)
    n_split = length(split_pairs)
    n_dummy = length(dummy_pairs)

    return ACHybridWorkspaceALADIN(;
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
    )
end
