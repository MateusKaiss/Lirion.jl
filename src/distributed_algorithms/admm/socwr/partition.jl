# ============================================================
# Tikhonov regularization
# ============================================================
#
# `_attach_pq_inj_variables_socwr!` lives in
# `distributed_algorithms/shared/sectioning_helpers.jl` (shared with
# ALADIN — pure JuMP plumbing, algorithm-agnostic).

function tikhonov_regularization_expr!(
    part::PartitionSOCWR;
    kappa_v2::Float64 = 1e-6,
    kappa_pq::Float64 = 1e-8,
    use_bus::Bool = true,
    use_pq::Bool = true,
)
    m = part.local_model::JuMP.Model
    quad = zero(JuMP.QuadExpr)

    if use_bus && haskey(m, :v_magnitude_sq)
        v2 = m[:v_magnitude_sq]
        part.v2_var = v2

        for I in keys(v2)
            x = v2[I]
            x_ref = (JuMP.start_value(x) === nothing) ? 1.0 : JuMP.start_value(x)
            quad += 0.5 * kappa_v2 * (x - x_ref)^2
        end
    end

    if use_pq && !isempty(part.zone_edges)
        @inbounds for e in part.zone_edges
            l = part.branch_idx_full_to_zone[e]
            l == 0 && continue

            p, q = _pq_flow_from_exprs(m, l)
            quad += 0.5 * kappa_pq * (p^2 + q^2)
        end
    end

    part.tikhonov_reg_expr = quad
    return quad
end

# ============================================================
# Per-partition initializers
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function init_state_partition!(
    part::PartitionSOCWR,
    sys_zone::SystemData,
    optimizer_factory,
    interzonal_branches,
)
    build_local_model!(part, sys_zone, State, optimizer_factory, interzonal_branches)
    return nothing
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function init_power_partition!(
    part::PartitionSOCWR,
    sys_zone::SystemData,
    optimizer_factory,
    interzonal_branches,
)
    build_local_model!(part, sys_zone, Power, optimizer_factory, interzonal_branches)
    if is_gurobi(part.local_model)
        tikhonov_regularization_expr!(
            part;
            kappa_v2 = 1e-6,
            kappa_pq = 1e-8,
            use_bus = true,
            use_pq = true,
        )
    else
        part.tikhonov_reg_expr = zero(JuMP.QuadExpr)
    end
    return nothing
end

# ------------------------------------------------------------
# Hybrid — TieLine
# ------------------------------------------------------------

function init_hybrid_partition!(
    part::PartitionSOCWR,
    sys_zone::SystemData,
    optimizer_factory,
    ::TieLine,
    interzonal_branches,
)
    build_local_model!(part, sys_zone, Hybrid, optimizer_factory, interzonal_branches)
    return nothing
end

# ------------------------------------------------------------
# Hybrid — DummyBus
# ------------------------------------------------------------

function init_hybrid_partition!(
    part::PartitionSOCWR,
    sys_zone::SystemData,
    optimizer_factory,
    ::DummyBus,
    interzonal_branches,
)
    build_local_model!(part, sys_zone, Hybrid, optimizer_factory, interzonal_branches)
    return nothing
end

# ------------------------------------------------------------
# Hybrid — SplitBus
# ------------------------------------------------------------

function init_hybrid_partition!(
    part::PartitionSOCWR,
    sys_zone::SystemData,
    optimizer_factory,
    ::SplitBus,
    interzonal_branches,
    bus_pairs::Vector{NTuple{3, Int}},
    p_inj_store::Dict{Int, JuMP.VariableRef},
    q_inj_store::Dict{Int, JuMP.VariableRef},
)
    build_local_model!(part, sys_zone, Hybrid, optimizer_factory, interzonal_branches)
    _attach_pq_inj_variables_socwr!(
        part, bus_pairs, p_inj_store, q_inj_store,
    )
    return nothing
end

# ============================================================
# Batch entry points
# ============================================================

# ------------------------------------------------------------
# State
# ------------------------------------------------------------

function init_state_partition!(
    partitions::Partitions{PartitionSOCWR},
    sys::DistributedSystem,
    optimizer_factory,
    interzonal_branches,
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n
    @inbounds for i in 1:n
        init_state_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory,
            interzonal_branches,
        )
    end
    return nothing
end

# ------------------------------------------------------------
# Power
# ------------------------------------------------------------

function init_power_partition!(
    partitions::Partitions{PartitionSOCWR},
    sys::DistributedSystem,
    optimizer_factory,
    interzonal_branches,
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n
    @inbounds for i in 1:n
        init_power_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory,
            interzonal_branches,
        )
    end
    return nothing
end

# ------------------------------------------------------------
# Hybrid — TieLine
# ------------------------------------------------------------

function init_hybrid_partition!(
    partitions::Partitions{PartitionSOCWR},
    sys::DistributedSystem,
    optimizer_factory,
    sectioning::TieLine,
    interzonal_branches,
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n
    @inbounds for i in 1:n
        init_hybrid_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory,
            sectioning, interzonal_branches,
        )
    end
    return nothing
end

# ------------------------------------------------------------
# Hybrid — DummyBus
# ------------------------------------------------------------

function init_hybrid_partition!(
    partitions::Partitions{PartitionSOCWR},
    sys::DistributedSystem,
    optimizer_factory,
    sectioning::DummyBus,
    interzonal_branches,
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n
    @inbounds for i in 1:n
        init_hybrid_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory,
            sectioning, interzonal_branches,
        )
    end
    return nothing
end

# ------------------------------------------------------------
# Hybrid — SplitBus
# ------------------------------------------------------------

function init_hybrid_partition!(
    partitions::Partitions{PartitionSOCWR},
    sys::DistributedSystem,
    optimizer_factory,
    sectioning::SplitBus,
    interzonal_branches,
    bus_pairs::Vector{NTuple{3, Int}},
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}},
    q_inj_vars::Vector{Dict{Int, JuMP.VariableRef}},
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n
    @assert length(p_inj_vars) == n
    @assert length(q_inj_vars) == n
    @inbounds for i in 1:n
        init_hybrid_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory,
            sectioning, interzonal_branches,
            bus_pairs, p_inj_vars[i], q_inj_vars[i],
        )
    end
    return nothing
end
