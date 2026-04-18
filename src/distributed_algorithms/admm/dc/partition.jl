# ============================================================
# DC PARTITION INITIALIZATION — ADMM
# ------------------------------------------------------------

function init_state_partition!(
    part::PartitionDC,
    sys_zone::SystemData,
    optimizer_factory,
    rho_theta::Float64,
)
    build_local_model!(part, sys_zone, State, optimizer_factory)
    build_quad_bus_expr!(part, rho_theta)
    part.obj_static = part.obj_base + part.quad_bus_expr
    @objective(part.model, Min, part.obj_static)
    return nothing
end

function init_power_partition!(
    part::PartitionDC,
    sys_zone::SystemData,
    optimizer_factory,
    rho::Float64,
    interzonal_branches::Vector{Int},
    branch_pairs::Vector{NTuple{4, Int}} = NTuple{4, Int}[],
)
    build_local_model!(part, sys_zone, Power, optimizer_factory, interzonal_branches)
    build_quad_branch_expr!(part, rho, part.zone_edges, branch_pairs)
    part.obj_static = part.obj_base + part.quad_branch_expr
    @objective(part.model, Min, part.obj_static)
    return nothing
end

function init_hybrid_partition!(
    part::PartitionDC,
    sys_zone::SystemData,
    optimizer_factory,
    rho_theta::Float64,
    rho::Float64,
    interzonal_branches::Vector{Int},
    branch_pairs::Vector{NTuple{4, Int}} = NTuple{4, Int}[],
    bus_pairs::Vector{NTuple{3, Int}} = NTuple{3, Int}[],
    p_inj_store::Union{Nothing, Dict{Int, JuMP.VariableRef}} = nothing,
)
    build_local_model!(part, sys_zone, Hybrid, optimizer_factory, interzonal_branches)

    # SplitBus (Method 3): add one free variable p_inj[m] per split bus present
    # in this zone and add it into the active power balance at that bus via
    # `set_normalized_coefficient`. The bus itself is already marked
    # non-virtual in the partitioning step, so the balance constraint exists.
    if !isempty(bus_pairs) && p_inj_store !== nothing
        _attach_p_inj_variables!(part, bus_pairs, p_inj_store)
    end

    build_quad_bus_expr!(part, rho_theta)
    build_quad_branch_expr!(part, rho, part.zone_edges, branch_pairs)
    build_quad_p_inj_expr!(part, rho, p_inj_store)

    part.obj_static =
        part.obj_base + part.quad_bus_expr + part.quad_branch_expr +
        part.quad_p_inj_expr
    @objective(part.model, Min, part.obj_static)
    return nothing
end

# Build (or rebuild) the quadratic term 0.5·ρ · Σ p_inj[b]² for this zone.
# Stored on `part.quad_p_inj_expr` so that rho refresh can regenerate obj_static
# consistently with the existing `quad_bus_expr`/`quad_branch_expr` pattern.
function build_quad_p_inj_expr!(
    part::PartitionDC,
    rho::Float64,
    p_inj_store::Union{Nothing, Dict{Int, JuMP.VariableRef}},
)
    quad = zero(JuMP.QuadExpr)
    if p_inj_store !== nothing
        for (_, v) in p_inj_store
            quad += 0.5 * rho * v^2
        end
    end
    part.quad_p_inj_expr = quad
    return quad
end

# ------------------------------------------------------------
# Helper: create p_inj variables in this partition and inject them into the
# existing active power balance constraints at split buses.
# ------------------------------------------------------------
function _attach_p_inj_variables!(
    part::PartitionDC,
    bus_pairs::Vector{NTuple{3, Int}},
    p_inj_store::Dict{Int, JuMP.VariableRef},
)
    model = part.model::JuMP.Model
    pos = part.position

    # Buses in *this* zone that are split points (appear in bus_pairs either
    # as zone_a_pos or zone_b_pos equal to this partition's position).
    zone_split_bus_ids = Int[]
    for (fbus, za, zb) in bus_pairs
        if za == pos || zb == pos
            push!(zone_split_bus_ids, fbus)
        end
    end
    isempty(zone_split_bus_ids) && return nothing

    bal_con = model[:active_power_balance]

    for b in zone_split_bus_ids
        v = @variable(model, base_name = "p_inj_$(b)")
        JuMP.set_start_value(v, 0.0)
        p_inj_store[b] = v
        JuMP.set_normalized_coefficient(bal_con[b, 1], v, 1.0)
    end
    return nothing
end

# ============================================================
# BATCH
# ============================================================

function init_state_partition!(
    partitions::Partitions{PartitionDC},
    sys::DistributedSystem,
    optimizer_factory,
    rho_theta::Float64,
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n
    @inbounds for i in 1:n
        init_state_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory, rho_theta,
        )
    end
    return nothing
end

function init_power_partition!(
    partitions::Partitions{PartitionDC},
    sys::DistributedSystem,
    optimizer_factory,
    rho::Float64,
    interzonal_branches::Vector{Int},
    branch_pairs::Vector{NTuple{4, Int}} = NTuple{4, Int}[],
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n
    @inbounds for i in 1:n
        init_power_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory,
            rho, interzonal_branches, branch_pairs,
        )
    end
    return nothing
end

function init_hybrid_partition!(
    partitions::Partitions{PartitionDC},
    sys::DistributedSystem,
    optimizer_factory,
    rho_theta::Float64,
    rho::Float64,
    interzonal_branches::Vector{Int},
    branch_pairs::Vector{NTuple{4, Int}} = NTuple{4, Int}[],
    bus_pairs::Vector{NTuple{3, Int}} = NTuple{3, Int}[],
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}} =
    Dict{Int, JuMP.VariableRef}[],
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n
    have_p_inj = length(p_inj_vars) == n
    @inbounds for i in 1:n
        init_hybrid_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory,
            rho_theta, rho, interzonal_branches, branch_pairs,
            bus_pairs,
            have_p_inj ? p_inj_vars[i] : nothing,
        )
    end
    return nothing
end
