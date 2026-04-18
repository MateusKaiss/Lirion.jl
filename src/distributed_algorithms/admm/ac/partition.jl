# ============================================================
# AC PARTITION INITIALIZATION — ADMM
# ------------------------------------------------------------
# See admm/dc/partition.jl header.
# ============================================================

function init_state_partition!(
    part::PartitionAC,
    sys_zone::SystemData,
    optimizer_factory,
    rho_theta::Float64,
    rho_v_mag::Float64,
)
    build_local_model!(part, sys_zone, State, optimizer_factory)
    build_quad_bus_expr!(part, rho_theta, rho_v_mag)
    @objective(part.model, Min, part.obj_base + part.quad_bus_expr)
    return nothing
end

function init_power_partition!(
    part::PartitionAC,
    sys_zone::SystemData,
    optimizer_factory,
    rho_p::Float64,
    rho_q::Float64,
    interzonal_branches,
    branch_pairs::Vector{NTuple{4, Int}} = NTuple{4, Int}[],
)
    build_local_model!(part, sys_zone, Power, optimizer_factory, interzonal_branches)
    build_quad_branch_expr!(part, rho_p, rho_q, part.zone_edges, branch_pairs)
    @objective(part.model, Min, part.obj_base + part.quad_branch_expr)
    return nothing
end

function init_hybrid_partition!(
    part::PartitionAC,
    sys_zone::SystemData,
    optimizer_factory,
    rho_theta::Float64,
    rho_v_mag::Float64,
    rho_p::Float64,
    rho_q::Float64,
    interzonal_branches,
    branch_pairs::Vector{NTuple{4, Int}} = NTuple{4, Int}[],
    bus_pairs::Vector{NTuple{3, Int}} = NTuple{3, Int}[],
    p_inj_store::Union{Nothing, Dict{Int, JuMP.VariableRef}} = nothing,
    q_inj_store::Union{Nothing, Dict{Int, JuMP.VariableRef}} = nothing,
)
    build_local_model!(part, sys_zone, Hybrid, optimizer_factory, interzonal_branches)

    # SplitBus (Method 3): for every split bus in this zone, introduce free
    # p_inj[m] and q_inj[m] scalars and fold them into the AC active /
    # reactive balance at bus m. The AC balance is NONLINEAR (cos/sin/v²),
    # so `set_normalized_coefficient` is not applicable — we delete the
    # existing balance constraints at split buses and rebuild them with the
    # virtual injection added.
    if !isempty(bus_pairs) && p_inj_store !== nothing && q_inj_store !== nothing
        _attach_pq_inj_variables_ac!(
            part, sys_zone, bus_pairs, p_inj_store, q_inj_store,
        )
    end

    build_quad_bus_expr!(part, rho_theta, rho_v_mag)
    build_quad_branch_expr!(part, rho_p, rho_q, part.zone_edges, branch_pairs)
    build_quad_p_inj_expr_ac!(part, rho_p, p_inj_store)
    build_quad_q_inj_expr_ac!(part, rho_q, q_inj_store)

    # Initial @objective: only include the extra terms if this zone actually
    # hosts split buses (avoids expression bloat on pure TieLine runs).
    if p_inj_store !== nothing && !isempty(p_inj_store)
        @objective(
            part.model, Min,
            part.obj_base + part.quad_bus_expr + part.quad_branch_expr +
            part.quad_p_inj_expr + part.quad_q_inj_expr,
        )
    else
        @objective(
            part.model, Min,
            part.obj_base + part.quad_bus_expr + part.quad_branch_expr,
        )
    end
    return nothing
end

# ------------------------------------------------------------
# Delete + rebuild the AC active/reactive balance at each split bus with
# virtual p_inj / q_inj variables folded in.
# ------------------------------------------------------------
function _attach_pq_inj_variables_ac!(
    part::PartitionAC,
    sys_zone::SystemData,
    bus_pairs::Vector{NTuple{3, Int}},
    p_inj_store::Dict{Int, JuMP.VariableRef},
    q_inj_store::Dict{Int, JuMP.VariableRef},
)
    model = part.model::JuMP.Model
    pos = part.position
    pf = sys_zone.powerflow
    bus_idx_map = pf.bus_idx_map
    gs = pf.shunt_conductance
    bs = pf.shunt_susceptance

    zone_split_bus_ids = Int[]
    for (fbus, za, zb) in bus_pairs
        if za == pos || zb == pos
            push!(zone_split_bus_ids, fbus)
        end
    end
    isempty(zone_split_bus_ids) && return nothing

    bal_p = model[:active_power_balance]
    bal_q = model[:reactive_power_balance]
    v_mag = model[:v_magnitude]

    for b in zone_split_bus_ids
        vp = @variable(model, base_name = "p_inj_$(b)")
        JuMP.set_start_value(vp, 0.0)
        p_inj_store[b] = vp

        vq = @variable(model, base_name = "q_inj_$(b)")
        JuMP.set_start_value(vq, 0.0)
        q_inj_store[b] = vq

        JuMP.delete(model, bal_p[b, 1])
        @constraint(
            model,
            compute_active_power_balance(model, sys_zone, b, 1) -
            gs[bus_idx_map[b]] * v_mag[b, 1]^2 + vp == 0
        )

        JuMP.delete(model, bal_q[b, 1])
        @constraint(
            model,
            compute_reactive_power_balance(model, sys_zone, b, 1) +
            bs[bus_idx_map[b]] * v_mag[b, 1]^2 + vq == 0
        )
    end
    return nothing
end

function build_quad_p_inj_expr_ac!(
    part::PartitionAC,
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

function build_quad_q_inj_expr_ac!(
    part::PartitionAC,
    rho::Float64,
    q_inj_store::Union{Nothing, Dict{Int, JuMP.VariableRef}},
)
    quad = zero(JuMP.QuadExpr)
    if q_inj_store !== nothing
        for (_, v) in q_inj_store
            quad += 0.5 * rho * v^2
        end
    end
    part.quad_q_inj_expr = quad
    return quad
end

# ============================================================
# BATCH
# ============================================================

function init_state_partition!(
    partitions::Partitions{PartitionAC},
    sys::DistributedSystem,
    optimizer_factory,
    rho_theta::Float64,
    rho_v_mag::Float64,
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n
    @inbounds for i in 1:n
        init_state_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory,
            rho_theta, rho_v_mag,
        )
    end
    return nothing
end

function init_power_partition!(
    partitions::Partitions{PartitionAC},
    sys::DistributedSystem,
    optimizer_factory,
    rho_p::Float64,
    rho_q::Float64,
    interzonal_branches,
    branch_pairs::Vector{NTuple{4, Int}} = NTuple{4, Int}[],
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n
    @inbounds for i in 1:n
        init_power_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory,
            rho_p, rho_q, interzonal_branches, branch_pairs,
        )
    end
    return nothing
end

function init_hybrid_partition!(
    partitions::Partitions{PartitionAC},
    sys::DistributedSystem,
    optimizer_factory,
    rho_theta::Float64,
    rho_v_mag::Float64,
    rho_p::Float64,
    rho_q::Float64,
    interzonal_branches,
    branch_pairs::Vector{NTuple{4, Int}} = NTuple{4, Int}[],
    bus_pairs::Vector{NTuple{3, Int}} = NTuple{3, Int}[],
    p_inj_vars::Vector{Dict{Int, JuMP.VariableRef}} =
    Dict{Int, JuMP.VariableRef}[],
    q_inj_vars::Vector{Dict{Int, JuMP.VariableRef}} =
    Dict{Int, JuMP.VariableRef}[],
)
    n = length(partitions.parts)
    @assert length(sys.parts) == n
    have_p_inj = length(p_inj_vars) == n
    have_q_inj = length(q_inj_vars) == n
    @inbounds for i in 1:n
        init_hybrid_partition!(
            partitions.parts[i], sys.parts[i], optimizer_factory,
            rho_theta, rho_v_mag, rho_p, rho_q, interzonal_branches,
            branch_pairs,
            bus_pairs,
            have_p_inj ? p_inj_vars[i] : nothing,
            have_q_inj ? q_inj_vars[i] : nothing,
        )
    end
    return nothing
end
