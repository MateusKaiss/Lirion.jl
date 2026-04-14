abstract type AbstractPartition end

struct Partitions{P <: AbstractPartition}
    parts::Vector{P}
    zones::Vector{Int}
end

@inline Base.length(ps::Partitions) = length(ps.parts)

@inline function Base.iterate(ps::Partitions, st::Int = 1)
    st > length(ps.parts) && return nothing
    @inbounds return (ps.parts[st], st + 1)
end

Base.@kwdef mutable struct PartitionDC <: AbstractPartition
    id::Int
    position::Int = 0

    bus::BusData
    branch::BranchData

    zone_branch_idx_in_full::Vector{Int} = Int[]
    zone_shared_bus_ids::Vector{Int} = Int[]
    branch_idx_full_to_zone::Vector{Int} = Int[]

    shared_bus_local_idx::Vector{Int} = Int[]
    shared_bus_consensus_idx::Vector{Int} = Int[]

    theta_buf::Vector{Float64} = Float64[]
    theta_shared_buf::Vector{Float64} = Float64[]

    model::Union{Nothing, JuMP.Model} = nothing
    theta_var::Union{Nothing, JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}} = nothing

    obj_base::Union{Nothing, Float64, JuMP.AffExpr, JuMP.QuadExpr} = nothing
    obj_static::JuMP.QuadExpr = zero(JuMP.QuadExpr)
    quad_bus_expr::JuMP.QuadExpr = zero(JuMP.QuadExpr)
    quad_branch_expr::JuMP.QuadExpr = zero(JuMP.QuadExpr)

    objective_params::ObjectiveData
    has_solution::Bool = false

    generation_cost::Float64 = 0.0

    zone_edges::Vector{Int} = Int[]
end

Base.@kwdef mutable struct PartitionAC <: AbstractPartition
    id::Int
    position::Int = 0

    bus::BusData
    branch::BranchData

    zone_branch_idx_in_full::Vector{Int} = Int[]
    zone_shared_bus_ids::Vector{Int} = Int[]
    branch_idx_full_to_zone::Vector{Int} = Int[]

    shared_bus_local_idx::Vector{Int} = Int[]
    shared_bus_consensus_idx::Vector{Int} = Int[]

    theta_buf::Vector{Float64} = Float64[]
    v_mag_buf::Vector{Float64} = Float64[]
    theta_shared_buf::Vector{Float64} = Float64[]
    v_mag_shared_buf::Vector{Float64} = Float64[]

    model::Union{Nothing, JuMP.Model} = nothing
    theta_var::Union{Nothing, JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}} = nothing
    v_mag_var::Union{Nothing, JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}} = nothing

    obj_base::Union{Nothing, Float64, JuMP.AffExpr, JuMP.QuadExpr, JuMP.NonlinearExpr} =
        nothing
    quad_bus_expr::Union{Nothing, JuMP.QuadExpr, JuMP.NonlinearExpr} = nothing
    quad_branch_expr::Union{Nothing, JuMP.QuadExpr, JuMP.NonlinearExpr} = nothing

    objective_params::ObjectiveData
    has_solution::Bool = false

    generation_cost::Float64 = 0.0

    zone_edges::Vector{Int} = Int[]
end

Base.@kwdef mutable struct PartitionSOCWR <: AbstractPartition
    id::Int
    position::Int = 0

    bus::BusData
    branch::BranchData

    zone_branch_idx_in_full::Vector{Int} = Int[]
    zone_shared_bus_ids::Vector{Int} = Int[]
    branch_idx_full_to_zone::Vector{Int} = Int[]

    shared_bus_local_idx::Vector{Int} = Int[]
    shared_bus_consensus_idx::Vector{Int} = Int[]

    v2_buf::Vector{Float64} = Float64[]
    v2_shared_buf::Vector{Float64} = Float64[]

    model::Union{Nothing, JuMP.Model} = nothing
    v2_var::Union{Nothing, JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}} = nothing

    obj_base::Union{Nothing, Float64, JuMP.AffExpr, JuMP.QuadExpr, JuMP.NonlinearExpr} =
        nothing
    obj_static::JuMP.QuadExpr = zero(JuMP.QuadExpr)

    quad_bus_expr::JuMP.QuadExpr = zero(JuMP.QuadExpr)
    quad_rt_expr::JuMP.QuadExpr = zero(JuMP.QuadExpr)
    quad_pq_expr::JuMP.QuadExpr = zero(JuMP.QuadExpr)

    soft_anchor_expr::JuMP.QuadExpr = zero(JuMP.QuadExpr)

    objective_params::ObjectiveData
    has_solution::Bool = false

    generation_cost::Float64 = 0.0

    zone_edges::Vector{Int} = Int[]
end

# =======================================================================
# ADMM Partitions Builder
# -----------------------------------------------------------------------
# Precomputes `shared_bus_local_idx` and `shared_bus_consensus_idx` at
# build time using the dense Vector index map already available in
# `BusData` and the canonical `full_shared_bus_ids` order.
# =======================================================================

function build_partitions(::Type{P}, sys::DistributedSystem) where {P <: AbstractPartition}
    n = length(sys.zones)
    parts = Vector{P}(undef, n)

    # Global lookup: bus_id -> position in full_shared_bus_ids (0 if absent).
    full_shared = sys.full_shared_bus_ids
    n_shared_global = length(full_shared)
    max_shared_id = n_shared_global == 0 ? 0 : maximum(full_shared)
    full_shared_lookup = zeros(Int, max_shared_id)
    @inbounds for (j, b) in enumerate(full_shared)
        full_shared_lookup[b] = j
    end

    @inbounds for i in 1:n
        zone = sys.zones[i]
        s = sys.parts[i]
        zone_shared_ids = sys.zone_shared_bus_ids[i]

        n_sh = length(zone_shared_ids)
        shared_local_idx = Vector{Int}(undef, n_sh)
        shared_cons_idx = Vector{Int}(undef, n_sh)
        bus_id_to_idx = s.bus.id_to_idx
        for k in 1:n_sh
            b = zone_shared_ids[k]
            shared_local_idx[k] = bus_id_to_idx[b]
            shared_cons_idx[k] = full_shared_lookup[b]
        end

        parts[i] = P(;
            id = zone,
            position = i,
            bus = s.bus,
            branch = s.branch,
            zone_branch_idx_in_full = sys.zone_branch_idx_in_full[i],
            branch_idx_full_to_zone = sys.branch_idx_full_to_zone[i],
            zone_shared_bus_ids = zone_shared_ids,
            shared_bus_local_idx = shared_local_idx,
            shared_bus_consensus_idx = shared_cons_idx,
            objective_params = s.objective,
        )
    end

    return Partitions{P}(parts, sys.zones)
end
