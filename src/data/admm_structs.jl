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

# =======================================================================
# ADMM Partitions structs
# =======================================================================

Base.@kwdef mutable struct PartitionDC <: AbstractPartition
    id::Int

    bus::BusData
    branch::BranchData

    zone_branch_idx_in_full::Vector{Int} = Int[]
    zone_shared_bus_ids::Vector{Int} = Int[]
    branch_idx_full_to_zone::Vector{Int} = Int[]

    bus_local_index::Dict{Int, Int} = Dict{Int, Int}()
    consensus_bus_local_idxs::Vector{Int} = Int[]

    theta_buf::Vector{Float64} = Float64[]
    theta_shared_buf::Vector{Float64} = Float64[]

    model::Union{Nothing, JuMP.Model} = nothing
    theta_var::Union{Nothing, JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}} = nothing

    obj_base::Union{Nothing, Float64, JuMP.AffExpr, JuMP.QuadExpr} = nothing
    quad_bus_expr::JuMP.QuadExpr = zero(JuMP.QuadExpr)
    quad_branch_expr::JuMP.QuadExpr = zero(JuMP.QuadExpr)

    objective_params::ObjectiveData
    has_solution::Bool = false

    generation_cost::Float64 = 0.0

    zone_edges::Vector{Int} = Int[]
end

Base.@kwdef mutable struct PartitionAC <: AbstractPartition
    id::Int

    bus::BusData
    branch::BranchData

    zone_branch_idx_in_full::Vector{Int} = Int[]
    zone_shared_bus_ids::Vector{Int} = Int[]
    branch_idx_full_to_zone::Vector{Int} = Int[]

    bus_local_index::Dict{Int, Int} = Dict{Int, Int}()
    consensus_bus_local_idxs::Vector{Int} = Int[]

    theta_buf::Vector{Float64} = Float64[]
    v_mag_buf::Vector{Float64} = Float64[]
    theta_shared_buf::Vector{Float64} = Float64[]
    v_mag_shared_buf::Vector{Float64} = Float64[]

    model::Union{Nothing, JuMP.Model} = nothing
    theta_var::Union{Nothing, JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}} = nothing
    v_mag_var::Union{Nothing, JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}} = nothing

    obj_base::Union{Nothing, Float64, JuMP.AffExpr, JuMP.QuadExpr, JuMP.NonlinearExpr} =
        nothing
    quad_bus_expr_ac::Union{Nothing, JuMP.NonlinearExpr} = nothing
    quad_branch_expr_ac::Union{Nothing, JuMP.QuadExpr, JuMP.NonlinearExpr} = nothing

    objective_params::ObjectiveData
    has_solution::Bool = false

    generation_cost::Float64 = 0.0

    zone_edges::Vector{Int} = Int[]
end

Base.@kwdef mutable struct PartitionSOCP <: AbstractPartition
    id::Int

    bus::BusData
    branch::BranchData

    zone_branch_idx_in_full::Vector{Int} = Int[]
    zone_shared_bus_ids::Vector{Int} = Int[]
    branch_idx_full_to_zone::Vector{Int} = Int[]

    bus_local_index::Dict{Int, Int} = Dict{Int, Int}()
    consensus_bus_local_idxs::Vector{Int} = Int[]

    v2_buf::Vector{Float64} = Float64[]
    v2_shared_buf::Vector{Float64} = Float64[]

    model::Union{Nothing, JuMP.Model} = nothing
    v2_var::Union{Nothing, JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}} = nothing

    obj_base::Union{Nothing, Float64, JuMP.AffExpr, JuMP.QuadExpr, JuMP.NonlinearExpr} =
        nothing
    obj_static::Union{JuMP.AffExpr, Float64, JuMP.QuadExpr} = 0.0

    quad_bus_expr_socp::JuMP.QuadExpr = zero(JuMP.QuadExpr)
    quad_rt_expr_socp::JuMP.QuadExpr = zero(JuMP.QuadExpr)
    quad_pq_expr_socp::JuMP.QuadExpr = zero(JuMP.QuadExpr)

    soft_anchor_expr::JuMP.QuadExpr = zero(JuMP.QuadExpr)

    objective_params::ObjectiveData
    has_solution::Bool = false

    generation_cost::Float64 = 0.0

    zone_edges::Vector{Int} = Int[]
end

# =======================================================================
# ADMM Partitions Builder
# =======================================================================

function build_partitions(::Type{P}, sys::DistributedSystem) where {P <: AbstractPartition}
    n = length(sys.zones)
    parts = Vector{P}(undef, n)

    @inbounds for i in 1:n
        zone = sys.zones[i]
        s = sys.parts[i]

        parts[i] = P(;
            id = zone,
            bus = s.bus,
            branch = s.branch,
            zone_branch_idx_in_full = sys.zone_branch_idx_in_full[i],
            branch_idx_full_to_zone = sys.branch_idx_full_to_zone[i],
            zone_shared_bus_ids = sys.zone_shared_bus_ids[i],
            objective_params = s.objective,
        )
    end

    return Partitions{P}(parts, sys.zones)
end
