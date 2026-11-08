# ============================================================
# Abstract type
# ============================================================

abstract type AbstractPartition end

# ============================================================
# Container
# ============================================================

struct Partitions{P <: AbstractPartition}
    parts::Vector{P}
    zones::Vector{Int}
end

@inline Base.length(ps::Partitions) = length(ps.parts)

@inline function Base.iterate(ps::Partitions, st::Int = 1)
    st > length(ps.parts) && return nothing
    @inbounds return (ps.parts[st], st + 1)
end

# ============================================================
# Partition structs
# ============================================================

# ------------------------------------------------------------
# DC
# ------------------------------------------------------------

Base.@kwdef mutable struct PartitionDC <: AbstractPartition
    # identity
    id::Int
    position::Int = 0

    # source data
    bus::BusData
    branch::BranchData
    objective_params::ObjectiveData

    # index maps
    zone_branch_idx_in_full::Vector{Int} = Int[]
    branch_idx_full_to_zone::Vector{Int} = Int[]
    zone_shared_bus_ids::Vector{Int} = Int[]
    shared_bus_local_idx::Vector{Int} = Int[]
    shared_bus_consensus_idx::Vector{Int} = Int[]
    zone_edges::Vector{Int} = Int[]

    # buffers
    theta_buf::Vector{Float64} = Float64[]
    theta_shared_buf::Vector{Float64} = Float64[]

    # local subproblem (LP)
    local_model::Union{Nothing, JuMP.Model} = nothing
    theta_var::Union{Nothing, JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}} = nothing

    # objective expressions
    obj_base::Union{Nothing, Float64, JuMP.AffExpr, JuMP.QuadExpr} = nothing
    obj_static::JuMP.QuadExpr = zero(JuMP.QuadExpr)
    quad_bus_expr::JuMP.QuadExpr = zero(JuMP.QuadExpr)
    quad_branch_expr::JuMP.QuadExpr = zero(JuMP.QuadExpr)
    quad_p_inj_expr::JuMP.QuadExpr = zero(JuMP.QuadExpr)

    # solution status
    generation_cost::Float64 = 0.0
end

# ------------------------------------------------------------
# AC
# ------------------------------------------------------------

Base.@kwdef mutable struct PartitionAC <: AbstractPartition
    # identity
    id::Int
    position::Int = 0

    # source data
    bus::BusData
    branch::BranchData
    objective_params::ObjectiveData

    # index maps
    zone_branch_idx_in_full::Vector{Int} = Int[]
    branch_idx_full_to_zone::Vector{Int} = Int[]
    zone_shared_bus_ids::Vector{Int} = Int[]
    shared_bus_local_idx::Vector{Int} = Int[]
    shared_bus_consensus_idx::Vector{Int} = Int[]
    zone_edges::Vector{Int} = Int[]

    # buffers
    theta_buf::Vector{Float64} = Float64[]
    theta_shared_buf::Vector{Float64} = Float64[]
    v_mag_buf::Vector{Float64} = Float64[]
    v_mag_shared_buf::Vector{Float64} = Float64[]

    # local subproblem (NLP)
    local_model::Union{Nothing, JuMP.Model} = nothing
    theta_var::Union{Nothing, JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}} = nothing
    v_mag_var::Union{Nothing, JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}} = nothing

    # objective expressions
    obj_base::Union{Nothing, Float64, JuMP.AffExpr, JuMP.QuadExpr, JuMP.NonlinearExpr} =
        nothing
    quad_bus_expr::Union{Nothing, JuMP.QuadExpr, JuMP.NonlinearExpr} = nothing
    quad_branch_expr::Union{Nothing, JuMP.QuadExpr, JuMP.NonlinearExpr} = nothing
    quad_p_inj_expr::JuMP.QuadExpr = zero(JuMP.QuadExpr)
    quad_q_inj_expr::JuMP.QuadExpr = zero(JuMP.QuadExpr)

    # solution status
    generation_cost::Float64 = 0.0
end

# ------------------------------------------------------------
# SOCWR
# ------------------------------------------------------------

Base.@kwdef mutable struct PartitionSOCWR <: AbstractPartition
    # identity
    id::Int
    position::Int = 0

    # W-space model for the zone subproblem (SOCWR or SOCWRNL)
    zone_model::Type{<:OPFModel} = SOCWR

    # source data
    bus::BusData
    branch::BranchData
    objective_params::ObjectiveData

    # index maps
    zone_branch_idx_in_full::Vector{Int} = Int[]
    branch_idx_full_to_zone::Vector{Int} = Int[]
    zone_shared_bus_ids::Vector{Int} = Int[]
    shared_bus_local_idx::Vector{Int} = Int[]
    shared_bus_consensus_idx::Vector{Int} = Int[]
    zone_edges::Vector{Int} = Int[]

    # buffers
    v2_buf::Vector{Float64} = Float64[]
    v2_shared_buf::Vector{Float64} = Float64[]

    # local subproblem (SOCP)
    local_model::Union{Nothing, JuMP.Model} = nothing
    v2_var::Union{Nothing, JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}} = nothing

    # objective expressions
    obj_base::Union{Nothing, Float64, JuMP.AffExpr, JuMP.QuadExpr, JuMP.NonlinearExpr} =
        nothing
    obj_static::JuMP.QuadExpr = zero(JuMP.QuadExpr)
    quad_bus_expr::JuMP.QuadExpr = zero(JuMP.QuadExpr)
    quad_rt_expr::JuMP.QuadExpr = zero(JuMP.QuadExpr)
    quad_pq_expr::JuMP.QuadExpr = zero(JuMP.QuadExpr)
    quad_p_inj_expr::JuMP.QuadExpr = zero(JuMP.QuadExpr)
    quad_q_inj_expr::JuMP.QuadExpr = zero(JuMP.QuadExpr)
    tikhonov_reg_expr::JuMP.QuadExpr = zero(JuMP.QuadExpr)

    # solution status
    generation_cost::Float64 = 0.0
end

# ============================================================
# Builder
# ============================================================

function build_partitions(::Type{P}, sys::DistributedSystem) where {P <: AbstractPartition}
    n_zones = length(sys.zones)
    parts = Vector{P}(undef, n_zones)

    full_shared = sys.full_shared_bus_ids
    n_shared_global = length(full_shared)
    max_shared_id = n_shared_global == 0 ? 0 : maximum(full_shared)
    bus_id_to_consensus_idx = zeros(Int, max_shared_id)
    @inbounds for (j, b) in enumerate(full_shared)
        bus_id_to_consensus_idx[b] = j
    end

    @inbounds for i in 1:n_zones
        zone = sys.zones[i]
        s = sys.parts[i]
        zone_shared_ids = sys.zone_shared_bus_ids[i]

        n_sh = length(zone_shared_ids)
        shared_local_idx = Vector{Int}(undef, n_sh)
        shared_cons_idx = Vector{Int}(undef, n_sh)
        bus_id_to_local_idx = s.bus.id_to_idx
        for k in 1:n_sh
            b = zone_shared_ids[k]
            shared_local_idx[k] = bus_id_to_local_idx[b]
            shared_cons_idx[k] = bus_id_to_consensus_idx[b]
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
