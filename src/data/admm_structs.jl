
struct Partitions{P}
    parts::Vector{P}
    zones::Vector{Int}
    zone_to_idx::Dict{Int, Int}
end

@inline Base.length(ps::Partitions) = length(ps.parts)

@inline function Base.iterate(ps::Partitions, st::Int = 1)
    st > length(ps.parts) && return nothing
    @inbounds return (ps.parts[st], st + 1)
end

@inline get_partition(ps::Partitions, zone::Int) =
    (@inbounds ps.parts[ps.zone_to_idx[zone]])

# =======================================================================
# ADMM Partitions structs
# =======================================================================

Base.@kwdef mutable struct PartitionDC
    id::Int

    bus::BusData
    branch::BranchData

    zone_branch_idx_in_full::Vector{Int} = Int[]
    shared_buses_zone::Vector{Int} = Int[]
    branch_idx_full_to_zone::Dict{Int, Int} = Dict{Int, Int}()

    bus_local_index::Dict{Int, Int} = Dict{Int, Int}()
    consensus_bus_local_idxs::Vector{Int} = Int[]

    theta_buf::Vector{Float64} = Float64[]
    theta_shared_buf::Vector{Float64} = Float64[]

    model::Union{Nothing, JuMP.Model} = nothing
    theta_var::Union{Nothing, JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}} = nothing

    obj_base::Union{Nothing, Float64, JuMP.AffExpr, JuMP.QuadExpr} = nothing
    quad_bus_expr::JuMP.QuadExpr = zero(JuMP.QuadExpr)
    rho_bus_cached::Float64 = 0.0
    quad_branch_expr::JuMP.QuadExpr = zero(JuMP.QuadExpr)
    rho_branch_cached::Float64 = 0.0

    objective_params::ObjectiveData
    has_solution::Bool = false

    generation_cost_model_expr::Float64 = 0.0
    real_generation_cost::Float64 = 0.0

    zone_edges::Vector{Int} = Int[]
end

Base.@kwdef mutable struct PartitionAC
    id::Int

    bus::BusData
    branch::BranchData

    zone_branch_idx_in_full::Vector{Int} = Int[]
    shared_buses_zone::Vector{Int} = Int[]
    branch_idx_full_to_zone::Dict{Int, Int} = Dict{Int, Int}()

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

    rho_bus_theta_cached::Float64 = 0.0
    rho_bus_v_cached::Float64 = 0.0

    quad_branch_expr_ac::Union{Nothing, JuMP.QuadExpr, JuMP.NonlinearExpr} = nothing
    rho_branch_p_cached::Float64 = 0.0
    rho_branch_q_cached::Float64 = 0.0

    objective_params::ObjectiveData
    has_solution::Bool = false

    generation_cost_model_expr::Float64 = 0.0
    real_generation_cost::Float64 = 0.0

    zone_edges::Vector{Int} = Int[]
end

Base.@kwdef mutable struct PartitionSOCP
    id::Int

    bus::BusData
    branch::BranchData

    zone_branch_idx_in_full::Vector{Int} = Int[]
    shared_buses_zone::Vector{Int} = Int[]
    branch_idx_full_to_zone::Dict{Int, Int} = Dict{Int, Int}()

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
    rho_bus_v2_cached::Float64 = 0.0

    quad_rt_expr_socp::JuMP.QuadExpr = zero(JuMP.QuadExpr)
    rho_branch_r_cached::Float64 = 0.0
    rho_branch_t_cached::Float64 = 0.0

    quad_pq_expr_socp::JuMP.QuadExpr = zero(JuMP.QuadExpr)
    rho_branch_p_cached::Float64 = 0.0
    rho_branch_q_cached::Float64 = 0.0

    soft_anchor_expr::JuMP.QuadExpr = zero(JuMP.QuadExpr)

    objective_params::ObjectiveData
    has_solution::Bool = false

    generation_cost_model_expr::Float64 = 0.0
    real_generation_cost::Float64 = 0.0

    zone_edges::Vector{Int} = Int[]
end

# =======================================================================
# ADMM Partitions Builders (meta vem do sys direto)
# =======================================================================

function build_dc_partitions(sys::DistributedSystem)::Partitions{PartitionDC}
    zones = sys.zones
    zone_to_idx = sys.zone_to_idx
    n_zones = length(zones)

    parts = Vector{PartitionDC}(undef, n_zones)
    @inbounds for i in 1:n_zones
        zone = zones[i]
        s = sys.parts[i]

        parts[i] = PartitionDC(;
            id = zone,
            bus = s.bus,
            branch = s.branch,
            zone_branch_idx_in_full = sys.zone_branch_idx_in_full[i],
            branch_idx_full_to_zone = sys.branch_idx_full_to_zone[i],
            shared_buses_zone = sys.shared_buses_zone[i],
            objective_params = s.objective,
        )
    end

    return Partitions{PartitionDC}(parts, zones, zone_to_idx)
end

function build_ac_partitions(sys::DistributedSystem)::Partitions{PartitionAC}
    zones = sys.zones
    zone_to_idx = sys.zone_to_idx
    n_zones = length(zones)

    parts = Vector{PartitionAC}(undef, n_zones)
    @inbounds for i in 1:n_zones
        zone = zones[i]
        s = sys.parts[i]

        parts[i] = PartitionAC(;
            id = zone,
            bus = s.bus,
            branch = s.branch,
            zone_branch_idx_in_full = sys.zone_branch_idx_in_full[i],
            branch_idx_full_to_zone = sys.branch_idx_full_to_zone[i],
            shared_buses_zone = sys.shared_buses_zone[i],
            objective_params = s.objective,
        )
    end

    return Partitions{PartitionAC}(parts, zones, zone_to_idx)
end

function build_socp_partitions(sys::DistributedSystem)::Partitions{PartitionSOCP}
    zones = sys.zones
    zone_to_idx = sys.zone_to_idx
    n_zones = length(zones)

    parts = Vector{PartitionSOCP}(undef, n_zones)
    @inbounds for i in 1:n_zones
        zone = zones[i]
        s = sys.parts[i]

        parts[i] = PartitionSOCP(;
            id = zone,
            bus = s.bus,
            branch = s.branch,
            zone_branch_idx_in_full = sys.zone_branch_idx_in_full[i],
            branch_idx_full_to_zone = sys.branch_idx_full_to_zone[i],
            shared_buses_zone = sys.shared_buses_zone[i],
            objective_params = s.objective,
        )
    end

    return Partitions{PartitionSOCP}(parts, zones, zone_to_idx)
end
