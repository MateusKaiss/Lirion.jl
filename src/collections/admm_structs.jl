struct Partitions{P}
    parts::Vector{P}
    zones::Vector{Int}
    zone_to_idx::Dict{Int, Int}
end

Base.length(ps::Partitions) = length(ps.parts)
Base.iterate(ps::Partitions, st::Int = 1) =
    st > length(ps.parts) ? nothing : (ps.parts[st], st + 1)
Base.getindex(ps::Partitions, i::Int) = ps.parts[i]

@inline function getzone(ps::Partitions, zone::Int)
    return ps.parts[ps.zone_to_idx[zone]]
end

@inline function zoneindex(ps::Partitions, zone::Int)
    return ps.zone_to_idx[zone]
end

# =======================================================================
# ADMM Partitions structs
# =======================================================================
mutable struct PartitionDC
    id::Int
    bus::BusData
    branch::BranchData
    gen::GenData
    branch_local_indices::Vector{Int}
    shared_buses_local::Vector{Int}
    edge_to_local::Dict{Int, Int}
    bus_local_index::Dict{Int, Int}
    shared_pos_local::Vector{Int}
    theta_buf::Vector{Float64}
    theta_shared_buf::Vector{Float64}
    model::Union{Nothing, JuMP.Model}
    theta_var::Union{Nothing, JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}}
    obj_base::Union{Nothing, Float64, JuMP.AffExpr, JuMP.QuadExpr}
    quad_bus_expr::JuMP.QuadExpr
    rho_bus_cached::Float64
    quad_branch_expr::JuMP.QuadExpr
    rho_branch_cached::Float64
    objective_params::ObjectiveData
    has_solution::Bool
    generation_cost_model_expr::Float64
    real_generation_cost::Float64
    zone_edges::Vector{Int}
end

mutable struct PartitionAC
    id::Int
    bus::BusData
    branch::BranchData
    gen::GenData
    branch_local_indices::Vector{Int}
    shared_buses_local::Vector{Int}
    edge_to_local::Dict{Int, Int}
    bus_local_index::Dict{Int, Int}
    shared_pos_local::Vector{Int}
    theta_buf::Vector{Float64}
    vmag_buf::Vector{Float64}
    theta_shared_buf::Vector{Float64}
    vmag_shared_buf::Vector{Float64}
    model::Union{Nothing, JuMP.Model}
    theta_var::Union{Nothing, JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}}
    vmag_var::Union{Nothing, JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}}
    obj_base::Union{Nothing, Float64, JuMP.AffExpr, JuMP.QuadExpr, JuMP.NonlinearExpr}
    quad_bus_expr_ac::JuMP.NonlinearExpr
    rho_bus_theta_cached::Float64
    rho_bus_v_cached::Float64
    quad_branch_expr_ac::Union{JuMP.QuadExpr, JuMP.NonlinearExpr}
    rho_branch_p_cached::Float64
    rho_branch_q_cached::Float64
    objective_params::ObjectiveData
    has_solution::Bool
    generation_cost_model_expr::Float64
    real_generation_cost::Float64
    zone_edges::Vector{Int}
end

mutable struct PartitionSOCP
    id::Int
    bus::BusData
    branch::BranchData
    gen::GenData
    branch_local_indices::Vector{Int}
    shared_buses_local::Vector{Int}
    edge_to_local::Dict{Int, Int}
    bus_local_index::Dict{Int, Int}
    shared_pos_local::Vector{Int}
    v2_buf::Vector{Float64}
    v2_shared_buf::Vector{Float64}
    model::Union{Nothing, JuMP.Model}
    v2_var::Union{Nothing, JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}}
    obj_base::Union{Nothing, Float64, JuMP.AffExpr, JuMP.QuadExpr, JuMP.NonlinearExpr}
    obj_static::Union{JuMP.AffExpr, Float64, JuMP.QuadExpr}
    quad_bus_expr_socp::JuMP.QuadExpr
    rho_bus_v2_cached::Float64
    quad_rt_expr_socp::JuMP.QuadExpr
    rho_branch_r_cached::Float64
    rho_branch_t_cached::Float64
    quad_pq_expr_socp::JuMP.QuadExpr
    rho_branch_p_cached::Float64
    rho_branch_q_cached::Float64
    soft_anchor_expr::JuMP.QuadExpr
    objective_params::ObjectiveData
    has_solution::Bool
    generation_cost_model_expr::Float64
    real_generation_cost::Float64
    zone_edges::Vector{Int}
end

# =======================================================================
# ADMM Partitions Builders
# =======================================================================

function build_ac_partitions(sys::ADMMSystem)::Partitions{PartitionAC}
    zones = sort!(collect(keys(sys.parts)))
    parts = Vector{PartitionAC}(undef, length(zones))
    zone_to_idx = Dict{Int, Int}(zone => i for (i, zone) in pairs(zones))

    for (i, zone) in pairs(zones)
        s = sys.parts[zone]
        bus, branch, gen = s.bus, s.branch, s.gen

        raw = s.raw
        branch_local_indices = raw["branch_local_indices"]::Vector{Int}
        edge_to_local = raw["edge_to_local"]::Dict{Int, Int}
        shared_buses_local = raw["shared_buses_local"]::Vector{Int}

        zone_edges = Int[]

        parts[i] = PartitionAC(
            zone,
            bus,
            branch,
            gen,
            branch_local_indices,
            shared_buses_local,
            edge_to_local,
            Dict{Int, Int}(),
            Int[],
            Float64[],
            Float64[],
            Float64[],
            Float64[],
            nothing,
            nothing,
            nothing,
            nothing,
            zero(JuMP.QuadExpr),
            0.0,
            0.0,
            zero(JuMP.QuadExpr),
            0.0,
            0.0,
            s.objective,
            false,
            0.0,
            0.0,
            zone_edges,
        )
    end

    return Partitions{PartitionAC}(parts, zones, zone_to_idx)
end

function build_dc_partitions(sys::ADMMSystem)::Partitions{PartitionDC}
    zones = sort!(collect(keys(sys.parts)))
    parts = Vector{PartitionDC}(undef, length(zones))
    zone_to_idx = Dict{Int, Int}(zone => i for (i, zone) in pairs(zones))

    for (i, zone) in pairs(zones)
        s = sys.parts[zone]
        bus, branch, gen = s.bus, s.branch, s.gen

        raw = s.raw
        branch_local_indices = raw["branch_local_indices"]::Vector{Int}
        edge_to_local = raw["edge_to_local"]::Dict{Int, Int}
        shared_buses_local = raw["shared_buses_local"]::Vector{Int}

        zone_edges = Int[]

        parts[i] = PartitionDC(
            zone,
            bus,
            branch,
            gen,
            branch_local_indices,
            shared_buses_local,
            edge_to_local,
            Dict{Int, Int}(),
            Int[],
            Float64[],
            Float64[],
            nothing,
            nothing,
            nothing,
            zero(JuMP.QuadExpr),
            0.0,
            zero(JuMP.QuadExpr),
            0.0,
            s.objective,
            false,
            0.0,
            0.0,
            zone_edges,
        )
    end

    return Partitions{PartitionDC}(parts, zones, zone_to_idx)
end

function build_socp_partitions(sys::ADMMSystem)::Partitions{PartitionSOCP}
    zones = sort!(collect(keys(sys.parts)))
    parts = Vector{PartitionSOCP}(undef, length(zones))
    zone_to_idx = Dict{Int, Int}(zone => i for (i, zone) in pairs(zones))

    for (i, zone) in pairs(zones)
        s = sys.parts[zone]
        bus, branch, gen = s.bus, s.branch, s.gen

        raw = s.raw
        branch_local_indices = raw["branch_local_indices"]::Vector{Int}
        edge_to_local = raw["edge_to_local"]::Dict{Int, Int}
        shared_buses_local = raw["shared_buses_local"]::Vector{Int}

        zone_edges = Int[]

        parts[i] = PartitionSOCP(
            zone,
            bus,
            branch,
            gen,
            branch_local_indices,
            shared_buses_local,
            edge_to_local,
            Dict{Int, Int}(),
            Int[],
            Float64[],
            Float64[],
            nothing,
            nothing,
            nothing,
            zero(JuMP.AffExpr),
            zero(JuMP.QuadExpr),
            0.0,
            zero(JuMP.QuadExpr),
            0.0,
            0.0,
            zero(JuMP.QuadExpr),
            0.0,
            0.0,
            zero(JuMP.QuadExpr),
            s.objective,
            false,
            0.0,
            0.0,
            zone_edges,
        )
    end

    return Partitions{PartitionSOCP}(parts, zones, zone_to_idx)
end
