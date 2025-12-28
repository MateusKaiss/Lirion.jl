struct Partitions{P}
    parts::Vector{P}
    zones::Vector{Int}
    zone_to_idx::Dict{Int, Int}
end

Base.length(ps::Partitions) = length(ps.parts)
Base.iterate(ps::Partitions, st::Int = 1) =
    st > length(ps.parts) ? nothing : (ps.parts[st], st + 1)
Base.getindex(ps::Partitions, i::Int) = ps.parts[i]

@inline getzone(ps::Partitions, zone::Int) = ps.parts[ps.zone_to_idx[zone]]
@inline zoneindex(ps::Partitions, zone::Int) = ps.zone_to_idx[zone]

# =======================================================================
# ADMM Partitions structs
# =======================================================================

Base.@kwdef mutable struct PartitionDC
    id::Int

    bus::BusData
    branch::BranchData
    gen::GenData

    branch_local_indices::Vector{Int} = Int[]
    shared_buses_local::Vector{Int} = Int[]
    edge_to_local::Dict{Int, Int} = Dict{Int, Int}()

    bus_local_index::Dict{Int, Int} = Dict{Int, Int}()
    shared_pos_local::Vector{Int} = Int[]

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
    gen::GenData

    branch_local_indices::Vector{Int} = Int[]
    shared_buses_local::Vector{Int} = Int[]
    edge_to_local::Dict{Int, Int} = Dict{Int, Int}()

    bus_local_index::Dict{Int, Int} = Dict{Int, Int}()
    shared_pos_local::Vector{Int} = Int[]

    theta_buf::Vector{Float64} = Float64[]
    vmag_buf::Vector{Float64} = Float64[]
    theta_shared_buf::Vector{Float64} = Float64[]
    vmag_shared_buf::Vector{Float64} = Float64[]

    model::Union{Nothing, JuMP.Model} = nothing
    theta_var::Union{Nothing, JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}} = nothing
    vmag_var::Union{Nothing, JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}} = nothing

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
    gen::GenData

    branch_local_indices::Vector{Int} = Int[]
    shared_buses_local::Vector{Int} = Int[]
    edge_to_local::Dict{Int, Int} = Dict{Int, Int}()

    bus_local_index::Dict{Int, Int} = Dict{Int, Int}()
    shared_pos_local::Vector{Int} = Int[]

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
# ADMM Partitions Builders
# =======================================================================

function build_dc_partitions(sys::ADMMSystem)::Partitions{PartitionDC}
    zones = sort!(collect(keys(sys.parts)))
    parts = Vector{PartitionDC}(undef, length(zones))
    zone_to_idx = Dict{Int, Int}(zone => i for (i, zone) in pairs(zones))

    for (i, zone) in pairs(zones)
        s = sys.parts[zone]
        raw = s.raw

        parts[i] = PartitionDC(;
            id = zone,
            bus = s.bus,
            branch = s.branch,
            gen = s.gen,
            branch_local_indices = raw["branch_local_indices"]::Vector{Int},
            edge_to_local = raw["edge_to_local"]::Dict{Int, Int},
            shared_buses_local = raw["shared_buses_local"]::Vector{Int},
            objective_params = s.objective,
        )
    end

    return Partitions{PartitionDC}(parts, zones, zone_to_idx)
end

function build_ac_partitions(sys::ADMMSystem)::Partitions{PartitionAC}
    zones = sort!(collect(keys(sys.parts)))
    parts = Vector{PartitionAC}(undef, length(zones))
    zone_to_idx = Dict{Int, Int}(zone => i for (i, zone) in pairs(zones))

    for (i, zone) in pairs(zones)
        s = sys.parts[zone]
        raw = s.raw

        parts[i] = PartitionAC(;
            id = zone,
            bus = s.bus,
            branch = s.branch,
            gen = s.gen,
            branch_local_indices = raw["branch_local_indices"]::Vector{Int},
            edge_to_local = raw["edge_to_local"]::Dict{Int, Int},
            shared_buses_local = raw["shared_buses_local"]::Vector{Int},
            objective_params = s.objective,
        )
    end

    return Partitions{PartitionAC}(parts, zones, zone_to_idx)
end

function build_socp_partitions(sys::ADMMSystem)::Partitions{PartitionSOCP}
    zones = sort!(collect(keys(sys.parts)))
    parts = Vector{PartitionSOCP}(undef, length(zones))
    zone_to_idx = Dict{Int, Int}(zone => i for (i, zone) in pairs(zones))

    for (i, zone) in pairs(zones)
        s = sys.parts[zone]
        raw = s.raw

        parts[i] = PartitionSOCP(;
            id = zone,
            bus = s.bus,
            branch = s.branch,
            gen = s.gen,
            branch_local_indices = raw["branch_local_indices"]::Vector{Int},
            edge_to_local = raw["edge_to_local"]::Dict{Int, Int},
            shared_buses_local = raw["shared_buses_local"]::Vector{Int},
            objective_params = s.objective,
        )
    end

    return Partitions{PartitionSOCP}(parts, zones, zone_to_idx)
end