# ------------------------------------------------------------
# Utilities
# ------------------------------------------------------------

function sorted_ids_and_id2idx(ids::Vector{Int})
    sorted = sort(ids)
    id2idx = Dict{Int, Int}(id => i for (i, id) in enumerate(sorted))
    return sorted, id2idx
end

@inline function _zone_cut_edges(part, interzonal_branches::AbstractVector{Int})
    return [e for e in interzonal_branches if haskey(part.branch_global_to_local, e)]
end

# ------------------------------------------------------------
# Hotpaths preparation
# ------------------------------------------------------------

@inline function _prepare_bus_maps!(part)
    ids = part.bus.id
    n_bus = length(ids)

    part.bus_local_index = Dict{Int, Int}(ids[i] => i for i in 1:n_bus)

    shared_bus_ids = part.shared_buses_zone
    n_shared = length(shared_bus_ids)

    part.shared_pos_local = Vector{Int}(undef, n_shared)
    @inbounds for i in 1:n_shared
        part.shared_pos_local[i] = part.bus_local_index[shared_bus_ids[i]]
    end

    return nothing
end

@inline _bus_dims(part) = (length(part.bus.id), length(part.shared_buses_zone))

function prepare_hotpaths!(part::PartitionAC)
    _prepare_bus_maps!(part)
    n_bus, n_shared = _bus_dims(part)
    part.theta_buf = zeros(n_bus)
    part.v_mag_buf = ones(n_bus)
    part.theta_shared_buf = zeros(n_shared)
    part.v_mag_shared_buf = ones(n_shared)
    return nothing
end

function prepare_hotpaths!(part::PartitionDC)
    _prepare_bus_maps!(part)
    n_bus, n_shared = _bus_dims(part)
    part.theta_buf = zeros(n_bus)
    part.theta_shared_buf = zeros(n_shared)
    return nothing
end

function prepare_hotpaths!(part::PartitionSOCP)
    _prepare_bus_maps!(part)
    n_bus, n_shared = _bus_dims(part)
    part.v2_buf = zeros(n_bus)
    part.v2_shared_buf = zeros(n_shared)
    return nothing
end

# ------------------------------------------------------------
# Value helpers (JuMP)
# ------------------------------------------------------------

@inline function _val_or_start(
    m::JuMP.Model,
    x::JuMP.VariableRef,
    default::Float64,
    has_vals::Bool,
)
    if has_vals
        return JuMP.value(x)
    end
    sv = JuMP.start_value(x)
    return sv === nothing ? default : sv
end

@inline _val_at(
    m::JuMP.Model,
    A::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2},
    b::Int,
    t::Int,
    default::Float64,
    has_vals::Bool,
) = _val_or_start(m, A[b, t], default, has_vals)

# ------------------------------------------------------------
# Snapshots
# ------------------------------------------------------------

@inline function _snapshot_bus_var!(
    buf::Vector{Float64},
    m::JuMP.Model,
    A::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2},
    ids::Vector{Int},
    default::Float64,
)
    has_vals = JuMP.has_values(m)
    @inbounds for i in eachindex(ids)
        buf[i] = _val_at(m, A, ids[i], 1, default, has_vals)
    end
    return buf
end

@inline function _snapshot_shared_var!(
    buf::Vector{Float64},
    m::JuMP.Model,
    A::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2},
    ids::Vector{Int},
    shared_pos_local::Vector{Int},
    default::Float64,
)
    has_vals = JuMP.has_values(m)
    @inbounds for i in eachindex(shared_pos_local)
        bidx = shared_pos_local[i]
        buf[i] = _val_at(m, A, ids[bidx], 1, default, has_vals)
    end
    return buf
end

# -------- AC --------

@inline snapshot_theta_ac!(p::PartitionAC) =
    _snapshot_bus_var!(p.theta_buf, p.model, p.theta_var, p.bus.id, 0.0)

@inline snapshot_v_mag!(p::PartitionAC) =
    _snapshot_bus_var!(p.v_mag_buf, p.model, p.v_mag_var, p.bus.id, 1.0)

@inline snapshot_theta_ac_shared!(p::PartitionAC) =
    _snapshot_shared_var!(
        p.theta_shared_buf, p.model, p.theta_var,
        p.bus.id, p.shared_pos_local, 0.0,
    )

@inline snapshot_v_magshared!(p::PartitionAC) =
    _snapshot_shared_var!(
        p.v_mag_shared_buf, p.model, p.v_mag_var,
        p.bus.id, p.shared_pos_local, 1.0,
    )

# -------- DC --------

@inline snapshot_theta!(p::PartitionDC) =
    _snapshot_bus_var!(p.theta_buf, p.model, p.theta_var, p.bus.id, 0.0)

@inline snapshot_theta_shared!(p::PartitionDC) =
    _snapshot_shared_var!(
        p.theta_shared_buf, p.model, p.theta_var,
        p.bus.id, p.shared_pos_local, 0.0,
    )

# -------- SOCP --------

@inline snapshot_v2!(p::PartitionSOCP) =
    _snapshot_bus_var!(p.v2_buf, p.model, p.v2_var, p.bus.id, 1.0)

@inline snapshot_v2_shared!(p::PartitionSOCP) =
    _snapshot_shared_var!(
        p.v2_shared_buf, p.model, p.v2_var,
        p.bus.id, p.shared_pos_local, 1.0,
    )

# ------------------------------------------------------------
# Flow helpers
# ------------------------------------------------------------

@inline function _dc_flow_val_model(
    m::JuMP.Model,
    part::PartitionDC,
    l_local::Int;
    t::Int = 1,
)
    theta = part.theta_var
    f = part.branch.from[l_local]
    tb = part.branch.to[l_local]
    x = part.branch.x[l_local]
    @assert abs(x) > 1e-8 "x≈0 no ramo local $l_local (DC)."

    has_vals = JuMP.has_values(m)
    θf = _val_at(m, theta, f, t, 0.0, has_vals)
    θtb = _val_at(m, theta, tb, t, 0.0, has_vals)

    return (θf - θtb) / x
end

@inline _pq_flow_from_exprs(m::JuMP.Model, l::Int, t::Int = 1) =
    (m[:p_flow_branches][l, t], m[:q_flow_branches][l, t])

@inline _pq_flow_from_exprs(p::PartitionAC, l::Int; t::Int = 1) =
    _pq_flow_from_exprs(p.model, l, t)

@inline _pq_flow_from_exprs(p::PartitionSOCP, l::Int; t::Int = 1) =
    _pq_flow_from_exprs(p.model, l, t)

@inline _rt_from_exprs(m::JuMP.Model, l::Int, t::Int = 1) =
    (m[:R][l, t], m[:T_var][l, t])
