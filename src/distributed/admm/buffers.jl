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

    part.consensus_bus_local_idxs = Vector{Int}(undef, n_shared)
    @inbounds for i in 1:n_shared
        part.consensus_bus_local_idxs[i] = part.bus_local_index[shared_bus_ids[i]]
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
    model::JuMP.Model,
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
    model::JuMP.Model,
    bus_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2},
    b::Int,
    t::Int,
    default::Float64,
    has_vals::Bool,
) = _val_or_start(model, bus_var[b, t], default, has_vals)

# ------------------------------------------------------------
# Snapshots
# ------------------------------------------------------------

@inline function _snapshot_bus_var!(
    buf::Vector{Float64},
    model::JuMP.Model,
    bus_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2},
    bus_ids::Vector{Int},
    default::Float64,
)
    has_vals = JuMP.has_values(model)
    @inbounds for i in eachindex(bus_ids)
        buf[i] = _val_at(model, bus_var, bus_ids[i], 1, default, has_vals)
    end
    return buf
end

# -------- AC --------

@inline snapshot_theta_ac!(part::PartitionAC) =
    _snapshot_bus_var!(part.theta_buf, part.model, part.theta_var, part.bus.id, 0.0)

@inline snapshot_v_mag!(part::PartitionAC) =
    _snapshot_bus_var!(part.v_mag_buf, part.model, part.v_mag_var, part.bus.id, 1.0)

# -------- DC --------

@inline snapshot_theta!(part::PartitionDC) =
    _snapshot_bus_var!(part.theta_buf, part.model, part.theta_var, part.bus.id, 0.0)

# -------- SOCP --------

@inline snapshot_v2!(part::PartitionSOCP) =
    _snapshot_bus_var!(part.v2_buf, part.model, part.v2_var, part.bus.id, 1.0)

@inline function _snapshot_shared_var!(
    buf::Vector{Float64},
    model::JuMP.Model,
    bus_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2},
    zone_bus_ids::Vector{Int},
    consensus_bus_local_idxs::Vector{Int},
    default::Float64,
)
    has_vals = JuMP.has_values(model)
    @inbounds for i in eachindex(consensus_bus_local_idxs)
        consensus_bus_local_idx = consensus_bus_local_idxs[i]
        buf[i] =
            _val_at(
                model,
                bus_var,
                zone_bus_ids[consensus_bus_local_idx],
                1,
                default,
                has_vals,
            )
    end
    return buf
end

# -------- AC: shared --------

@inline snapshot_theta_shared!(part::PartitionAC) =
    _snapshot_shared_var!(
        part.theta_shared_buf, part.model, part.theta_var,
        part.bus.id, part.consensus_bus_local_idxs, 0.0,
    )

@inline snapshot_v_mag_shared!(part::PartitionAC) =
    _snapshot_shared_var!(
        part.v_mag_shared_buf, part.model, part.v_mag_var,
        part.bus.id, part.consensus_bus_local_idxs, 1.0,
    )

# -------- DC: shared --------

@inline snapshot_theta_shared!(part::PartitionDC) =
    _snapshot_shared_var!(
        part.theta_shared_buf, part.model, part.theta_var,
        part.bus.id, part.consensus_bus_local_idxs, 0.0,
    )

# -------- SOCP: shared --------
@inline snapshot_v2_shared!(part::PartitionSOCP) =
    _snapshot_shared_var!(
        part.v2_shared_buf, part.model, part.v2_var,
        part.bus.id, part.consensus_bus_local_idxs, 1.0,
    )

# ------------------------------------------------------------
# Flow helpers
# ------------------------------------------------------------


# ---------------------------------------------------------------------
# NOTE / TODO (DC flow helpers / indexing mismatch — investigate)

# Try to make DC follow the same pattern as AC/SOCP (use model[:p_flow_branches] via _p_flow_from_exprs), 
# but keep _dc_flow_val_model for now to avoid wrong indexing/results to the paper;
# ---------------------------------------------------------------------

@inline function _dc_flow_val_model(
    model::JuMP.Model,
    part::PartitionDC,
    l_local::Int;
    t::Int = 1,
)
    theta = part.theta_var
    f = part.branch.from[l_local]
    tb = part.branch.to[l_local]
    x = part.branch.x[l_local]
    @assert abs(x) > 1e-8 "x≈0 no ramo local $l_local (DC)."

    has_vals = JuMP.has_values(model)
    θf = _val_at(model, theta, f, t, 0.0, has_vals)
    θtb = _val_at(model, theta, tb, t, 0.0, has_vals)

    return (θf - θtb) / x
end

@inline _pq_flow_from_exprs(model::JuMP.Model, l::Int, t::Int = 1) =
    (model[:p_flow_branches][l, t], model[:q_flow_branches][l, t])

@inline _pq_flow_from_exprs(part::PartitionAC, l::Int; t::Int = 1) =
    _pq_flow_from_exprs(part.model, l, t)

@inline _pq_flow_from_exprs(part::PartitionSOCP, l::Int; t::Int = 1) =
    _pq_flow_from_exprs(part.model, l, t)

@inline _rt_from_exprs(model::JuMP.Model, l::Int, t::Int = 1) =
    (model[:R][l, t], model[:T_var][l, t])
