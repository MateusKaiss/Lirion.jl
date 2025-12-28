# ------------------------------------------------------------
# Utilities
# ------------------------------------------------------------
function sorted_ids_and_id2idx(ids::Vector{Int})
    sorted = sort(ids)
    id2idx = Dict{Int, Int}(id => i for (i, id) in pairs(sorted))
    return sorted, id2idx
end

@inline function _zone_cut_edges(
    part::Union{PartitionDC, PartitionAC, PartitionSOCP},
    cut_edges::Vector{Int},
)
    return [e for e in cut_edges if haskey(part.edge_to_local, e)]
end

# ------------------------------------------------------------
# Value helpers (JuMP)
# ------------------------------------------------------------
@inline function _val_or_start(m::JuMP.Model, x::JuMP.VariableRef; default::Float64 = 0.0)
    if JuMP.has_values(m)
        return JuMP.value(x)
    else
        sv = JuMP.start_value(x)
        return sv === nothing ? default : sv
    end
end

@inline _vmag_val(
    m::JuMP.Model,
    V::JuMP.Containers.DenseAxisArray,
    b::Int,
    t::Int,
) = _val_or_start(m, V[b, t]; default = 1.0)

@inline _v2_val(
    m::JuMP.Model,
    V2::JuMP.Containers.DenseAxisArray,
    b::Int,
    t::Int,
) = _val_or_start(m, V2[b, t]; default = 1.0)

@inline function _theta_val(
    m::JuMP.Model,
    theta::JuMP.Containers.DenseAxisArray,
    b::Int,
    t::Int,
)
    return _val_or_start(m, theta[b, t]; default = 0.0)
end

# ------------------------------------------------------------
# Hotpaths preparation (allocate buffers + mapping)
# ------------------------------------------------------------

function prepare_hotpaths!(part::PartitionDC)
    ids = part.bus.id
    nbus = length(ids)

    # sempre (re)gera para ficar consistente com part.bus atual
    bus_local_index = Dict{Int, Int}(ids[i] => i for i in 1:nbus)
    part.bus_local_index = bus_local_index

    lbs = part.shared_buses_local
    nshared = length(lbs)

    shared_pos_local = Vector{Int}(undef, nshared)
    @inbounds for i in 1:nshared
        shared_pos_local[i] = bus_local_index[lbs[i]]
    end
    part.shared_pos_local = shared_pos_local

    part.theta_buf = zeros(Float64, nbus)
    part.theta_shared_buf = zeros(Float64, nshared)

    return nothing
end

function prepare_hotpaths!(part::PartitionSOCP)
    ids = part.bus.id
    nbus = length(ids)

    # sempre (re)gera
    bus_local_index = Dict{Int, Int}(ids[i] => i for i in 1:nbus)
    part.bus_local_index = bus_local_index

    lbs = part.shared_buses_local
    nshared = length(lbs)

    shared_pos_local = Vector{Int}(undef, nshared)
    @inbounds for i in 1:nshared
        shared_pos_local[i] = bus_local_index[lbs[i]]
    end
    part.shared_pos_local = shared_pos_local

    part.v2_buf = zeros(Float64, nbus)
    part.v2_shared_buf = zeros(Float64, nshared)

    return nothing
end

function prepare_hotpaths_ac!(part::PartitionAC)
    ids = part.bus.id
    nbus = length(ids)

    bus_local_index = Dict{Int, Int}(ids[i] => i for i in 1:nbus)
    part.bus_local_index = bus_local_index

    lbs = part.shared_buses_local
    nshared = length(lbs)

    shared_pos_local = Vector{Int}(undef, nshared)
    @inbounds for i in 1:nshared
        shared_pos_local[i] = bus_local_index[lbs[i]]
    end
    part.shared_pos_local = shared_pos_local

    part.theta_buf = zeros(Float64, nbus)
    part.vmag_buf = ones(Float64, nbus)
    part.theta_shared_buf = zeros(Float64, nshared)
    part.vmag_shared_buf = ones(Float64, nshared)

    return nothing
end

# ------------------------------------------------------------
# Snapshots (copy current model solution into buffers)
# ------------------------------------------------------------

# ---------------- DC ----------------

@inline function snapshot_theta!(part::PartitionDC)
    m = part.model::JuMP.Model
    theta = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    buf = part.theta_buf
    ids = part.bus.id
    @inbounds for i in eachindex(ids)
        buf[i] = _theta_val(m, theta, ids[i], 1)
    end
    return buf
end

@inline function snapshot_theta_shared!(part::PartitionDC)
    m = part.model::JuMP.Model
    theta = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    buf = part.theta_shared_buf
    pos = part.shared_pos_local
    ids = part.bus.id
    @inbounds for i in eachindex(pos)
        bidx = pos[i]
        buf[i] = _theta_val(m, theta, ids[bidx], 1)
    end
    return buf
end

# ---------------- AC ----------------

@inline function snapshot_theta_ac!(part::PartitionAC)
    m = part.model::JuMP.Model
    theta = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    buf = part.theta_buf
    ids = part.bus.id
    @inbounds for i in eachindex(ids)
        buf[i] = _theta_val(m, theta, ids[i], 1)
    end
    return buf
end

@inline function snapshot_vmag!(part::PartitionAC)
    m = part.model::JuMP.Model
    V = part.vmag_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    buf = part.vmag_buf
    ids = part.bus.id
    @inbounds for i in eachindex(ids)
        buf[i] = _vmag_val(m, V, ids[i], 1)
    end
    return buf
end

@inline function snapshot_theta_ac_shared!(part::PartitionAC)
    m = part.model::JuMP.Model
    theta = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    buf = part.theta_shared_buf
    pos = part.shared_pos_local
    ids = part.bus.id
    @inbounds for i in eachindex(pos)
        bidx = pos[i]
        buf[i] = _theta_val(m, theta, ids[bidx], 1)
    end
    return buf
end

@inline function snapshot_vmag_shared!(part::PartitionAC)
    m = part.model::JuMP.Model
    V = part.vmag_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    buf = part.vmag_shared_buf
    pos = part.shared_pos_local
    ids = part.bus.id
    @inbounds for i in eachindex(pos)
        bidx = pos[i]
        buf[i] = _vmag_val(m, V, ids[bidx], 1)
    end
    return buf
end

# ---------------- SOCP (V^2) ----------------

@inline function snapshot_v2!(sub::PartitionSOCP)
    m = sub.model::JuMP.Model
    V2 = sub.v2_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    buf = sub.v2_buf
    ids = sub.bus.id
    @inbounds for i in eachindex(ids)
        buf[i] = _v2_val(m, V2, ids[i], 1)
    end
    return buf
end

@inline function snapshot_v2_shared!(sub::PartitionSOCP)
    m = sub.model::JuMP.Model
    V2 = sub.v2_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    buf = sub.v2_shared_buf
    pos = sub.shared_pos_local
    ids = sub.bus.id
    @inbounds for i in eachindex(pos)
        bidx = pos[i]
        buf[i] = _v2_val(m, V2, ids[bidx], 1)
    end
    return buf
end

# ------------------------------------------------------------
# Flow helpers
# ------------------------------------------------------------

@inline function _dc_flow_val_model(
    m::JuMP.Model,
    part::PartitionDC,
    l_local::Int;
    t::Int = 1,
)
    theta = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    f = part.branch.from[l_local]
    tb = part.branch.to[l_local]
    x = part.branch.x[l_local]
    @assert abs(x) > 1e-8 "x≈0 no ramo local $l_local (DC)."
    return (_theta_val(m, theta, f, t) - _theta_val(m, theta, tb, t)) / x
end

@inline function _ac_flow_from_exprs_existing(
    part::PartitionAC,
    l_local::Int;
    t::Int = 1,
)
    m = part.model::JuMP.Model
    return m[:p_flow_branches][l_local, t], m[:q_flow_branches][l_local, t]
end

@inline function _socp_flow_from_exprs_existing(
    model::JuMP.Model,
    l_local::Int;
    t::Int = 1,
)
    return model[:p_flow_branches][l_local, t], model[:q_flow_branches][l_local, t]
end

@inline function _rt_from_exprs_socp(
    model::JuMP.Model,
    l_local::Int;
    t::Int = 1,
)
    return model[:R][l_local, t], model[:T_var][l_local, t]
end
