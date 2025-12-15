function prepare_dense_maps(shared_buses::Vector{Int})
    shared_buses_sorted = sort(shared_buses)
    b2pos = Dict{Int, Int}(b => i for (i, b) in enumerate(shared_buses_sorted))
    return shared_buses_sorted, b2pos
end

function prepare_hotpaths!(part::PartitionDC)
    # mapa id de barra -> posição local no vetor (1..n_local)

    # pega o bus_id e salva em bus_local index sua posição
    bus_ids = part.bus.id
    bus_local_index = Dict{Int, Int}()
    @inbounds for (i, b) in enumerate(bus_ids)
        bus_local_index[b] = i
    end
    part.bus_local_index = bus_local_index

    # vetor das posições locais das barras compartilhadas

    # pega o bus_compartilhado e salva em bus_local index sua posição

    local_shared = part.shared_buses_local
    shared_pos_local = Vector{Int}(undef, length(local_shared))
    @inbounds for (k, b) in enumerate(local_shared)
        shared_pos_local[k] = bus_local_index[b]
    end
    part.shared_pos_local = shared_pos_local

    # buffers para snapshots
    #guarda os buffers
    part.theta_buf = zeros(Float64, length(bus_ids))
    part.theta_shared_buf = zeros(Float64, length(shared_pos_local))
    return nothing
end


# ----------------- Hot-path / buffers (bus) -----------------
function prepare_socp_hotpaths!(sub::PartitionSOCP, shared_buses_sorted::Vector{Int})
    # mapas de índice local para acesso rápido
    bus_ids = sub.bus.id
    bus_local_index = Dict{Int, Int}()
    @inbounds for (i, b) in enumerate(bus_ids)
        bus_local_index[b] = i
    end
    sub.bus_local_index = bus_local_index

    # filtra shared da zona
    local_shared = sub.shared_buses_local
    if isempty(local_shared)
        busset = Set(bus_ids)
        sub.shared_buses_local = [b for b in shared_buses_sorted if b in busset]
        local_shared = sub.shared_buses_local
    end

    # posições locais das compartilhadas
    shared_pos_local = Vector{Int}(undef, length(local_shared))
    @inbounds for (k, b) in enumerate(local_shared)
        shared_pos_local[k] = bus_local_index[b]
    end
    sub.shared_pos_local = shared_pos_local

    # buffers
    sub.v2_buf = zeros(Float64, length(bus_ids))
    sub.v2_shared_buf = zeros(Float64, length(shared_pos_local))
    return nothing
end

function prepare_hotpaths_ac!(part::PartitionAC)
    bus_ids = part.bus.id
    bus_local_index = Dict{Int, Int}()
    @inbounds for (i, b) in enumerate(bus_ids)
        bus_local_index[b] = i
    end
    part.bus_local_index = bus_local_index

    local_shared = part.shared_buses_local
    shared_pos_local = Vector{Int}(undef, length(local_shared))
    @inbounds for (k, b) in enumerate(local_shared)
        shared_pos_local[k] = bus_local_index[b]
    end
    part.shared_pos_local = shared_pos_local

    part.theta_buf = zeros(Float64, length(bus_ids))
    part.vmag_buf = ones(Float64, length(bus_ids))
    part.theta_shared_buf = zeros(Float64, length(shared_pos_local))
    part.vmag_shared_buf = ones(Float64, length(shared_pos_local))
    return nothing
end


@inline function _zone_cut_edges(part::PartitionAC, cut_edges::Vector{Int})
    return [e for e in cut_edges if haskey(part.edge_to_local, e)]
end

@inline function _zone_cut_edges(part::PartitionDC, cut_edges::Vector{Int})
    return [e for e in cut_edges if haskey(part.edge_to_local, e)]
end

@inline function _zone_cut_edges(
    sub::Union{PartitionDC, PartitionAC, PartitionSOCP},
    cut_edges::Vector{Int},
)
    return [e for e in cut_edges if haskey(sub.edge_to_local, e)]
end


@inline function _val_or_start_ac(m::JuMP.Model, x::JuMP.VariableRef)
    if JuMP.has_values(m)
        return JuMP.value(x)
    else
        sv = JuMP.start_value(x)
        return sv === nothing ? 0.0 : sv
    end
end

@inline _theta_val_ac(
    m::JuMP.Model,
    θ::JuMP.Containers.DenseAxisArray,
    b::Int,
    t::Int,
) = _val_or_start_ac(m, θ[b, t])

@inline _vmag_val(
    m::JuMP.Model,
    V::JuMP.Containers.DenseAxisArray,
    b::Int,
    t::Int,
) = _val_or_start_ac(m, V[b, t])

@inline function _val_or_start(m::JuMP.Model, x::JuMP.VariableRef)
    if JuMP.has_values(m)
        return JuMP.value(x)
    else
        sv = JuMP.start_value(x)
        return sv === nothing ? 0.0 : sv
    end
end

@inline function _theta_val(
    m::JuMP.Model,
    θ::JuMP.Containers.DenseAxisArray,
    b::Int,
    t::Int,
)
    return _val_or_start(m, θ[b, t])
end

@inline function snapshot_theta!(part::PartitionDC)
    # θ "completo" (todas as barras da zona)
    m = part.model::JuMP.Model
    θ = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    buf = part.theta_buf
    bus_ids = part.bus.id
    @inbounds for i in eachindex(bus_ids)
        buf[i] = _theta_val(m, θ, bus_ids[i], 1)
    end
    return buf
end

@inline function snapshot_theta_shared!(part::PartitionDC)
    # θ só nas barras compartilhadas da zona
    m = part.model::JuMP.Model
    θ = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    buf = part.theta_shared_buf
    pos = part.shared_pos_local
    bus_ids = part.bus.id
    @inbounds for i in eachindex(pos)
        bidx = pos[i]
        buf[i] = _theta_val(m, θ, bus_ids[bidx], 1)
    end
    return buf
end

@inline function snapshot_theta_ac!(part::PartitionAC)
    m = part.model::JuMP.Model
    θ = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    buf = part.theta_buf
    bus_ids = part.bus.id
    @inbounds for i in eachindex(bus_ids)
        buf[i] = _theta_val_ac(m, θ, bus_ids[i], 1)
    end
    return buf
end

@inline function snapshot_vmag!(part::PartitionAC)
    m = part.model::JuMP.Model
    V = part.vmag_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    buf = part.vmag_buf
    bus_ids = part.bus.id
    @inbounds for i in eachindex(bus_ids)
        buf[i] = _vmag_val(m, V, bus_ids[i], 1)
    end
    return buf
end

@inline function snapshot_theta_ac_shared!(part::PartitionAC)
    m = part.model::JuMP.Model
    θ = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    buf = part.theta_shared_buf
    pos = part.shared_pos_local
    bus_ids = part.bus.id
    @inbounds for i in eachindex(pos)
        bidx = pos[i]
        buf[i] = _theta_val_ac(m, θ, bus_ids[bidx], 1)
    end
    return buf
end

@inline function snapshot_vmag_shared!(part::PartitionAC)
    m = part.model::JuMP.Model
    V = part.vmag_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    buf = part.vmag_shared_buf
    pos = part.shared_pos_local
    bus_ids = part.bus.id
    @inbounds for i in eachindex(pos)
        bidx = pos[i]
        buf[i] = _vmag_val(m, V, bus_ids[bidx], 1)
    end
    return buf
end


@inline function snapshot_v2!(sub::PartitionSOCP)
    m = sub.model::JuMP.Model
    V2 = m[:v_magnitude_sq]
    buf = sub.v2_buf
    ids = sub.bus.id
    @inbounds for i in eachindex(ids)
        x = V2[ids[i], 1]
        buf[i] =
            JuMP.has_values(m) ? value(x) :
            (
                JuMP.start_value(x) === nothing ? 1.0 :
                JuMP.start_value(x)
            )
    end
    return buf
end

@inline function snapshot_v2_shared!(sub::PartitionSOCP)
    m = sub.model::JuMP.Model
    V2 = m[:v_magnitude_sq]
    buf = sub.v2_shared_buf
    pos = sub.shared_pos_local
    ids = sub.bus.id
    @inbounds for i in eachindex(pos)
        bidx = pos[i]
        x = V2[ids[bidx], 1]
        buf[i] =
            JuMP.has_values(m) ? value(x) :
            (
                JuMP.start_value(x) === nothing ? 1.0 :
                JuMP.start_value(x)
            )
    end
    return buf
end

@inline function _dc_flow_val_model(
    m::JuMP.Model,
    part::PartitionDC,
    l_local::Int;
    t::Int = 1,
)
    θ = part.theta_var::JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 2}
    f = part.branch.from[l_local]
    tb = part.branch.to[l_local]
    x = part.branch.x[l_local]
    @assert abs(x) > 1e-8 "x≈0 no ramo local $l_local (DC)."
    return (_theta_val(m, θ, f, t) - _theta_val(m, θ, tb, t)) / x
end

@inline function _ac_flow_from_exprs_existing(
    part::PartitionAC,
    l_local::Int;
    t::Int = 1,
)
    m = part.model::JuMP.Model
    return m[:p_flow_branches][l_local, t], m[:q_flow_branches][l_local, t]
end

