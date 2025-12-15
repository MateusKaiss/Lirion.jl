#############################
# MAPAS DENSOS E CONSENSO   #
#############################

function prepare_dense_maps(shared_buses::Vector{Int})
    shared_buses_sorted = sort(shared_buses)
    b2pos = Dict{Int, Int}(b => i for (i, b) in enumerate(shared_buses_sorted))
    return shared_buses_sorted, b2pos
end

#############################
# PARTITIONDC: HOTPATHS     #
# SNAPSHOTS                 #
#############################
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


#############################
# HOTPATHS & SNAPSHOTS AC   #
#############################

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