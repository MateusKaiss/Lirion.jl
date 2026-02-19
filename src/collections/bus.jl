struct BusData{T <: Real}
    n::Int
    id::Vector{Int}
    bus_type::Vector{Int}
    pd::Vector{T}
    qd::Vector{T}
    gs::Vector{T}
    bs::Vector{T}
    area::Vector{Int}
    vm::Vector{T}
    va::Vector{T}
    base_kv::Vector{T}
    zone::Vector{Int}
    vmax::Vector{T}
    vmin::Vector{T}

    reference_bus::Union{Nothing, Int}
    id_to_idx::Dict{Int, Int}
    vstart::Vector{Union{Nothing, T}}
    vastart::Vector{Union{Nothing, T}}
end

function get_bus_parameters(bus::BusRaw{T}) where {T <: Real}
    bus_ids = bus.bus_i
    btype = bus.bus_type
    vm = bus.vm
    va = bus.va
    vmin = bus.vmin
    vmax = bus.vmax

    reference_bus_idx = findfirst(==(3), btype)
    reference_bus = reference_bus_idx === nothing ? nothing : bus_ids[reference_bus_idx]

    n = length(bus_ids)
    vstart = Vector{Union{Nothing, T}}(undef, n)
    vastart = Vector{Union{Nothing, T}}(undef, n)

    @inbounds for i in 1:n
        v = vm[i]
        a = va[i]
        vstart[i] = isnan(v) ? nothing : v
        vastart[i] = isnan(a) ? nothing : a
    end

    return (
        bus_ids = bus_ids,
        reference_bus = reference_bus,
        vmin = vmin,
        vmax = vmax,
        vstart = vstart,
        initial_voltage_angle = vastart,
    )
end

function build_bus_data(raw::CaseRaw)::BusData{Float64}
    @assert raw.bus !== nothing "CaseRaw.bus is nothing"
    bus = raw.bus::BusRaw{Float64}

    params = get_bus_parameters(bus)
    bus_ids = params.bus_ids
    n = length(bus_ids)

    id_to_idx = Dict{Int, Int}()
    @inbounds for (i, b) in pairs(bus_ids)
        id_to_idx[Int(b)] = i
    end

    return BusData{Float64}(
        n,
        bus_ids,
        bus.bus_type,
        bus.pd,
        bus.qd,
        bus.gs,
        bus.bs,
        bus.area,
        bus.vm,
        bus.va,
        bus.base_kv,
        bus.zone,
        bus.vmax,
        bus.vmin,
        params.reference_bus,
        id_to_idx,
        params.vstart,
        params.initial_voltage_angle,
    )
end
