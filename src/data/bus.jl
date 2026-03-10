# ============================================================
# Bus data
# ============================================================

struct BusData
    n::Int
    id::Vector{Int}
    bus_type::Vector{Int}
    pd::Vector{Float64}
    qd::Vector{Float64}
    gs::Vector{Float64}
    bs::Vector{Float64}
    area::Vector{Int}
    vm::Vector{Float64}
    va::Vector{Float64}
    base_kv::Vector{Float64}
    zone::Vector{Int}
    vmax::Vector{Float64}
    vmin::Vector{Float64}

    reference_bus::Union{Nothing, Int}
    id_to_idx::Dict{Int, Int}
    vstart::Vector{Union{Nothing, Float64}}
    vastart::Vector{Union{Nothing, Float64}}
end

function get_bus_parameters(case_raw::CaseRaw)
    @assert case_raw.bus !== nothing "CaseRaw.bus is nothing"

    bus = case_raw.bus::BusRaw{Float64}

    bus_ids = bus.bus_i
    btype = bus.bus_type
    vm = bus.vm
    va = bus.va

    reference_bus_idx = findfirst(==(3), btype)
    reference_bus = reference_bus_idx === nothing ? nothing : bus_ids[reference_bus_idx]

    n = length(bus_ids)

    vstart = Vector{Union{Nothing, Float64}}(undef, n)
    vastart = Vector{Union{Nothing, Float64}}(undef, n)

    @inbounds for i in 1:n
        v = vm[i]
        a = va[i]
        vstart[i] = isnan(v) ? nothing : v
        vastart[i] = isnan(a) ? nothing : a
    end

    id_to_idx = Dict{Int, Int}()
    @inbounds for (i, b) in pairs(bus_ids)
        id_to_idx[b] = i
    end

    return (
        n = n,
        id = bus_ids,
        bus_type = bus.bus_type,
        pd = bus.pd,
        qd = bus.qd,
        gs = bus.gs,
        bs = bus.bs,
        area = bus.area,
        vm = bus.vm,
        va = bus.va,
        base_kv = bus.base_kv,
        zone = bus.zone,
        vmax = bus.vmax,
        vmin = bus.vmin,
        reference_bus = reference_bus,
        id_to_idx = id_to_idx,
        vstart = vstart,
        vastart = vastart,
    )
end

function build_bus_data(case_raw::CaseRaw)::BusData
    p = get_bus_parameters(case_raw)

    return BusData(
        p.n,
        p.id,
        p.bus_type,
        p.pd,
        p.qd,
        p.gs,
        p.bs,
        p.area,
        p.vm,
        p.va,
        p.base_kv,
        p.zone,
        p.vmax,
        p.vmin,
        p.reference_bus,
        p.id_to_idx,
        p.vstart,
        p.vastart,
    )
end
