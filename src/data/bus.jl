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
    v_max::Vector{Float64}
    v_min::Vector{Float64}
    reference_bus::Union{Nothing, Int}

    id_to_idx::Dict{Int, Int}
    vm_start::Vector{Union{Nothing, Float64}}
    va_start::Vector{Union{Nothing, Float64}}
end

function get_bus_parameters(case_raw::CaseRaw)
    @assert case_raw.bus !== nothing "CaseRaw.bus is nothing"

    bus = case_raw.bus::BusRaw{Float64}

    bus_ids = bus.bus_i
    bus_type = bus.bus_type
    reference_bus_idx = findfirst(==(3), bus_type)
    reference_bus = reference_bus_idx === nothing ? nothing : bus_ids[reference_bus_idx]

    n = length(bus_ids)
    vm = bus.vm
    va = bus.va
    vm_start = Vector{Union{Nothing, Float64}}(undef, n)
    va_start = Vector{Union{Nothing, Float64}}(undef, n)

    @inbounds for i in 1:n
        v = vm[i]
        a = va[i]
        vm_start[i] = isnan(v) ? nothing : v
        va_start[i] = isnan(a) ? nothing : a
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
        v_max = bus.vmax,
        v_min = bus.vmin,
        reference_bus = reference_bus,
        
        id_to_idx = id_to_idx,
        vm_start = vm_start,
        va_start = va_start,
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
        p.v_max,
        p.v_min,
        p.reference_bus,
        p.id_to_idx,
        p.vm_start,
        p.va_start,
    )
end
