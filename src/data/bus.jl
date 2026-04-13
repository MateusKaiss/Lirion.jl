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
    reference_bus::Int
    vm_start::Vector{Float64}
    va_start::Vector{Float64}

    id_to_idx::Vector{Int}
end

@inline has_reference_bus(bus::BusData)::Bool = bus.reference_bus != 0

function _build_id_to_idx(bus_ids::Vector{Int})::Vector{Int}
    n = length(bus_ids)
    n == 0 && return Int[]
    max_id = 0
    @inbounds for i in 1:n
        id = bus_ids[i]
        id > max_id && (max_id = id)
    end
    v = zeros(Int, max_id)
    @inbounds for i in 1:n
        v[bus_ids[i]] = i
    end
    return v
end

function get_bus_parameters(case_raw::CaseRaw)
    @assert case_raw.bus !== nothing "CaseRaw.bus is nothing"

    bus_raw = case_raw.bus::BusRaw{Float64}

    bus_ids = bus_raw.bus_i
    bus_type = bus_raw.bus_type

    reference_bus_idx = findfirst(==(3), bus_type)
    reference_bus = reference_bus_idx === nothing ? 0 : bus_ids[reference_bus_idx]

    n = length(bus_ids)

    vm = bus_raw.vm
    va = bus_raw.va

    vm_start = Vector{Float64}(undef, n)
    va_start = Vector{Float64}(undef, n)

    @inbounds for i in 1:n
        vm_start[i] = vm[i]
        va_start[i] = va[i]
    end

    id_to_idx = _build_id_to_idx(bus_ids)

    return (
        n = n,
        id = bus_ids,
        bus_type = bus_type,
        pd = bus_raw.pd,
        qd = bus_raw.qd,
        gs = bus_raw.gs,
        bs = bus_raw.bs,
        area = bus_raw.area,
        vm = vm,
        va = va,
        base_kv = bus_raw.base_kv,
        zone = bus_raw.zone,
        v_max = bus_raw.vmax,
        v_min = bus_raw.vmin,
        reference_bus = reference_bus,
        vm_start = vm_start,
        va_start = va_start,
        id_to_idx = id_to_idx,
    )
end

function build_bus_data(case_raw::CaseRaw)::BusData
    params = get_bus_parameters(case_raw)

    return BusData(
        params.n,
        params.id,
        params.bus_type,
        params.pd,
        params.qd,
        params.gs,
        params.bs,
        params.area,
        params.vm,
        params.va,
        params.base_kv,
        params.zone,
        params.v_max,
        params.v_min,
        params.reference_bus,
        params.vm_start,
        params.va_start,
        params.id_to_idx,
    )
end
