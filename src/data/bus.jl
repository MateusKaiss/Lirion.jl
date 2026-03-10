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

    bus_raw = case_raw.bus::BusRaw{Float64}

    bus_ids = bus_raw.bus_i
    bus_type = bus_raw.bus_type

    reference_bus_idx = findfirst(==(3), bus_type)
    reference_bus = reference_bus_idx === nothing ? nothing : bus_ids[reference_bus_idx]

    n = length(bus_ids)

    vm = bus_raw.vm
    va = bus_raw.va

    vm_start = Vector{Union{Nothing, Float64}}(undef, n)
    va_start = Vector{Union{Nothing, Float64}}(undef, n)

    @inbounds for i in 1:n
        v = vm[i]
        a = va[i]
        vm_start[i] = isnan(v) ? nothing : v
        va_start[i] = isnan(a) ? nothing : a
    end

    id_to_idx = Dict{Int, Int}()
    sizehint!(id_to_idx, n)

    @inbounds for i in 1:n
        id_to_idx[bus_ids[i]] = i
    end

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
        id_to_idx = id_to_idx,
        vm_start = vm_start,
        va_start = va_start,
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
        params.id_to_idx,
        params.vm_start,
        params.va_start,
    )
end
