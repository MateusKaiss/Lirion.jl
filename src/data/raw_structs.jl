# ============================================================
# Static tables
# ============================================================

abstract type AbstractRawTable end

struct BusNameRaw <: AbstractRawTable
    name::Vector{String}
end

struct BusRaw{T <: Real} <: AbstractRawTable
    bus_i::Vector{Int}
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
    lam_p::Vector{T}
    lam_q::Vector{T}
    mu_vmax::Vector{T}
    mu_vmin::Vector{T}
    is_virtual::Vector{Bool}
end

struct GenRaw{T <: Real} <: AbstractRawTable
    gen_bus::Vector{Int}
    pg::Vector{T}
    qg::Vector{T}
    qmax::Vector{T}
    qmin::Vector{T}
    vg::Vector{T}
    mbase::Vector{T}
    gen_status::Vector{Int}
    pmax::Vector{T}
    pmin::Vector{T}
    pc1::Vector{T}
    pc2::Vector{T}
    qc1min::Vector{T}
    qc1max::Vector{T}
    qc2min::Vector{T}
    qc2max::Vector{T}
    ramp_agc::Vector{T}
    ramp_10::Vector{T}
    ramp_30::Vector{T}
    ramp_q::Vector{T}
    apf::Vector{T}
    mu_pmax::Vector{T}
    mu_pmin::Vector{T}
    mu_qmax::Vector{T}
    mu_qmin::Vector{T}
    plant_type::Vector{String}
    curtailment_cost::Vector{T}
end

struct BranchRaw{T <: Real} <: AbstractRawTable
    f_bus::Vector{Int}
    t_bus::Vector{Int}
    br_r::Vector{T}
    br_x::Vector{T}
    br_b::Vector{T}
    rate_a::Vector{T}
    rate_b::Vector{T}
    rate_c::Vector{T}
    tap::Vector{T}
    shift::Vector{T}
    br_status::Vector{Int}
    angmin::Vector{T}
    angmax::Vector{T}
    pf::Vector{T}
    qf::Vector{T}
    pt::Vector{T}
    qt::Vector{T}
    mu_sf::Vector{T}
    mu_st::Vector{T}
    mu_angmin::Vector{T}
    mu_angmax::Vector{T}
end

struct DCLineRaw{T <: Real} <: AbstractRawTable
    f_bus::Vector{Int}
    t_bus::Vector{Int}
    br_status::Vector{Int}
    pf::Vector{T}
    pt::Vector{T}
    qf::Vector{T}
    qt::Vector{T}
    vf::Vector{T}
    vt::Vector{T}
    pmin::Vector{T}
    pmax::Vector{T}
    qminf::Vector{T}
    qmaxf::Vector{T}
    qmint::Vector{T}
    qmaxt::Vector{T}
    loss0::Vector{T}
    loss1::Vector{T}
    mu_pmin::Vector{T}
    mu_pmax::Vector{T}
    mu_qminf::Vector{T}
    mu_qmaxf::Vector{T}
    mu_qmint::Vector{T}
    mu_qmaxt::Vector{T}
end

struct StorageRaw{T <: Real} <: AbstractRawTable
    storage_bus::Vector{Int}
    ps::Vector{T}
    qs::Vector{T}
    energy::Vector{T}
    energy_rating::Vector{T}
    charge_rating::Vector{T}
    discharge_rating::Vector{T}
    charge_efficiency::Vector{T}
    discharge_efficiency::Vector{T}
    thermal_rating::Vector{T}
    qmin::Vector{T}
    qmax::Vector{T}
    r::Vector{T}
    x::Vector{T}
    p_loss::Vector{T}
    q_loss::Vector{T}
    status::Vector{Int}
end

struct SwitchRaw{T <: Real} <: AbstractRawTable
    f_bus::Vector{Int}
    t_bus::Vector{Int}
    psw::Vector{T}
    qsw::Vector{T}
    state::Vector{Int}
    thermal_rating::Vector{T}
    status::Vector{Int}
end

struct GencostRaw{T <: Real} <: AbstractRawTable
    model::Vector{Int}
    startup::Vector{T}
    shutdown::Vector{T}
    ncost::Vector{Int}
    cost_coefficients::Vector{Vector{T}}
end

# ============================================================
# Timeseries tables
# ============================================================

struct BusTimeseriesRaw{T <: Real} <: AbstractRawTable
    time::Vector{Int}
    bus::Vector{Int}
    pd::Vector{T}
    qd::Vector{T}
    gs::Vector{T}
    bs::Vector{T}
end

struct GencostTimeseriesRaw{T <: Real} <: AbstractRawTable
    time::Vector{Int}
    gen_id::Vector{Int}
    model::Vector{Int}
    startup::Vector{T}
    shutdown::Vector{T}
    ncost::Vector{Int}
    cost_coefficients::Vector{Vector{T}}
    tou::Vector{T}
    feed_in_tariff::Vector{T}
end

struct GenTimeseriesRaw{T <: Real} <: AbstractRawTable
    time::Vector{Int}
    gen_id::Vector{Int}
    pmax::Vector{T}
    pmin::Vector{T}
    qmax::Vector{T}
    qmin::Vector{T}
end

# ============================================================
# Empty constructors
# ============================================================

@inline _empty_vec(::Type{T}) where {T} = T[]

function empty_bus_raw(::Type{T} = Float64) where {T <: Real}
    return BusRaw{T}(
        Int[], Int[],
        _empty_vec(T), _empty_vec(T), _empty_vec(T), _empty_vec(T),
        Int[],
        _empty_vec(T), _empty_vec(T), _empty_vec(T),
        Int[],
        _empty_vec(T), _empty_vec(T),
        _empty_vec(T), _empty_vec(T),
        _empty_vec(T), _empty_vec(T),
        Bool[],
    )
end

function empty_gen_raw(::Type{T} = Float64) where {T <: Real}
    return GenRaw{T}(
        Int[],
        _empty_vec(T), _empty_vec(T), _empty_vec(T), _empty_vec(T),
        _empty_vec(T), _empty_vec(T),
        Int[],
        _empty_vec(T), _empty_vec(T),
        _empty_vec(T), _empty_vec(T),
        _empty_vec(T), _empty_vec(T), _empty_vec(T), _empty_vec(T),
        _empty_vec(T), _empty_vec(T), _empty_vec(T), _empty_vec(T), _empty_vec(T),
        _empty_vec(T), _empty_vec(T), _empty_vec(T), _empty_vec(T),
        String[],
        _empty_vec(T),
    )
end

function empty_branch_raw(::Type{T} = Float64) where {T <: Real}
    return BranchRaw{T}(
        Int[], Int[],
        _empty_vec(T), _empty_vec(T), _empty_vec(T),
        _empty_vec(T), _empty_vec(T), _empty_vec(T),
        _empty_vec(T), _empty_vec(T),
        Int[],
        _empty_vec(T), _empty_vec(T),
        _empty_vec(T), _empty_vec(T), _empty_vec(T), _empty_vec(T),
        _empty_vec(T), _empty_vec(T), _empty_vec(T), _empty_vec(T),
    )
end

function empty_gencost_raw(::Type{T} = Float64) where {T <: Real}
    return GencostRaw{T}(Int[], _empty_vec(T), _empty_vec(T), Int[], Vector{Vector{T}}())
end

# ============================================================
# Index maps and distributed metadata
# ============================================================

mutable struct CaseIndexMaps
    bus_idx_map::Vector{Int}
    gen_idx_map::Vector{Int}

    bus_to_gens::Vector{Vector{Int}}
    bus_to_branches_from::Vector{Vector{Int}}
    bus_to_branches_to::Vector{Vector{Int}}
    bus_to_storages::Vector{Vector{Int}}

    function CaseIndexMaps()
        return new(
            Int[],
            Int[],
            Vector{Vector{Int}}(),
            Vector{Vector{Int}}(),
            Vector{Vector{Int}}(),
            Vector{Vector{Int}}(),
        )
    end
end

mutable struct DistributedRawMeta
    zone_branch_idx_in_full::Vector{Int}
    branch_idx_full_to_zone::Vector{Int}
    zone_shared_bus_ids::Vector{Int}
    full_shared_bus_ids::Vector{Int}
    branch_pairs::Vector{NTuple{4, Int}}
    bus_pairs::Vector{NTuple{3, Int}}
    # SplitBus N-way groups: (bus_id, [zone_pos...]) with the owner zone first.
    # Generalizes `bus_pairs` to a bus shared by more than two zones.
    bus_groups::Vector{Tuple{Int, Vector{Int}}}

    function DistributedRawMeta()
        return new(
            Int[],
            Int[],
            Int[],
            Int[],
            NTuple{4, Int}[],
            NTuple{3, Int}[],
            Tuple{Int, Vector{Int}}[],
        )
    end
end

# ============================================================
# CaseRaw container
# ============================================================

mutable struct CaseRaw
    name::String
    source_type::String
    source_version::String
    baseMVA::Float64
    T::Int

    bus::BusRaw{Float64}
    gen::GenRaw{Float64}
    branch::BranchRaw{Float64}
    gencost::GencostRaw{Float64}

    dcline::Union{Nothing, DCLineRaw{Float64}}
    bus_name::Union{Nothing, BusNameRaw}
    storage::Union{Nothing, StorageRaw{Float64}}
    switch::Union{Nothing, SwitchRaw{Float64}}

    bus_timeseries::Union{Nothing, BusTimeseriesRaw{Float64}}
    gencost_timeseries::Union{Nothing, GencostTimeseriesRaw{Float64}}
    gen_timeseries::Union{Nothing, GenTimeseriesRaw{Float64}}

    maps::CaseIndexMaps
    dist::DistributedRawMeta

    function CaseRaw(; name = "no_name_found", source_type = "matpower",
        source_version = "0.0.0+",
        baseMVA = 1.0, T = 1)
        return new(
            String(name),
            String(source_type),
            String(source_version),
            Float64(baseMVA),
            Int(T),
            empty_bus_raw(Float64),
            empty_gen_raw(Float64),
            empty_branch_raw(Float64),
            empty_gencost_raw(Float64),
            nothing, nothing, nothing, nothing,
            nothing, nothing, nothing,
            CaseIndexMaps(),
            DistributedRawMeta(),
        )
    end
end

# ============================================================
# Timeseries slicing
# ============================================================

@inline function _slice_vec(v::AbstractVector, idx::AbstractVector{<:Integer})
    return v[idx]
end

slice_bus_ts(ts::BusTimeseriesRaw{T}, order::Vector{Int}) where {T} =
    BusTimeseriesRaw{T}(
        _slice_vec(ts.time, order),
        _slice_vec(ts.bus, order),
        _slice_vec(ts.pd, order),
        _slice_vec(ts.qd, order),
        _slice_vec(ts.gs, order),
        _slice_vec(ts.bs, order),
    )

slice_gencost_ts(ts::GencostTimeseriesRaw{T}, order::Vector{Int}) where {T} =
    GencostTimeseriesRaw{T}(
        _slice_vec(ts.time, order),
        _slice_vec(ts.gen_id, order),
        _slice_vec(ts.model, order),
        _slice_vec(ts.startup, order),
        _slice_vec(ts.shutdown, order),
        _slice_vec(ts.ncost, order),
        _slice_vec(ts.cost_coefficients, order),
        _slice_vec(ts.tou, order),
        _slice_vec(ts.feed_in_tariff, order),
    )

slice_gen_ts(ts::GenTimeseriesRaw{T}, order::Vector{Int}) where {T} =
    GenTimeseriesRaw{T}(
        _slice_vec(ts.time, order),
        _slice_vec(ts.gen_id, order),
        _slice_vec(ts.pmax, order),
        _slice_vec(ts.pmin, order),
        _slice_vec(ts.qmax, order),
        _slice_vec(ts.qmin, order),
    )
