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

struct ConsumerRaw{T <: Real} <: AbstractRawTable
    bus::Vector{Int}
    demand::Vector{T}
    q_demand::Vector{T}
end

struct ProsumerRaw{T <: Real} <: AbstractRawTable
    bus::Vector{Int}
    generation::Vector{T}
    demand::Vector{T}
    q_generation::Vector{T}
    q_demand::Vector{T}
    cost_injection::Vector{T}
    inverter_capacity::Vector{T}
end

# -------------------------
# Timeseries case_raw tables
# -------------------------

struct BusTimeseriesRaw{T <: Real} <: AbstractRawTable
    time::Vector{Int}
    bus::Vector{Int}
    pd::Vector{T}
    qd::Vector{T}
    gs::Vector{T}
    bs::Vector{T}
end

struct ConsumerTimeseriesRaw{T <: Real} <: AbstractRawTable
    time::Vector{Int}
    bus::Vector{Int}
    demand::Vector{T}
    q_demand::Vector{T}
end

struct ProsumerTimeseriesRaw{T <: Real} <: AbstractRawTable
    time::Vector{Int}
    bus::Vector{Int}
    p_pv::Vector{T}
    p_load::Vector{T}
    p_curtail::Vector{T}
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

# -------------------------
# Index maps container
# -------------------------

mutable struct CaseIndexMaps
    bus_idx_map::Dict{Int, Int}
    gen_idx_map::Vector{Int}

    bus_to_gens::Dict{Int, Vector{Int}}
    bus_to_branches_from::Dict{Int, Vector{Int}}
    bus_to_branches_to::Dict{Int, Vector{Int}}
    bus_to_neighbors::Dict{Int, Vector{Int}}
    bus_to_prosumers::Dict{Int, Vector{Int}}
    bus_to_batteries::Dict{Int, Vector{Int}}

    function CaseIndexMaps()
        return new(
            Dict{Int, Int}(),
            Int[],
            Dict{Int, Vector{Int}}(),
            Dict{Int, Vector{Int}}(),
            Dict{Int, Vector{Int}}(),
            Dict{Int, Vector{Int}}(),
            Dict{Int, Vector{Int}}(),
            Dict{Int, Vector{Int}}(),
        )
    end
end

mutable struct DistributedRawMeta
    zone_branch_idx_in_full::Vector{Int}
    full_to_zone_branch_idx::Vector{Int}

    zone_shared_bus_ids::Vector{Int}
    full_shared_bus_ids::Vector{Int}

    bus_set::Set{Int}

    function DistributedRawMeta()
        return new(
            Int[],
            Int[],
            Int[],
            Int[],
            Set{Int}(),
        )
    end
end

# -------------------------
# CaseRaw
# -------------------------

mutable struct CaseRaw
    # core meta
    name::String
    source_type::String
    source_version::String
    baseMVA::Float64
    T::Int

    bus::BusRaw{Float64}
    gen::GenRaw{Float64}
    branch::BranchRaw{Float64}
    gencost::GencostRaw{Float64}

    # optional tables
    dcline::Union{Nothing, DCLineRaw{Float64}}
    bus_name::Union{Nothing, BusNameRaw}
    storage::Union{Nothing, StorageRaw{Float64}}
    switch::Union{Nothing, SwitchRaw{Float64}}
    prosumer::Union{Nothing, ProsumerRaw{Float64}}
    consumer::Union{Nothing, ConsumerRaw{Float64}}

    # timeseries
    bus_timeseries::Union{Nothing, BusTimeseriesRaw{Float64}}
    prosumer_timeseries::Union{Nothing, ProsumerTimeseriesRaw{Float64}}
    consumer_timeseries::Union{Nothing, ConsumerTimeseriesRaw{Float64}}
    gencost_timeseries::Union{Nothing, GencostTimeseriesRaw{Float64}}

    # maps/meta
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
            nothing, nothing, nothing, nothing, nothing, nothing,
            nothing, nothing, nothing, nothing,
            CaseIndexMaps(),
            DistributedRawMeta(),
        )
    end
end

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

slice_consumer_ts(ts::ConsumerTimeseriesRaw{T}, order::Vector{Int}) where {T} =
    ConsumerTimeseriesRaw{T}(
        _slice_vec(ts.time, order),
        _slice_vec(ts.bus, order),
        _slice_vec(ts.demand, order),
        _slice_vec(ts.q_demand, order),
    )

slice_prosumer_ts(ts::ProsumerTimeseriesRaw{T}, order::Vector{Int}) where {T} =
    ProsumerTimeseriesRaw{T}(
        _slice_vec(ts.time, order),
        _slice_vec(ts.bus, order),
        _slice_vec(ts.p_pv, order),
        _slice_vec(ts.p_load, order),
        _slice_vec(ts.p_curtail, order),
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
