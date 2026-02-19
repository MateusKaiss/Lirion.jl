struct GenData{T <: Real}
    n::Int
    bus::Vector{Int}

    pg::Vector{T}
    qg::Vector{T}
    pmin::Vector{T}
    pmax::Vector{T}
    qmin::Vector{T}
    qmax::Vector{T}

    status::Vector{Int}
    existing::Vector{Bool}

    initial_pg::Vector{Union{Nothing, T}}
    initial_qg::Vector{Union{Nothing, T}}

    plant_type::Vector{String}
    curtailment_cost::Vector{T}
end

function get_generator_parameters(raw::CaseRaw)
    if raw.gen === nothing
        return (
            n = 0,
            bus = Int[],
            pg = Float64[], qg = Float64[],
            pmin = Float64[], pmax = Float64[],
            qmin = Float64[], qmax = Float64[],
            status = Int[],
            existing = Bool[],
            initial_pg = Union{Nothing, Float64}[],
            initial_qg = Union{Nothing, Float64}[],
            plant_type = String[],
            curtailment_cost = Float64[],
        )
    end

    gen = raw.gen::GenRaw{Float64}
    n = length(gen.gen_bus)

    pg = gen.pg
    qg = gen.qg
    pmin = gen.pmin
    pmax = gen.pmax
    qmin = gen.qmin
    qmax = gen.qmax

    status = gen.gen_status
    existing = status .== 1

    initial_pg = Vector{Union{Nothing, Float64}}(undef, n)
    initial_qg = Vector{Union{Nothing, Float64}}(undef, n)
    @inbounds for g in 1:n
        initial_pg[g] = isnan(pg[g]) ? nothing : pg[g]
        initial_qg[g] = isnan(qg[g]) ? nothing : qg[g]
    end

    plant_type = gen.plant_type
    curtailment_cost = gen.curtailment_cost

    return (
        n = n,
        bus = gen.gen_bus,
        pg = pg, qg = qg,
        pmin = pmin, pmax = pmax,
        qmin = qmin, qmax = qmax,
        status = status,
        existing = existing,
        initial_pg = initial_pg,
        initial_qg = initial_qg,
        plant_type = plant_type,
        curtailment_cost = curtailment_cost,
    )
end

function build_gen_data(raw::CaseRaw)::GenData{Float64}
    p = get_generator_parameters(raw)
    return GenData{Float64}(
        p.n,
        p.bus,
        p.pg,
        p.qg,
        p.pmin,
        p.pmax,
        p.qmin,
        p.qmax,
        p.status,
        p.existing,
        p.initial_pg,
        p.initial_qg,
        p.plant_type,
        p.curtailment_cost,
    )
end
