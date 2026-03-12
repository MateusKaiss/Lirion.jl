# ============================================================
# Generator data
# ============================================================

struct GenData
    n::Int
    bus::Vector{Int}

    pg::Vector{Float64}
    qg::Vector{Float64}
    pmin::Vector{Float64}
    pmax::Vector{Float64}
    qmin::Vector{Float64}
    qmax::Vector{Float64}

    status::Vector{Int}
    existing::Vector{Bool}

    pg_start::Vector{Union{Nothing, Float64}}
    qg_start::Vector{Union{Nothing, Float64}}

    plant_type::Vector{String}
    curtailment_cost::Vector{Float64}
end

function get_gen_parameters(case_raw::CaseRaw)
    gen_raw = case_raw.gen::GenRaw{Float64}

    pg = gen_raw.pg
    qg = gen_raw.qg
    status = gen_raw.gen_status

    n = length(gen_raw.gen_bus)

    pg_start = Vector{Union{Nothing, Float64}}(undef, n)
    qg_start = Vector{Union{Nothing, Float64}}(undef, n)

    @inbounds for i in 1:n
        pg_i = pg[i]
        qg_i = qg[i]
        pg_start[i] = isnan(pg_i) ? nothing : pg_i
        qg_start[i] = isnan(qg_i) ? nothing : qg_i
    end

    return (
        n = n,
        bus = gen_raw.gen_bus,
        pg = pg,
        qg = qg,
        pmin = gen_raw.pmin,
        pmax = gen_raw.pmax,
        qmin = gen_raw.qmin,
        qmax = gen_raw.qmax,
        status = status,
        existing = status .== 1,
        pg_start = pg_start,
        qg_start = qg_start,
        plant_type = gen_raw.plant_type,
        curtailment_cost = gen_raw.curtailment_cost,
    )
end

function build_gen_data(case_raw::CaseRaw)::GenData
    params = get_gen_parameters(case_raw)

    return GenData(
        params.n,
        params.bus,
        params.pg,
        params.qg,
        params.pmin,
        params.pmax,
        params.qmin,
        params.qmax,
        params.status,
        params.existing,
        params.pg_start,
        params.qg_start,
        params.plant_type,
        params.curtailment_cost,
    )
end
