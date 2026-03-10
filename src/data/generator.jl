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

    initial_pg::Vector{Union{Nothing, Float64}}
    initial_qg::Vector{Union{Nothing, Float64}}

    plant_type::Vector{String}
    curtailment_cost::Vector{Float64}
end

# ---------------------------------------------------------------------
# Generator Parameters (CaseRaw)
# ---------------------------------------------------------------------
function get_generator_parameters(case_raw::CaseRaw)
    if case_raw.gen === nothing
        return (
            n = 0,
            bus = Int[],
            pg = Float64[],
            qg = Float64[],
            pmin = Float64[],
            pmax = Float64[],
            qmin = Float64[],
            qmax = Float64[],
            status = Int[],
            existing = Bool[],
            initial_pg = Union{Nothing, Float64}[],
            initial_qg = Union{Nothing, Float64}[],
            plant_type = String[],
            curtailment_cost = Float64[],
        )
    end

    gen = case_raw.gen::GenRaw{Float64}
    n = length(gen.gen_bus)

    initial_pg = Vector{Union{Nothing, Float64}}(undef, n)
    initial_qg = Vector{Union{Nothing, Float64}}(undef, n)

    @inbounds for i in 1:n
        pg = gen.pg[i]
        qg = gen.qg[i]
        initial_pg[i] = isnan(pg) ? nothing : pg
        initial_qg[i] = isnan(qg) ? nothing : qg
    end

    return (
        n = n,
        bus = gen.gen_bus,
        pg = gen.pg,
        qg = gen.qg,
        pmin = gen.pmin,
        pmax = gen.pmax,
        qmin = gen.qmin,
        qmax = gen.qmax,
        status = gen.gen_status,
        existing = gen.gen_status .== 1,
        initial_pg = initial_pg,
        initial_qg = initial_qg,
        plant_type = gen.plant_type,
        curtailment_cost = gen.curtailment_cost,
    )
end

function build_gen_data(case_raw::CaseRaw)::GenData
    params = get_generator_parameters(case_raw)

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
        params.initial_pg,
        params.initial_qg,
        params.plant_type,
        params.curtailment_cost,
    )
end
