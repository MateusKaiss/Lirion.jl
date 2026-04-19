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

    pg_start::Vector{Float64}
    qg_start::Vector{Float64}

    plant_type::Vector{String}
    curtailment_cost::Vector{Float64}

    curtailable::Vector{Int}

    # Time-varying availability ceiling for each generator (size n × T).
    # For a period t, pmax_ts[g, t] is the effective upper bound for
    # p_generation[g, t]. If no gen_timeseries is provided, pmax_ts
    # defaults to pmax[g] replicated across all periods.
    pmax_ts::Matrix{Float64}
end

function get_gen_parameters(case_raw::CaseRaw)
    gen_raw = case_raw.gen::GenRaw{Float64}

    pg = gen_raw.pg
    qg = gen_raw.qg
    status = gen_raw.gen_status

    n = length(gen_raw.gen_bus)

    pg_start = Vector{Float64}(undef, n)
    qg_start = Vector{Float64}(undef, n)

    @inbounds for i in 1:n
        pg_start[i] = pg[i]
        qg_start[i] = qg[i]
    end

    existing = status .== 1
    plant_type = gen_raw.plant_type

    curtailable = Int[]
    @inbounds for g in 1:n
        (existing[g] && plant_type[g] == "R") && push!(curtailable, g)
    end

    T = max(case_raw.T, 1)
    pmax_ts = Matrix{Float64}(undef, n, T)
    @inbounds for g in 1:n, t in 1:T
        pmax_ts[g, t] = gen_raw.pmax[g]
    end

    if case_raw.gen_timeseries !== nothing &&
       !isempty(case_raw.gen_timeseries.time)
        ts = case_raw.gen_timeseries
        @inbounds for i in eachindex(ts.time)
            g = ts.gen_id[i]
            t = ts.time[i]
            (g < 1 || g > n) && continue
            (t < 1 || t > T) && continue
            pmax_ts[g, t] = ts.pmax[i]
        end
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
        existing = existing,
        pg_start = pg_start,
        qg_start = qg_start,
        plant_type = plant_type,
        curtailment_cost = gen_raw.curtailment_cost,
        curtailable = curtailable,
        pmax_ts = pmax_ts,
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
        params.curtailable,
        params.pmax_ts,
    )
end
