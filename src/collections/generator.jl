struct GenData{T<:Real}
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

    initial_pg::Vector{Union{Nothing,T}}
    initial_qg::Vector{Union{Nothing,T}}

    plant_type::Vector{String}
    curtailment_cost::Vector{T}
end

# ---------------------------------------------------------------------
# Generator Parameters
# ---------------------------------------------------------------------
function get_generator_parameters(matpower_data::Dict{String,Any})
    gen_data    = matpower_data["gen"]

    n = length(gen_data["gen_bus"])
    # Concretiza tipos numéricos em Float64 (sem Vector{Any})
    pg   = Float64.(gen_data["pg"])
    qg   = Float64.(gen_data["qg"])
    pmin = Float64.(gen_data["pmin"])
    pmax = Float64.(gen_data["pmax"])
    qmin = Float64.(gen_data["qmin"])
    qmax = Float64.(gen_data["qmax"])

    gen_status_vec = Int.(gen_data["gen_status"])
    existing       = gen_status_vec .== 1

    # Valores iniciais por índice de gerador (1:n), não precisa de Dict
    initial_pg = Vector{Union{Nothing,Float64}}(undef, n)
    initial_qg = Vector{Union{Nothing,Float64}}(undef, n)

    @inbounds for g in 1:n
        pg_val = pg[g]
        qg_val = qg[g]
        initial_pg[g] = isnan(pg_val) ? nothing : pg_val
        initial_qg[g] = isnan(qg_val) ? nothing : qg_val
    end

    raw_pt     = gen_data["plant_type"]
    plant_type = String.(raw_pt)

    curtailment_cost = if haskey(gen_data, "curtailment_cost")
        Float64.(gen_data["curtailment_cost"])
    else
        zeros(Float64, n)
    end

    return (
        n               = n,
        bus             = Int.(gen_data["gen_bus"]),
        pg              = pg,
        qg              = qg,
        pmin            = pmin,
        pmax            = pmax,
        qmin            = qmin,
        qmax            = qmax,
        status          = gen_status_vec,
        existing        = existing,
        initial_pg      = initial_pg,
        initial_qg      = initial_qg,
        plant_type      = plant_type,
        curtailment_cost = curtailment_cost,
    )
end

# =====================================================================
# Builder
# =====================================================================
function build_gen_data(data::Dict{String,Any})::GenData{Float64}
    params = get_generator_parameters(data)

    return GenData{Float64}(
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