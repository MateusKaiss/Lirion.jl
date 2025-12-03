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
    initial_pg::Dict{Int,Union{Nothing,Float64}}
    initial_qg::Dict{Int,Union{Nothing,Float64}}
    plant_type::Vector{String}
    curtailment_cost::Vector{Float64}
end

function get_generator_parameters(gen_data)
    n = length(gen_data["gen_bus"])
    start_pg = Dict{Int,Union{Nothing,Float64}}()
    start_qg = Dict{Int,Union{Nothing,Float64}}()

    @inbounds for g in 1:n
        pg_val = gen_data["pg"][g]
        qg_val = gen_data["qg"][g]

        start_pg[g] = isnan(pg_val) ? nothing : pg_val
        start_qg[g] = isnan(qg_val) ? nothing : qg_val
    end

    gen_status_vec = Vector{Int}(gen_data["gen_status"])
    existing       = gen_status_vec .== 1

    raw_pt = gen_data["plant_type"]
    plant_type = String.(raw_pt)
    curtailment_cost = if haskey(gen_data, "curtailment_cost")
        Vector{Float64}(gen_data["curtailment_cost"])
    else
        zeros(Float64, n)
    end

    return (
        existing         = existing,
        initial_pg       = start_pg,
        initial_qg       = start_qg,
        plant_type       = plant_type,
        curtailment_cost = curtailment_cost,
        gen_status       = gen_status_vec,
    )
end

# =====================================================================
# Builder
# =====================================================================
function build_gen_data(data) :: GenData
    gen = data["gen"]
    n   = length(gen["gen_bus"])

    params = get_generator_parameters(gen)

    return GenData(
        n,
        Int.(gen["gen_bus"]),
        Float64.(gen["pg"]),
        Float64.(gen["qg"]),
        Float64.(gen["pmin"]),
        Float64.(gen["pmax"]),
        Float64.(gen["qmin"]),
        Float64.(gen["qmax"]),
        params.gen_status,         

        params.existing,           
        params.initial_pg,          
        params.initial_qg,          
        params.plant_type,
        params.curtailment_cost,     
    )
end
