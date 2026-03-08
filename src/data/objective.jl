struct ObjectiveData
    number_of_generators::Int
    number_of_batteries::Int
    generator_costs::Dict{Int, Dict{Int, Vector{Float64}}}   # t -> gen_id -> coeffs
    num_polynomial_coefs::Vector{Int}
    curtailment_cost::Vector{Float64}
    plant_type::Vector{String}
    generator_status::Vector{Int}
    number_of_prosumers::Int
    prosumer_costs::Vector{Float64}
    has_timeseries::Bool
    sbase::Float64
end

function get_objective_parameters(
    case_raw::CaseRaw;
    gencost_timeseries::Union{Nothing, GencostTimeseriesRaw{Float64}} = case_raw.gencost_timeseries,
)
    T = case_raw.T
    has_timeseries = T > 1

    ng = case_raw.gen === nothing ? 0 : length(case_raw.gen.gen_bus)
    nbatt = case_raw.storage === nothing ? 0 : length(case_raw.storage.storage_bus)
    npros = case_raw.prosumer === nothing ? 0 : length(case_raw.prosumer.bus)

    generator_costs = Dict{Int, Dict{Int, Vector{Float64}}}()

    if has_timeseries && gencost_timeseries !== nothing && !isempty(gencost_timeseries.time)
        ts = gencost_timeseries
        @inbounds for i in eachindex(ts.time)
            t = ts.time[i]
            g = ts.gen_id[i]                 # MUST be global gen index
            c = ts.cost_coefficients[i]
            get!(generator_costs, t, Dict{Int, Vector{Float64}}())[g] = c
        end
    else
        generator_costs[1] = Dict{Int, Vector{Float64}}()
        if case_raw.gencost !== nothing
            @inbounds for g in 1:length(case_raw.gencost.cost_coefficients)
                generator_costs[1][g] = case_raw.gencost.cost_coefficients[g]
            end
        end
    end

    num_poly = case_raw.gencost === nothing ? Int[] : case_raw.gencost.ncost
    curt = case_raw.gen === nothing ? Float64[] : case_raw.gen.curtailment_cost
    ptype = case_raw.gen === nothing ? String[] : case_raw.gen.plant_type
    gstat = case_raw.gen === nothing ? Int[] : case_raw.gen.gen_status
    pcost = case_raw.prosumer === nothing ? Float64[] : case_raw.prosumer.cost_injection

    return (
        number_of_generators = ng,
        number_of_batteries = nbatt,
        generator_costs = generator_costs,
        num_polynomial_coefs = num_poly,
        curtailment_cost = curt,
        plant_type = ptype,
        generator_status = gstat,
        number_of_prosumers = npros,
        prosumer_costs = pcost,
        has_timeseries = has_timeseries,
    )
end

function build_objective_data(case_raw::CaseRaw)::ObjectiveData
    p = get_objective_parameters(case_raw)
    return ObjectiveData(
        p.number_of_generators,
        p.number_of_batteries,
        p.generator_costs,
        p.num_polynomial_coefs,
        p.curtailment_cost,
        p.plant_type,
        p.generator_status,
        p.number_of_prosumers,
        p.prosumer_costs,
        p.has_timeseries,
        case_raw.baseMVA,
    )
end
