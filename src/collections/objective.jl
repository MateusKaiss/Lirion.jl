struct ObjectiveData
    number_of_generators::Int
    number_of_batteries::Int
    generator_costs::Dict{Int, Dict{Int, Vector{Float64}}}
    num_polynomial_coefs::Vector
    curtailment_cost::Vector
    plant_type::Vector
    generator_status::Vector
    number_of_prosumers::Int
    prosumer_costs::Vector
    has_timeseries::Bool
    sbase::Float64
end

function get_objective_parameters(
    data::Dict{String, Any},
    gen_data,
    battery_data,
    gencost_data,
    prosumer_data = Dict(),
    gencost_timeseries = nothing,
)
    T = data["T"]
    has_timeseries = T > 1

    generator_costs = Dict{Int, Dict{Int, Vector{Float64}}}()

    if has_timeseries && gencost_timeseries !== nothing
        gen_ids = gencost_timeseries["gen_id"]
        times = gencost_timeseries["time"]
        cost_coeffs = gencost_timeseries["cost_coefficients"]

        for i in eachindex(gen_ids)
            t = times[i]
            g = gen_ids[i]
            c = cost_coeffs[i]

            if !haskey(generator_costs, t)
                generator_costs[t] = Dict{Int, Vector{Float64}}()
            end
            generator_costs[t][g] = c
        end
    else
        generator_costs[1] = Dict{Int, Vector{Float64}}()
        for (g, coefs) in enumerate(getfieldsafe(gencost_data, "cost_coefficients", []))
            generator_costs[1][g] = coefs
        end
    end

    return (
        number_of_generators = getfieldsafe(gen_data, "gen_bus", []) |> length,
        number_of_batteries = getfieldsafe(battery_data, "storage_bus", []) |> length,
        generator_costs = generator_costs,
        num_polynomial_coefs = getfieldsafe(gencost_data, "ncost", []),
        curtailment_cost = getfieldsafe(gen_data, "curtailment_cost", []),
        plant_type = getfieldsafe(gen_data, "plant_type", []),
        generator_status = getfieldsafe(gen_data, "gen_status", []),
        number_of_prosumers = getfieldsafe(prosumer_data, "bus", []) |> length,
        prosumer_costs = getfieldsafe(prosumer_data, "cost_injection", []),
        has_timeseries = has_timeseries,
    )
end

function build_objective_data(data::Dict{String, Any})::ObjectiveData
    params = get_objective_parameters(
        data,
        data["gen"],
        get(data, "storage", Dict()),
        data["gencost"],
        get(data, "prosumer", Dict()),
        get(data, "gencost_timeseries", nothing),
    )

    return ObjectiveData(
        params.number_of_generators,
        params.number_of_batteries,
        params.generator_costs,
        params.num_polynomial_coefs,
        params.curtailment_cost,
        params.plant_type,
        params.generator_status,
        params.number_of_prosumers,
        params.prosumer_costs,
        params.has_timeseries,
        data["baseMVA"],
    )
end
