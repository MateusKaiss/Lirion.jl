# ============================================================
# ObjectiveData
# ============================================================

struct ObjectiveData
    number_of_generators::Int
    number_of_batteries::Int
    generator_costs::Vector{Vector{Vector{Float64}}}
    num_polynomial_coefs::Vector{Int}
    curtailment_cost::Vector{Float64}
    plant_type::Vector{String}
    generator_status::Vector{Int}
    has_timeseries::Bool
    sbase::Float64
end

@inline function _empty_generator_costs(T::Int, ng::Int)::Vector{Vector{Vector{Float64}}}
    return [[Float64[] for _ in 1:ng] for _ in 1:T]
end

function get_objective_parameters(
    case_raw::CaseRaw;
    gencost_timeseries::Union{Nothing, GencostTimeseriesRaw{Float64}} = case_raw.gencost_timeseries,
)
    T = case_raw.T
    has_timeseries = T > 1

    ng = case_raw.gen === nothing ? 0 : length(case_raw.gen.gen_bus)
    nbatt = case_raw.storage === nothing ? 0 : length(case_raw.storage.storage_bus)

    T_eff = max(T, 1)
    generator_costs = _empty_generator_costs(T_eff, ng)

    if has_timeseries && gencost_timeseries !== nothing && !isempty(gencost_timeseries.time)
        ts = gencost_timeseries
        @inbounds for i in eachindex(ts.time)
            t = ts.time[i]
            g = ts.gen_id[i]                 # 1-based global gen index
            (t < 1 || t > T_eff) && continue
            (g < 1 || g > ng) && continue
            generator_costs[t][g] = ts.cost_coefficients[i]
        end
    else
        if case_raw.gencost !== nothing
            coeffs = case_raw.gencost.cost_coefficients
            @inbounds for g in eachindex(coeffs)
                generator_costs[1][g] = coeffs[g]
            end
        end
    end

    num_poly = case_raw.gencost === nothing ? Int[] : case_raw.gencost.ncost
    curt = case_raw.gen === nothing ? Float64[] : case_raw.gen.curtailment_cost
    ptype = case_raw.gen === nothing ? String[] : case_raw.gen.plant_type
    gstat = case_raw.gen === nothing ? Int[] : case_raw.gen.gen_status

    return (
        number_of_generators = ng,
        number_of_batteries = nbatt,
        generator_costs = generator_costs,
        num_polynomial_coefs = num_poly,
        curtailment_cost = curt,
        plant_type = ptype,
        generator_status = gstat,
        has_timeseries = has_timeseries,
    )
end

function build_objective_data(case_raw::CaseRaw)::ObjectiveData
    params = get_objective_parameters(case_raw)
    return ObjectiveData(
        params.number_of_generators,
        params.number_of_batteries,
        params.generator_costs,
        params.num_polynomial_coefs,
        params.curtailment_cost,
        params.plant_type,
        params.generator_status,
        params.has_timeseries,
        case_raw.baseMVA,
    )
end
