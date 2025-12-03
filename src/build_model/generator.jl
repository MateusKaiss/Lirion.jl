# =====================================================================
# Helpers para variáveis de geração / curtailment
# =====================================================================

"Retorna se o gerador g é curtível (renovável) e está existente.
   Mesma lógica antiga: existing[g] == true && plant_type[g] == \"R\"."
function has_curtailment(gen::GenData, g::Int)
    return gen.existing[g] && gen.plant_type[g] == "R"
end

"Constrói vetor de start a partir do Dict de iniciais, caindo para fallback se nothing/NaN."
function _build_start_vector(
    initial_dict::Dict{Int,Union{Nothing,Float64}},
    fallback::AbstractVector{<:Real},
)
    n = length(fallback)
    starts = zeros(Float64, n)
    @inbounds for g in 1:n
        v = get(initial_dict, g, nothing)
        if v === nothing || isnan(v)
            starts[g] = fallback[g]
        else
            starts[g] = v
        end
    end
    return starts
end

# =====================================================================
# Variáveis de geração (p/q) e curtailment
# =====================================================================

function create_active_generation_variables!(model::JuMP.Model, gen::GenData)
    T = get(model.ext, :time_periods, 1)
    G = gen.n

    pg_start = _build_start_vector(gen.initial_pg, gen.pg)

    @variable(
        model,
        p_generation[g = 1:G, t = 1:T; gen.existing[g]],
        lower_bound = gen.pmin[g],
        upper_bound = gen.pmax[g],
        start       = pg_start[g],
    )

    model[:p_generation] = p_generation
    return nothing
end

function create_reactive_generation_variables!(model::JuMP.Model, gen::GenData)
    T = get(model.ext, :time_periods, 1)
    G = gen.n

    qg_start = _build_start_vector(gen.initial_qg, gen.qg)

    @variable(
        model,
        q_generation[g = 1:G, t = 1:T; gen.existing[g]],
        lower_bound = gen.qmin[g],
        upper_bound = gen.qmax[g],
        start       = qg_start[g],
    )

    model[:q_generation] = q_generation
    return nothing
end

function create_curtailment_variables!(model::JuMP.Model, gen::GenData)
    T = get(model.ext, :time_periods, 1)
    G = gen.n

    curtailable = [g for g in 1:G if has_curtailment(gen, g)]

    @variable(model, curtailment[g in curtailable, t in 1:T] >= 0.0)

    model[:curtailment] = curtailment
    return nothing
end

function create_curtailment_constraints!(model::JuMP.Model, gen::GenData)
    T = get(model.ext, :time_periods, 1)
    G = gen.n

    curtailable = [g for g in 1:G if has_curtailment(gen, g)]

    @constraint(model, [g in curtailable, t in 1:T],
        model[:p_generation][g, t] + model[:curtailment][g, t] == gen.pmax[g]
    )

    return nothing
end

# =======================================================================
# Dispatch by model type (DC / AC / SOCP)
# =======================================================================
function generators!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{DC},
    ::Type{SubproblemBuild},
)
    gen = sys.gen
    if gen.n == 0
        @info "No generators in SystemData (DC). Skipping generator modeling."
        return nothing
    end

    create_active_generation_variables!(model, gen)
    create_curtailment_variables!(model, gen)
    create_curtailment_constraints!(model, gen)
    return nothing
end

function generators!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{AC},
    ::Type{SubproblemBuild},
)
    gen = sys.gen
    if gen.n == 0
        @info "No generators in SystemData (AC). Skipping generator modeling."
        return nothing
    end

    create_active_generation_variables!(model, gen)
    create_reactive_generation_variables!(model, gen)
    create_curtailment_variables!(model, gen)
    create_curtailment_constraints!(model, gen)
    return nothing
end

function generators!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{SOCP},
    ::Type{SubproblemBuild},
)
    gen = sys.gen
    if gen.n == 0
        @info "No generators in SystemData (SOCP). Skipping generator modeling."
        return nothing
    end

    create_active_generation_variables!(model, gen)
    create_reactive_generation_variables!(model, gen)
    create_curtailment_variables!(model, gen)
    create_curtailment_constraints!(model, gen)
    return nothing
end