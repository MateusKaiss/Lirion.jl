@inline has_curtailment(gen::GenData, g::Int) =
    gen.existing[g] && gen.plant_type[g] == "R"

# =====================================================================
# VARIABLES
# =====================================================================

function add_active_generation_variables!(model::JuMP.Model, gen::GenData)
    T = get(model.ext, :time_periods, 1)
    G = gen.n
    pg_start = gen.pg_start
    pg_fallback = gen.pg
    pmin = gen.pmin
    pmax = gen.pmax

    @variable(
        model,
        p_generation[g = 1:G, t = 1:T; gen.existing[g]],
        lower_bound = pmin[g],
        upper_bound = pmax[g],
        start = isnan(pg_start[g]) ? pg_fallback[g] : pg_start[g],
    )

    model[:p_generation] = p_generation
    return nothing
end

function add_reactive_generation_variables!(model::JuMP.Model, gen::GenData)
    T = get(model.ext, :time_periods, 1)
    G = gen.n
    qg_start = gen.qg_start
    qg_fallback = gen.qg
    qmin = gen.qmin
    qmax = gen.qmax

    @variable(
        model,
        q_generation[g = 1:G, t = 1:T; gen.existing[g]],
        lower_bound = qmin[g],
        upper_bound = qmax[g],
        start = isnan(qg_start[g]) ? qg_fallback[g] : qg_start[g],
    )

    model[:q_generation] = q_generation
    return nothing
end

function add_curtailment_variables!(model::JuMP.Model, gen::GenData)
    T = get(model.ext, :time_periods, 1)

    @variable(model, curtailment[g in gen.curtailable, t in 1:T] >= 0.0)

    model[:curtailment] = curtailment
    return nothing
end

# =====================================================================
# IVR generator variables: current (crg, cig) + power-current link
# =====================================================================

function add_ivr_generator_variables!(model::JuMP.Model, gen::GenData, bus::BusData)
    T = get(model.ext, :time_periods, 1)
    G = gen.n

    # Generator current variables (rectangular)
    @variable(model, crg[g in 1:G, t in 1:T; gen.existing[g]])
    @variable(model, cig[g in 1:G, t in 1:T; gen.existing[g]])

    model[:crg] = crg
    model[:cig] = cig

    # Power-current link: S = V · I*
    #   pg = vr · crg + vi · cig
    #   qg = vi · crg - vr · cig
    vr = model[:v_real]
    vi = model[:v_imag]
    pg = model[:p_generation]
    qg = model[:q_generation]

    @inbounds for g in 1:G
        gen.existing[g] || continue
        bus_id = gen.bus[g]
        for t in 1:T
            @constraint(
                model,
                pg[g, t] == vr[bus_id, t] * crg[g, t] + vi[bus_id, t] * cig[g, t]
            )
            @constraint(
                model,
                qg[g, t] == vi[bus_id, t] * crg[g, t] - vr[bus_id, t] * cig[g, t]
            )
        end
    end

    return nothing
end
