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
    pmax_ts = gen.pmax_ts

    # We keep the static nameplate (pmax[g]) as the variable's hard upper
    # bound because JuMP bounds can only be scalars. The per-period
    # availability ceiling pmax_ts[g, t] is enforced by the curtailment
    # constraint (for curtailable gens) or by an extra inequality
    # constraint (for non-curtailable gens, only when it is strictly
    # tighter than pmax).
    @variable(
        model,
        p_generation[g = 1:G, t = 1:T; gen.existing[g]],
        lower_bound = pmin[g],
        upper_bound = pmax[g],
    )

    # Warm start clipped to the per-period availability — otherwise a pg
    # from the matpower `pg` field (> pmax_ts[g, t]) puts Ipopt far outside
    # the feasible region and can trigger LOCALLY_INFEASIBLE.
    @inbounds for g in 1:G
        gen.existing[g] || continue
        base = isnan(pg_start[g]) ? pg_fallback[g] : pg_start[g]
        for t in 1:T
            s = min(base, pmax_ts[g, t])
            s = max(s, pmin[g])
            JuMP.set_start_value(p_generation[g, t], s)
        end
    end

    model[:p_generation] = p_generation

    # Time-varying bounds on active power. pmax is the upper cap; pmin is
    # the lower cap (useful for must-run windows). We only add the
    # inequality when the time series value is STRICTLY tighter than the
    # scalar variable bound, to avoid inflating the constraint count.
    pmin_ts = gen.pmin_ts
    for g in 1:G
        gen.existing[g] || continue
        is_renewable = gen.plant_type[g] == "R"
        for t in 1:T
            # Skip renewables' pmax — they get a tight equality via
            # add_curtailment_constraints!.
            if !is_renewable && pmax_ts[g, t] < pmax[g] - 1e-12
                @constraint(model, p_generation[g, t] <= pmax_ts[g, t])
            end
            if pmin_ts[g, t] > pmin[g] + 1e-12
                @constraint(model, p_generation[g, t] >= pmin_ts[g, t])
            end
        end
    end

    return nothing
end

function add_reactive_generation_variables!(model::JuMP.Model, gen::GenData)
    T = get(model.ext, :time_periods, 1)
    G = gen.n
    qg_start = gen.qg_start
    qg_fallback = gen.qg
    qmin = gen.qmin
    qmax = gen.qmax
    qmax_ts = gen.qmax_ts
    qmin_ts = gen.qmin_ts

    @variable(
        model,
        q_generation[g = 1:G, t = 1:T; gen.existing[g]],
        lower_bound = qmin[g],
        upper_bound = qmax[g],
    )

    @inbounds for g in 1:G
        gen.existing[g] || continue
        base = isnan(qg_start[g]) ? qg_fallback[g] : qg_start[g]
        for t in 1:T
            s = min(max(base, qmin[g]), qmax[g])
            JuMP.set_start_value(q_generation[g, t], s)
        end
    end

    model[:q_generation] = q_generation

    # Time-varying reactive envelope. Only tighten when strictly narrower
    # than the variable's scalar bound.
    for g in 1:G
        gen.existing[g] || continue
        for t in 1:T
            if qmax_ts[g, t] < qmax[g] - 1e-12
                @constraint(model, q_generation[g, t] <= qmax_ts[g, t])
            end
            if qmin_ts[g, t] > qmin[g] + 1e-12
                @constraint(model, q_generation[g, t] >= qmin_ts[g, t])
            end
        end
    end

    return nothing
end

function add_curtailment_variables!(model::JuMP.Model, gen::GenData)
    T = get(model.ext, :time_periods, 1)
    pmax_ts = gen.pmax_ts

    @variable(model, curtailment[g in gen.curtailable, t in 1:T] >= 0.0)

    # Curtailment is capped at the per-period availability to keep the
    # feasible region aligned with `p_gen + curt == pmax_ts`.
    @inbounds for g in gen.curtailable, t in 1:T
        JuMP.set_upper_bound(curtailment[g, t], pmax_ts[g, t])
    end

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
