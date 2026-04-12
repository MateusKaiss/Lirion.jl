function add_ivr_current_balance_constraints!(
    model::JuMP.Model,
    sys::SystemData,
)
    T = get(model.ext, :time_periods, 1)
    pf = sys.powerflow
    bus = sys.bus
    gen = sys.gen

    vr = model[:v_real]
    vi = model[:v_imag]
    crg = model[:crg]
    cig = model[:cig]
    cr_fr = model[:cr_fr]
    ci_fr = model[:ci_fr]
    cr_to = model[:cr_to]
    ci_to = model[:ci_to]

    valid_bus_ids = pf.non_virtual_bus_ids
    bus_idx_map = pf.bus_idx_map
    gs = pf.shunt_conductance
    bs = pf.shunt_susceptance

    @inbounds for b in valid_bus_ids
        bus_idx = bus_idx_map[b]

        # Demand at this bus (static or timeseries snapshot)
        pd_val = pf.has_bus_ts ? pf.pd_ts[bus_idx, 1] : bus.pd[bus_idx]
        qd_val = pf.has_bus_ts ? pf.qd_ts[bus_idx, 1] : bus.qd[bus_idx]

        # Consumer demand
        pd_cons = pf.has_consumer_ts ? pf.consumer_pd_ts[bus_idx, 1] : 0.0
        qd_cons = pf.has_consumer_ts ? pf.consumer_qd_ts[bus_idx, 1] : 0.0

        pd_total = pd_val + pd_cons
        qd_total = qd_val + qd_cons

        gs_val = gs[bus_idx]
        bs_val = bs[bus_idx]

        gen_at_bus = pf.bus_to_generators[bus_idx]
        br_from = pf.bus_to_branches_from[bus_idx]
        br_to = pf.bus_to_branches_to[bus_idx]

        for t in 1:T
            vr_b = vr[b, t]
            vi_b = vi[b, t]

            # Sum generator currents at this bus
            sum_crg = sum(crg[g, t] for g in gen_at_bus if gen.existing[g]; init = 0.0)
            sum_cig = sum(cig[g, t] for g in gen_at_bus if gen.existing[g]; init = 0.0)

            # Sum branch currents leaving this bus
            sum_cr_out =
                sum(cr_fr[l, t] for l in br_from if pf.branch_status[l] == 1; init = 0.0) +
                sum(cr_to[l, t] for l in br_to if pf.branch_status[l] == 1; init = 0.0)
            sum_ci_out =
                sum(ci_fr[l, t] for l in br_from if pf.branch_status[l] == 1; init = 0.0) +
                sum(ci_to[l, t] for l in br_to if pf.branch_status[l] == 1; init = 0.0)

            # KCL real: Σ crg = Σ cr_branch + (pd·vr + qd·vi)/(vr²+vi²) + gs·vr - bs·vi
            @constraint(model,
                sum_crg ==
                sum_cr_out +
                (pd_total * vr_b + qd_total * vi_b) / (vr_b^2 + vi_b^2) +
                gs_val * vr_b - bs_val * vi_b
            )

            # KCL imag: Σ cig = Σ ci_branch + (pd·vi - qd·vr)/(vr²+vi²) + gs·vi + bs·vr
            @constraint(model,
                sum_cig ==
                sum_ci_out +
                (pd_total * vi_b - qd_total * vr_b) / (vr_b^2 + vi_b^2) +
                gs_val * vi_b + bs_val * vr_b
            )
        end
    end

    return nothing
end
