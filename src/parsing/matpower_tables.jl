const _mp_data_names = [
    "mpc.version", "mpc.baseMVA", "mpc.bus", "mpc.gen",
    "mpc.branch", "mpc.dcline", "mpc.gencost", "mpc.dclinecost",
    "mpc.bus_name", "mpc.storage", "mpc.switch",
    "mpc.bus_timeseries",
    "mpc.gencost_timeseries",
]

# ==============================================================================
# Token helpers
# ------------------------------------------------------------------------------
# All helpers are `@inline`, type-stable, and allocation-free on the hot path.
# `tryparse` handles every numeric case natively (incl. Inf/NaN), so no manual
# fast paths are needed.
# ==============================================================================

@inline _nrows(data::Vector{Vector{String}}) = length(data)

@inline function _tok(row::Vector{String}, j::Int, default::String = "")
    return (j <= length(row)) ? row[j] : default
end

@inline function _strip_token_quotes(s::String)::String
    isempty(s) && return s
    c1 = s[begin]
    if (c1 == '\'' || c1 == '"') && lastindex(s) > firstindex(s) && s[end] == c1
        inner = s[nextind(s, firstindex(s)):prevind(s, lastindex(s))]
        if c1 == '\''
            return replace(inner, "''" => "'")
        else
            return replace(inner, "\"\"" => "\"")
        end
    end
    return s
end

@inline function _as_str(s::String, default::String = "")::String
    isempty(s) && return default
    return _strip_token_quotes(s)
end

@inline function _as_int(s::String, default::Int = 0)::Int
    isempty(s) && return default
    v = tryparse(Int, s)
    return v === nothing ? default : v
end

@inline function _as_f64(s::String, default::Float64 = 0.0)::Float64
    isempty(s) && return default
    v = tryparse(Float64, s)
    return v === nothing ? default : v
end

function _bus_from_mat(data::Vector{Vector{String}})
    n = _nrows(data)

    bus_i = Vector{Int}(undef, n)
    bus_type = Vector{Int}(undef, n)
    pd = Vector{Float64}(undef, n)
    qd = Vector{Float64}(undef, n)
    gs = Vector{Float64}(undef, n)
    bs = Vector{Float64}(undef, n)
    area = Vector{Int}(undef, n)
    vm = Vector{Float64}(undef, n)
    va = Vector{Float64}(undef, n)
    base_kv = Vector{Float64}(undef, n)
    zone = Vector{Int}(undef, n)
    vmax = Vector{Float64}(undef, n)
    vmin = Vector{Float64}(undef, n)
    lam_p = Vector{Float64}(undef, n)
    lam_q = Vector{Float64}(undef, n)
    mu_vmax = Vector{Float64}(undef, n)
    mu_vmin = Vector{Float64}(undef, n)

    @inbounds for i in 1:n
        row = data[i]
        bus_i[i] = _as_int(_tok(row, 1, "0"), 0)
        bus_type[i] = _as_int(_tok(row, 2, "0"), 0)
        pd[i] = _as_f64(_tok(row, 3, ""), 0.0)
        qd[i] = _as_f64(_tok(row, 4, ""), 0.0)
        gs[i] = _as_f64(_tok(row, 5, ""), 0.0)
        bs[i] = _as_f64(_tok(row, 6, ""), 0.0)
        area[i] = _as_int(_tok(row, 7, "0"), 0)
        vm[i] = _as_f64(_tok(row, 8, ""), 0.0)
        va[i] = _as_f64(_tok(row, 9, ""), 0.0)
        base_kv[i] = _as_f64(_tok(row, 10, ""), 0.0)
        zone[i] = _as_int(_tok(row, 11, "0"), 0)
        vmax[i] = _as_f64(_tok(row, 12, ""), 0.0)
        vmin[i] = _as_f64(_tok(row, 13, ""), 0.0)
        lam_p[i] = _as_f64(_tok(row, 14, ""), 0.0)
        lam_q[i] = _as_f64(_tok(row, 15, ""), 0.0)
        mu_vmax[i] = _as_f64(_tok(row, 16, ""), 0.0)
        mu_vmin[i] = _as_f64(_tok(row, 17, ""), 0.0)
    end

    return BusRaw{Float64}(
        bus_i, bus_type,
        pd, qd, gs, bs,
        area, vm, va, base_kv, zone,
        vmax, vmin,
        lam_p, lam_q,
        mu_vmax, mu_vmin,
        falses(n),
    )
end

function _busname_from_mat(data::Vector{Vector{String}})
    n = _nrows(data)
    names = Vector{String}(undef, n)
    @inbounds for i in 1:n
        row = data[i]
        names[i] = _as_str(_tok(row, 1, ""), "")
    end
    return BusNameRaw(names)
end

function _gen_from_mat(data::Vector{Vector{String}})
    n = _nrows(data)

    gen_bus = Vector{Int}(undef, n)
    pg = Vector{Float64}(undef, n)
    qg = Vector{Float64}(undef, n)
    qmax = Vector{Float64}(undef, n)
    qmin = Vector{Float64}(undef, n)
    vg = Vector{Float64}(undef, n)
    mbase = Vector{Float64}(undef, n)
    status = Vector{Int}(undef, n)
    pmax = Vector{Float64}(undef, n)
    pmin = Vector{Float64}(undef, n)
    pc1 = Vector{Float64}(undef, n)
    pc2 = Vector{Float64}(undef, n)
    qc1min = Vector{Float64}(undef, n)
    qc1max = Vector{Float64}(undef, n)
    qc2min = Vector{Float64}(undef, n)
    qc2max = Vector{Float64}(undef, n)
    ramp_agc = Vector{Float64}(undef, n)
    ramp_10 = Vector{Float64}(undef, n)
    ramp_30 = Vector{Float64}(undef, n)
    ramp_q = Vector{Float64}(undef, n)
    apf = Vector{Float64}(undef, n)

    mu_pmax = Vector{Float64}(undef, n)
    mu_pmin = Vector{Float64}(undef, n)
    mu_qmax = Vector{Float64}(undef, n)
    mu_qmin = Vector{Float64}(undef, n)

    plant_type = Vector{String}(undef, n)
    curt_cost = Vector{Float64}(undef, n)

    @inbounds for i in 1:n
        row = data[i]
        gen_bus[i] = _as_int(_tok(row, 1, "0"), 0)
        pg[i] = _as_f64(_tok(row, 2, ""), 0.0)
        qg[i] = _as_f64(_tok(row, 3, ""), 0.0)
        qmax[i] = _as_f64(_tok(row, 4, ""), 0.0)
        qmin[i] = _as_f64(_tok(row, 5, ""), 0.0)
        vg[i] = _as_f64(_tok(row, 6, ""), 0.0)
        mbase[i] = _as_f64(_tok(row, 7, ""), 0.0)
        status[i] = _as_int(_tok(row, 8, "0"), 0)
        pmax[i] = _as_f64(_tok(row, 9, ""), 0.0)
        pmin[i] = _as_f64(_tok(row, 10, ""), 0.0)
        pc1[i] = _as_f64(_tok(row, 11, ""), 0.0)
        pc2[i] = _as_f64(_tok(row, 12, ""), 0.0)
        qc1min[i] = _as_f64(_tok(row, 13, ""), 0.0)
        qc1max[i] = _as_f64(_tok(row, 14, ""), 0.0)
        qc2min[i] = _as_f64(_tok(row, 15, ""), 0.0)
        qc2max[i] = _as_f64(_tok(row, 16, ""), 0.0)
        ramp_agc[i] = _as_f64(_tok(row, 17, ""), 0.0)
        ramp_10[i] = _as_f64(_tok(row, 18, ""), 0.0)
        ramp_30[i] = _as_f64(_tok(row, 19, ""), 0.0)
        ramp_q[i] = _as_f64(_tok(row, 20, ""), 0.0)
        apf[i] = _as_f64(_tok(row, 21, ""), 0.0)

        mu_pmax[i] = _as_f64(_tok(row, 22, ""), 0.0)
        mu_pmin[i] = _as_f64(_tok(row, 23, ""), 0.0)
        mu_qmax[i] = _as_f64(_tok(row, 24, ""), 0.0)
        mu_qmin[i] = _as_f64(_tok(row, 25, ""), 0.0)

        plant_type[i] = _as_str(_tok(row, 26, "'T'"), "T")
        curt_cost[i] = _as_f64(_tok(row, 27, ""), 0.0)
    end

    return GenRaw{Float64}(
        gen_bus, pg, qg, qmax, qmin, vg, mbase, status, pmax, pmin,
        pc1, pc2, qc1min, qc1max, qc2min, qc2max,
        ramp_agc, ramp_10, ramp_30, ramp_q, apf,
        mu_pmax, mu_pmin, mu_qmax, mu_qmin,
        plant_type, curt_cost,
    )
end

function _branch_from_mat(data::Vector{Vector{String}})
    n = _nrows(data)

    f_bus = Vector{Int}(undef, n)
    t_bus = Vector{Int}(undef, n)
    br_r = Vector{Float64}(undef, n)
    br_x = Vector{Float64}(undef, n)
    br_b = Vector{Float64}(undef, n)
    rate_a = Vector{Float64}(undef, n)
    rate_b = Vector{Float64}(undef, n)
    rate_c = Vector{Float64}(undef, n)
    tap = Vector{Float64}(undef, n)
    shift = Vector{Float64}(undef, n)
    status = Vector{Int}(undef, n)
    angmin = Vector{Float64}(undef, n)
    angmax = Vector{Float64}(undef, n)
    pf = Vector{Float64}(undef, n)
    qf = Vector{Float64}(undef, n)
    pt = Vector{Float64}(undef, n)
    qt = Vector{Float64}(undef, n)
    mu_sf = Vector{Float64}(undef, n)
    mu_st = Vector{Float64}(undef, n)
    mu_angmin = Vector{Float64}(undef, n)
    mu_angmax = Vector{Float64}(undef, n)

    @inbounds for i in 1:n
        row = data[i]
        f_bus[i] = _as_int(_tok(row, 1, "0"), 0)
        t_bus[i] = _as_int(_tok(row, 2, "0"), 0)
        br_r[i] = _as_f64(_tok(row, 3, ""), 0.0)
        br_x[i] = _as_f64(_tok(row, 4, ""), 0.0)
        br_b[i] = _as_f64(_tok(row, 5, ""), 0.0)
        rate_a[i] = _as_f64(_tok(row, 6, ""), 0.0)
        rate_b[i] = _as_f64(_tok(row, 7, ""), 0.0)
        rate_c[i] = _as_f64(_tok(row, 8, ""), 0.0)
        tap[i] = _as_f64(_tok(row, 9, ""), 0.0)
        shift[i] = _as_f64(_tok(row, 10, ""), 0.0)
        status[i] = _as_int(_tok(row, 11, "0"), 0)

        angmin[i] = _as_f64(_tok(row, 12, ""), 0.0)
        angmax[i] = _as_f64(_tok(row, 13, ""), 0.0)

        pf[i] = _as_f64(_tok(row, 14, ""), 0.0)
        qf[i] = _as_f64(_tok(row, 15, ""), 0.0)
        pt[i] = _as_f64(_tok(row, 16, ""), 0.0)
        qt[i] = _as_f64(_tok(row, 17, ""), 0.0)

        mu_sf[i] = _as_f64(_tok(row, 18, ""), 0.0)
        mu_st[i] = _as_f64(_tok(row, 19, ""), 0.0)
        mu_angmin[i] = _as_f64(_tok(row, 20, ""), 0.0)
        mu_angmax[i] = _as_f64(_tok(row, 21, ""), 0.0)
    end

    return BranchRaw{Float64}(
        f_bus, t_bus,
        br_r, br_x, br_b,
        rate_a, rate_b, rate_c,
        tap, shift,
        status,
        angmin, angmax,
        pf, qf, pt, qt,
        mu_sf, mu_st, mu_angmin, mu_angmax,
    )
end

function _dcline_from_mat(data::Vector{Vector{String}})
    n = _nrows(data)

    f_bus = Vector{Int}(undef, n)
    t_bus = Vector{Int}(undef, n)
    br_status = Vector{Int}(undef, n)
    pf = Vector{Float64}(undef, n)
    pt = Vector{Float64}(undef, n)
    qf = Vector{Float64}(undef, n)
    qt = Vector{Float64}(undef, n)
    vf = Vector{Float64}(undef, n)
    vt = Vector{Float64}(undef, n)
    pmin = Vector{Float64}(undef, n)
    pmax = Vector{Float64}(undef, n)
    qminf = Vector{Float64}(undef, n)
    qmaxf = Vector{Float64}(undef, n)
    qmint = Vector{Float64}(undef, n)
    qmaxt = Vector{Float64}(undef, n)
    loss0 = Vector{Float64}(undef, n)
    loss1 = Vector{Float64}(undef, n)
    mu_pmin = Vector{Float64}(undef, n)
    mu_pmax = Vector{Float64}(undef, n)
    mu_qminf = Vector{Float64}(undef, n)
    mu_qmaxf = Vector{Float64}(undef, n)
    mu_qmint = Vector{Float64}(undef, n)
    mu_qmaxt = Vector{Float64}(undef, n)

    @inbounds for i in 1:n
        row = data[i]
        f_bus[i] = _as_int(_tok(row, 1, "0"), 0)
        t_bus[i] = _as_int(_tok(row, 2, "0"), 0)
        br_status[i] = _as_int(_tok(row, 3, "0"), 0)
        pf[i] = _as_f64(_tok(row, 4, ""), 0.0)
        pt[i] = _as_f64(_tok(row, 5, ""), 0.0)
        qf[i] = _as_f64(_tok(row, 6, ""), 0.0)
        qt[i] = _as_f64(_tok(row, 7, ""), 0.0)
        vf[i] = _as_f64(_tok(row, 8, ""), 0.0)
        vt[i] = _as_f64(_tok(row, 9, ""), 0.0)
        pmin[i] = _as_f64(_tok(row, 10, ""), 0.0)
        pmax[i] = _as_f64(_tok(row, 11, ""), 0.0)
        qminf[i] = _as_f64(_tok(row, 12, ""), 0.0)
        qmaxf[i] = _as_f64(_tok(row, 13, ""), 0.0)
        qmint[i] = _as_f64(_tok(row, 14, ""), 0.0)
        qmaxt[i] = _as_f64(_tok(row, 15, ""), 0.0)
        loss0[i] = _as_f64(_tok(row, 16, ""), 0.0)
        loss1[i] = _as_f64(_tok(row, 17, ""), 0.0)
        mu_pmin[i] = _as_f64(_tok(row, 18, ""), 0.0)
        mu_pmax[i] = _as_f64(_tok(row, 19, ""), 0.0)
        mu_qminf[i] = _as_f64(_tok(row, 20, ""), 0.0)
        mu_qmaxf[i] = _as_f64(_tok(row, 21, ""), 0.0)
        mu_qmint[i] = _as_f64(_tok(row, 22, ""), 0.0)
        mu_qmaxt[i] = _as_f64(_tok(row, 23, ""), 0.0)
    end

    return DCLineRaw{Float64}(
        f_bus, t_bus, br_status,
        pf, pt, qf, qt, vf, vt,
        pmin, pmax, qminf, qmaxf, qmint, qmaxt,
        loss0, loss1,
        mu_pmin, mu_pmax, mu_qminf, mu_qmaxf, mu_qmint, mu_qmaxt,
    )
end

function _storage_from_mat(data::Vector{Vector{String}})
    n = _nrows(data)

    storage_bus = Vector{Int}(undef, n)
    ps = Vector{Float64}(undef, n)
    qs = Vector{Float64}(undef, n)
    energy = Vector{Float64}(undef, n)
    energy_rating = Vector{Float64}(undef, n)
    charge_rating = Vector{Float64}(undef, n)
    discharge_rating = Vector{Float64}(undef, n)
    charge_efficiency = Vector{Float64}(undef, n)
    discharge_efficiency = Vector{Float64}(undef, n)
    thermal_rating = Vector{Float64}(undef, n)
    qmin = Vector{Float64}(undef, n)
    qmax = Vector{Float64}(undef, n)
    r = Vector{Float64}(undef, n)
    x = Vector{Float64}(undef, n)
    p_loss = Vector{Float64}(undef, n)
    q_loss = Vector{Float64}(undef, n)
    status = Vector{Int}(undef, n)

    @inbounds for i in 1:n
        row = data[i]
        storage_bus[i] = _as_int(_tok(row, 1, "0"), 0)
        ps[i] = _as_f64(_tok(row, 2, ""), 0.0)
        qs[i] = _as_f64(_tok(row, 3, ""), 0.0)
        energy[i] = _as_f64(_tok(row, 4, ""), 0.0)
        energy_rating[i] = _as_f64(_tok(row, 5, ""), 0.0)
        charge_rating[i] = _as_f64(_tok(row, 6, ""), 0.0)
        discharge_rating[i] = _as_f64(_tok(row, 7, ""), 0.0)
        charge_efficiency[i] = _as_f64(_tok(row, 8, ""), 0.0)
        discharge_efficiency[i] = _as_f64(_tok(row, 9, ""), 0.0)
        thermal_rating[i] = _as_f64(_tok(row, 10, ""), 0.0)
        qmin[i] = _as_f64(_tok(row, 11, ""), 0.0)
        qmax[i] = _as_f64(_tok(row, 12, ""), 0.0)
        r[i] = _as_f64(_tok(row, 13, ""), 0.0)
        x[i] = _as_f64(_tok(row, 14, ""), 0.0)
        p_loss[i] = _as_f64(_tok(row, 15, ""), 0.0)
        q_loss[i] = _as_f64(_tok(row, 16, ""), 0.0)
        status[i] = _as_int(_tok(row, 17, "0"), 0)
    end

    return StorageRaw{Float64}(
        storage_bus, ps, qs, energy, energy_rating,
        charge_rating, discharge_rating, charge_efficiency, discharge_efficiency,
        thermal_rating, qmin, qmax, r, x, p_loss, q_loss, status,
    )
end

function _switch_from_mat(data::Vector{Vector{String}})
    n = _nrows(data)

    f_bus = Vector{Int}(undef, n)
    t_bus = Vector{Int}(undef, n)
    psw = Vector{Float64}(undef, n)
    qsw = Vector{Float64}(undef, n)
    state = Vector{Int}(undef, n)
    thermal_rating = Vector{Float64}(undef, n)
    status = Vector{Int}(undef, n)

    @inbounds for i in 1:n
        row = data[i]
        f_bus[i] = _as_int(_tok(row, 1, "0"), 0)
        t_bus[i] = _as_int(_tok(row, 2, "0"), 0)
        psw[i] = _as_f64(_tok(row, 3, ""), 0.0)
        qsw[i] = _as_f64(_tok(row, 4, ""), 0.0)
        state[i] = _as_int(_tok(row, 5, "0"), 0)
        thermal_rating[i] = _as_f64(_tok(row, 6, ""), 0.0)
        status[i] = _as_int(_tok(row, 7, "0"), 0)
    end

    return SwitchRaw{Float64}(f_bus, t_bus, psw, qsw, state, thermal_rating, status)
end

function _bus_ts_from_mat(data::Vector{Vector{String}})
    n = _nrows(data)

    time = Vector{Int}(undef, n)
    bus = Vector{Int}(undef, n)
    pd = Vector{Float64}(undef, n)
    qd = Vector{Float64}(undef, n)

    @inbounds for i in 1:n
        row = data[i]
        time[i] = _as_int(_tok(row, 1, "0"), 0)
        bus[i] = _as_int(_tok(row, 2, "0"), 0)
        pd[i] = _as_f64(_tok(row, 3, ""), 0.0)
        qd[i] = _as_f64(_tok(row, 4, ""), 0.0)
    end

    return BusTimeseriesRaw{Float64}(time, bus, pd, qd, fill(0.0, n), fill(0.0, n))
end

function _gen_ts_from_mat(data::Vector{Vector{String}})
    n = _nrows(data)
    time = Vector{Int}(undef, n)
    gen_id = Vector{Int}(undef, n)
    pmax = Vector{Float64}(undef, n)

    @inbounds for i in 1:n
        row = data[i]
        time[i] = _as_int(_tok(row, 1, "0"), 0)
        gen_id[i] = _as_int(_tok(row, 2, "0"), 0)
        pmax[i] = _as_f64(_tok(row, 3, ""), 0.0)
    end

    return GenTimeseriesRaw{Float64}(time, gen_id, pmax)
end

function _parse_gencost_table_raw(data::Vector{Vector{String}})
    n = _nrows(data)
    model = Vector{Int}(undef, n)
    startup = Vector{Float64}(undef, n)
    shutdown = Vector{Float64}(undef, n)
    ncost = Vector{Int}(undef, n)
    cost_coefficients = Vector{Vector{Float64}}(undef, n)

    @inbounds for i in 1:n
        row = data[i]
        model[i] = _as_int(_tok(row, 1, "0"), 0)
        startup[i] = _as_f64(_tok(row, 2, ""), 0.0)
        shutdown[i] = _as_f64(_tok(row, 3, ""), 0.0)
        ncost[i] = _as_int(_tok(row, 4, "0"), 0)

        k = ncost[i]
        if k > 0 && length(row) >= 4 + k
            coeff = Vector{Float64}(undef, k)
            @inbounds for j in 1:k
                coeff[j] = _as_f64(row[4+j], 0.0)
            end
            cost_coefficients[i] = coeff
        else
            cost_coefficients[i] = Float64[]
        end
    end

    return GencostRaw{Float64}(model, startup, shutdown, ncost, cost_coefficients)
end

function _parse_gencost_timeseries_table_raw(data::Vector{Vector{String}})
    n = _nrows(data)
    time = Vector{Int}(undef, n)
    gen_id = Vector{Int}(undef, n)
    model = Vector{Int}(undef, n)
    startup = Vector{Float64}(undef, n)
    shutdown = Vector{Float64}(undef, n)
    ncost = Vector{Int}(undef, n)
    cost_coefficients = Vector{Vector{Float64}}(undef, n)
    tou = Vector{Float64}(undef, n)
    feed_in_tariff = Vector{Float64}(undef, n)

    @inbounds for i in 1:n
        row = data[i]
        time[i] = _as_int(_tok(row, 1, "0"), 0)
        gen_id[i] = _as_int(_tok(row, 2, "0"), 0)
        model[i] = _as_int(_tok(row, 3, "0"), 0)
        startup[i] = _as_f64(_tok(row, 4, ""), 0.0)
        shutdown[i] = _as_f64(_tok(row, 5, ""), 0.0)
        ncost[i] = _as_int(_tok(row, 6, "0"), 0)

        k = ncost[i]
        if k > 0 && length(row) >= 6 + k
            coeff = Vector{Float64}(undef, k)
            @inbounds for j in 1:k
                coeff[j] = _as_f64(row[6+j], 0.0)
            end
            cost_coefficients[i] = coeff
        else
            cost_coefficients[i] = Float64[]
        end

        tou[i] = _as_f64(_tok(row, 7 + k, ""), 0.0)
        feed_in_tariff[i] = _as_f64(_tok(row, 8 + k, ""), 0.0)
    end

    return GencostTimeseriesRaw{Float64}(
        time, gen_id, model, startup, shutdown, ncost,
        cost_coefficients, tou, feed_in_tariff,
    )
end
