module CSVExt

using Lirion
using CSV
using DataFrames

# ============================================================
# DataFrame helpers
# ============================================================

@inline _sym(s::AbstractString) = Symbol(String(s))

function _read_csv(path::String)::DataFrame
    return CSV.read(path, DataFrame; normalizenames = false, silencewarnings = true)
end

@inline _hascol(df::DataFrame, col::AbstractString)::Bool =
    _sym(col) in propertynames(df)

function _icol(df::DataFrame, col::AbstractString, default::Int)::Vector{Int}
    s = _sym(col)
    n = nrow(df)
    out = Vector{Int}(undef, n)
    if s in propertynames(df)
        v = df[!, s]
        @inbounds for i in 1:n
            x = v[i]
            out[i] = (x === missing || x === nothing) ? default : Int(x)
        end
    else
        fill!(out, default)
    end
    return out
end

function _fcol(df::DataFrame, col::AbstractString, default::Float64)::Vector{Float64}
    s = _sym(col)
    n = nrow(df)
    out = Vector{Float64}(undef, n)
    if s in propertynames(df)
        v = df[!, s]
        @inbounds for i in 1:n
            x = v[i]
            out[i] = (x === missing || x === nothing) ? default : Float64(x)
        end
    else
        fill!(out, default)
    end
    return out
end

function _scol(df::DataFrame, col::AbstractString, default::String)::Vector{String}
    s = _sym(col)
    n = nrow(df)
    out = Vector{String}(undef, n)
    if s in propertynames(df)
        v = df[!, s]
        @inbounds for i in 1:n
            x = v[i]
            out[i] = (x === missing || x === nothing || (x isa AbstractString && isempty(x))) ?
                     default : String(x)
        end
    else
        fill!(out, default)
    end
    return out
end

# ============================================================
# Table → Raw struct converters
# ============================================================

function _bus_from_df(df::DataFrame)
    n = nrow(df)
    return Lirion.BusRaw{Float64}(
        _icol(df, "bus_i",    0),
        _icol(df, "bus_type", 0),
        _fcol(df, "pd",       0.0),
        _fcol(df, "qd",       0.0),
        _fcol(df, "gs",       0.0),
        _fcol(df, "bs",       0.0),
        _icol(df, "area",     0),
        _fcol(df, "vm",       1.0),
        _fcol(df, "va",       0.0),
        _fcol(df, "base_kv",  0.0),
        _icol(df, "zone",     0),
        _fcol(df, "vmax",     1.1),
        _fcol(df, "vmin",     0.9),
        _fcol(df, "lam_p",    0.0),
        _fcol(df, "lam_q",    0.0),
        _fcol(df, "mu_vmax",  0.0),
        _fcol(df, "mu_vmin",  0.0),
        falses(n),
    )
end

function _gen_from_df(df::DataFrame)
    n = nrow(df)

    plant_type = _scol(df, "plant_type", "T")
    @inbounds for i in 1:n
        isempty(plant_type[i]) && (plant_type[i] = "T")
    end

    return Lirion.GenRaw{Float64}(
        _icol(df, "gen_bus",          0),
        _fcol(df, "pg",               0.0),
        _fcol(df, "qg",               0.0),
        _fcol(df, "qmax",             0.0),
        _fcol(df, "qmin",             0.0),
        _fcol(df, "vg",               1.0),
        _fcol(df, "mbase",            100.0),
        _icol(df, "gen_status",       1),
        _fcol(df, "pmax",             0.0),
        _fcol(df, "pmin",             0.0),
        _fcol(df, "pc1",              0.0),
        _fcol(df, "pc2",              0.0),
        _fcol(df, "qc1min",           0.0),
        _fcol(df, "qc1max",           0.0),
        _fcol(df, "qc2min",           0.0),
        _fcol(df, "qc2max",           0.0),
        _fcol(df, "ramp_agc",         0.0),
        _fcol(df, "ramp_10",          0.0),
        _fcol(df, "ramp_30",          0.0),
        _fcol(df, "ramp_q",           0.0),
        _fcol(df, "apf",              0.0),
        _fcol(df, "mu_pmax",          0.0),
        _fcol(df, "mu_pmin",          0.0),
        _fcol(df, "mu_qmax",          0.0),
        _fcol(df, "mu_qmin",          0.0),
        plant_type,
        _fcol(df, "curtailment_cost", 0.0),
    )
end

function _branch_from_df(df::DataFrame)
    return Lirion.BranchRaw{Float64}(
        _icol(df, "f_bus",     0),
        _icol(df, "t_bus",     0),
        _fcol(df, "br_r",      0.0),
        _fcol(df, "br_x",      0.0),
        _fcol(df, "br_b",      0.0),
        _fcol(df, "rate_a",    0.0),
        _fcol(df, "rate_b",    0.0),
        _fcol(df, "rate_c",    0.0),
        _fcol(df, "tap",       0.0),
        _fcol(df, "shift",     0.0),
        _icol(df, "br_status", 1),
        _fcol(df, "angmin",    -60.0),
        _fcol(df, "angmax",     60.0),
        _fcol(df, "pf",        0.0),
        _fcol(df, "qf",        0.0),
        _fcol(df, "pt",        0.0),
        _fcol(df, "qt",        0.0),
        _fcol(df, "mu_sf",     0.0),
        _fcol(df, "mu_st",     0.0),
        _fcol(df, "mu_angmin", 0.0),
        _fcol(df, "mu_angmax", 0.0),
    )
end

function _dcline_from_df(df::DataFrame)
    return Lirion.DCLineRaw{Float64}(
        _icol(df, "f_bus",     0),
        _icol(df, "t_bus",     0),
        _icol(df, "br_status", 1),
        _fcol(df, "pf",        0.0),
        _fcol(df, "pt",        0.0),
        _fcol(df, "qf",        0.0),
        _fcol(df, "qt",        0.0),
        _fcol(df, "vf",        1.0),
        _fcol(df, "vt",        1.0),
        _fcol(df, "pmin",      0.0),
        _fcol(df, "pmax",      0.0),
        _fcol(df, "qminf",     0.0),
        _fcol(df, "qmaxf",     0.0),
        _fcol(df, "qmint",     0.0),
        _fcol(df, "qmaxt",     0.0),
        _fcol(df, "loss0",     0.0),
        _fcol(df, "loss1",     0.0),
        _fcol(df, "mu_pmin",   0.0),
        _fcol(df, "mu_pmax",   0.0),
        _fcol(df, "mu_qminf",  0.0),
        _fcol(df, "mu_qmaxf",  0.0),
        _fcol(df, "mu_qmint",  0.0),
        _fcol(df, "mu_qmaxt",  0.0),
    )
end

function _storage_from_df(df::DataFrame)
    return Lirion.StorageRaw{Float64}(
        _icol(df, "storage_bus",         0),
        _fcol(df, "ps",                  0.0),
        _fcol(df, "qs",                  0.0),
        _fcol(df, "energy",              0.0),
        _fcol(df, "energy_rating",       0.0),
        _fcol(df, "charge_rating",       0.0),
        _fcol(df, "discharge_rating",    0.0),
        _fcol(df, "charge_efficiency",   1.0),
        _fcol(df, "discharge_efficiency",1.0),
        _fcol(df, "thermal_rating",      0.0),
        _fcol(df, "qmin",                0.0),
        _fcol(df, "qmax",                0.0),
        _fcol(df, "r",                   0.0),
        _fcol(df, "x",                   0.0),
        _fcol(df, "p_loss",              0.0),
        _fcol(df, "q_loss",              0.0),
        _icol(df, "status",              1),
    )
end

function _switch_from_df(df::DataFrame)
    return Lirion.SwitchRaw{Float64}(
        _icol(df, "f_bus",          0),
        _icol(df, "t_bus",          0),
        _fcol(df, "psw",            0.0),
        _fcol(df, "qsw",            0.0),
        _icol(df, "state",          1),
        _fcol(df, "thermal_rating", 0.0),
        _icol(df, "status",         1),
    )
end

function _bus_ts_from_df(df::DataFrame)
    return Lirion.BusTimeseriesRaw{Float64}(
        _icol(df, "time", 1),
        _icol(df, "bus",  0),
        _fcol(df, "pd",   0.0),
        _fcol(df, "qd",   0.0),
        _fcol(df, "gs",   0.0),
        _fcol(df, "bs",   0.0),
    )
end

function _gen_ts_from_df(df::DataFrame)
    return Lirion.GenTimeseriesRaw{Float64}(
        _icol(df, "time", 1),
        _icol(df, "gen_id", 0),
        _fcol(df, "pmax", NaN),
        _fcol(df, "pmin", NaN),
        _fcol(df, "qmax", NaN),
        _fcol(df, "qmin", NaN),
    )
end

function _gencost_from_df(df::DataFrame)
    n = nrow(df)

    model    = _icol(df, "model",    0)
    startup  = _fcol(df, "startup",  0.0)
    shutdown = _fcol(df, "shutdown", 0.0)
    ncost    = _icol(df, "ncost",    0)
    coeffs   = Vector{Vector{Float64}}(undef, n)

    max_k = isempty(ncost) ? 0 : maximum(ncost)
    props = propertynames(df)
    coef_cols = Vector{Union{Vector{Float64}, Nothing}}(undef, max_k)
    @inbounds for j in 1:max_k
        sj = _sym("c$(j)")
        coef_cols[j] = sj in props ? _fcol(df, "c$(j)", 0.0) : nothing
    end

    @inbounds for i in 1:n
        k = ncost[i]
        if k <= 0
            coeffs[i] = Float64[]
            continue
        end
        c = Vector{Float64}(undef, k)
        @inbounds for j in 1:k
            col = coef_cols[j]
            c[j] = col === nothing ? 0.0 : col[i]
        end
        coeffs[i] = c
    end

    return Lirion.GencostRaw{Float64}(model, startup, shutdown, ncost, coeffs)
end

function _gencost_ts_from_df(df::DataFrame)
    n = nrow(df)

    ncost  = _icol(df, "ncost", 0)
    coeffs = Vector{Vector{Float64}}(undef, n)

    max_k = isempty(ncost) ? 0 : maximum(ncost)
    props = propertynames(df)
    coef_cols = Vector{Union{Vector{Float64}, Nothing}}(undef, max_k)
    @inbounds for j in 1:max_k
        sj = _sym("c$(j)")
        coef_cols[j] = sj in props ? _fcol(df, "c$(j)", 0.0) : nothing
    end

    @inbounds for i in 1:n
        k = ncost[i]
        if k <= 0
            coeffs[i] = Float64[]
            continue
        end
        c = Vector{Float64}(undef, k)
        @inbounds for j in 1:k
            col = coef_cols[j]
            c[j] = col === nothing ? 0.0 : col[i]
        end
        coeffs[i] = c
    end

    stou = _sym("tou");          has_tou = stou in propertynames(df)
    sfit = _sym("feed_in_tariff"); has_fit = sfit in propertynames(df)
    tou = has_tou ? _fcol(df, "tou",           0.0) : zeros(n)
    fit = has_fit ? _fcol(df, "feed_in_tariff", 0.0) : zeros(n)

    return Lirion.GencostTimeseriesRaw{Float64}(
        _icol(df, "time",     1),
        _icol(df, "gen_id",   0),
        _icol(df, "model",    0),
        _fcol(df, "startup",  0.0),
        _fcol(df, "shutdown", 0.0),
        ncost, coeffs, tou, fit,
    )
end

# ============================================================
# _csv_parse_case — extension hook registered in case_parser.jl
# ============================================================

function Lirion._csv_parse_case(dir::String)::Lirion.CaseRaw
    d = abspath(String(dir))

    case = Lirion.CaseRaw(;
        name           = splitdir(d)[end],
        source_type    = "csv",
        source_version = "2",
        baseMVA        = 100.0,
        T              = 1,
    )

    let p = joinpath(d, "baseMVA.csv")
        if isfile(p)
            df = _read_csv(p)
            if nrow(df) > 0 && _hascol(df, "baseMVA")
                x = df[1, _sym("baseMVA")]
                case.baseMVA = (x === missing) ? 100.0 : Float64(x)
            else
                @warn "baseMVA.csv empty or missing baseMVA column. Using 100.0"
            end
        else
            @warn "baseMVA.csv not found. Using 100.0"
        end
    end

    let p = joinpath(d, "bus.csv")
        isfile(p) && (case.bus = _bus_from_df(_read_csv(p)))
    end

    let p = joinpath(d, "gen.csv")
        isfile(p) && (case.gen = _gen_from_df(_read_csv(p)))
    end

    let p = joinpath(d, "branch.csv")
        isfile(p) && (case.branch = _branch_from_df(_read_csv(p)))
    end

    let p = joinpath(d, "dcline.csv")
        isfile(p) && (case.dcline = _dcline_from_df(_read_csv(p)))
    end

    let p = joinpath(d, "storage.csv")
        isfile(p) && (case.storage = _storage_from_df(_read_csv(p)))
    end

    let p = joinpath(d, "switch.csv")
        isfile(p) && (case.switch = _switch_from_df(_read_csv(p)))
    end

    let p = joinpath(d, "bus_name.csv")
        if isfile(p)
            df = _read_csv(p)
            case.bus_name = Lirion.BusNameRaw(_scol(df, "name", ""))
            if case.bus !== nothing && length(case.bus_name.name) != length(case.bus.bus_i)
                @error "Incorrect CSV case: number of bus names != number of buses."
            end
        end
    end

    let p = joinpath(d, "gencost.csv")
        isfile(p) && (case.gencost = _gencost_from_df(_read_csv(p)))
    end

    tdir = joinpath(d, "time")
    if isdir(tdir)
        let p = joinpath(tdir, "bus_timeseries.csv")
            if isfile(p)
                case.bus_timeseries = _bus_ts_from_df(_read_csv(p))
                !isempty(case.bus_timeseries.time) &&
                    (case.T = maximum(case.bus_timeseries.time))
            end
        end

        let p = joinpath(tdir, "gencost_timeseries.csv")
            if isfile(p)
                case.gencost_timeseries = _gencost_ts_from_df(_read_csv(p))
                !isempty(case.gencost_timeseries.time) &&
                    (case.T = max(case.T, maximum(case.gencost_timeseries.time)))
            end
        end

        let p = joinpath(tdir, "gen_timeseries.csv")
            if isfile(p)
                case.gen_timeseries = _gen_ts_from_df(_read_csv(p))
                !isempty(case.gen_timeseries.time) &&
                    (case.T = max(case.T, maximum(case.gen_timeseries.time)))
            end
        end
    end

    return case
end

# ============================================================
# LirionResult → CSV dump
# ============================================================

function Lirion.to_csvs(result::Lirion.LirionResult, outdir::AbstractString)
    mkpath(outdir)

    summary_df = DataFrame(
        algorithm = [String(result.algorithm)],
        base_model = [String(result.base_model)],
        coupling = [String(result.coupling)],
        sectioning = [String(result.sectioning)],
        objective = [result.objective],
        status = [String(result.status)],
        solve_time = [result.solve_time],
        iter_done = [result.iter_done === nothing ? -1 : result.iter_done],
    )
    CSV.write(joinpath(outdir, "summary.csv"), summary_df)

    nrow(result.bus) > 0 && CSV.write(joinpath(outdir, "bus.csv"), result.bus)
    nrow(result.gen) > 0 && CSV.write(joinpath(outdir, "gen.csv"), result.gen)
    nrow(result.branch) > 0 && CSV.write(joinpath(outdir, "branch.csv"), result.branch)

    if result.consensus !== nothing && nrow(result.consensus) > 0
        CSV.write(joinpath(outdir, "consensus.csv"), result.consensus)
    end
    if result.zones !== nothing
        for z in result.zones
            zdir = joinpath(outdir, "zone_$(z.id)")
            mkpath(zdir)
            nrow(z.bus) > 0 && CSV.write(joinpath(zdir, "bus.csv"), z.bus)
            nrow(z.gen) > 0 && CSV.write(joinpath(zdir, "gen.csv"), z.gen)
            nrow(z.branch) > 0 && CSV.write(joinpath(zdir, "branch.csv"), z.branch)
        end
    end

    return outdir
end

end 
