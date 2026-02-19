# -------------------------
# helpers
# -------------------------

@inline _sym(s::AbstractString) = Symbol(String(s))

function _read_csv(path::String)::DataFrame
    return CSV.read(path, DataFrame; normalizenames = false, silencewarnings = true)
end

@inline function _hascol(df::DataFrame, col::AbstractString)::Bool
    return _sym(col) in propertynames(df)
end

@inline function _nrows_df(df::DataFrame)::Int
    return nrow(df)
end

@inline function _getcol(df::DataFrame, col::AbstractString)
    s = _sym(col)
    return (s in propertynames(df)) ? df[!, s] : nothing
end

function _col_or_default(df::DataFrame, col::AbstractString, ::Type{Int}, default::Int)
    v = _getcol(df, col)
    n = _nrows_df(df)
    out = Vector{Int}(undef, n)

    if v === nothing
        fill!(out, default)
        return out
    end

    @inbounds for i in 1:n
        x = v[i]
        out[i] = (x === missing) ? default : Int(x)
    end
    return out
end

function _col_or_default(
    df::DataFrame,
    col::AbstractString,
    ::Type{Float64},
    default::Float64,
)
    v = _getcol(df, col)
    n = _nrows_df(df)
    out = Vector{Float64}(undef, n)

    if v === nothing
        fill!(out, default)
        return out
    end

    @inbounds for i in 1:n
        x = v[i]
        out[i] = (x === missing) ? default : Float64(x)
    end
    return out
end

function _col_or_default(
    df::DataFrame,
    col::AbstractString,
    ::Type{String},
    default::String,
)
    v = _getcol(df, col)
    n = _nrows_df(df)
    out = Vector{String}(undef, n)

    if v === nothing
        fill!(out, default)
        return out
    end

    @inbounds for i in 1:n
        x = v[i]
        out[i] = (x === missing) ? default : String(x)
    end
    return out
end

function _col_or_default(
    df::DataFrame,
    col::AbstractString,
    ::Type{Union{String, SubString{String}}},
    default::String,
)
    return _col_or_default(df, col, String, default)
end

function _parse_csv_table_dict(path::String, columns::Vector{Tuple{String, Type}})
    df = _read_csv(path)
    tbl = Dict{String, Any}()
    for (name, typ) in columns
        if typ == Int
            tbl[name] = _col_or_default(df, name, Int, 0)
        elseif typ == Float64
            tbl[name] = _col_or_default(df, name, Float64, 0.0)
        elseif typ == String || typ == Union{String, SubString{String}}
            tbl[name] = _col_or_default(df, name, String, "")
        else
            tbl[name] = _col_or_default(df, name, String, "")
        end
    end
    return tbl
end

# -------------------------
# gencost + gencost_timeseries (colunas c1,c2,...)
# -------------------------

@inline function _safe_int(x, default::Int = 0)::Int
    return (x === missing || x === nothing) ? default : Int(x)
end

@inline function _safe_f64(x, default::Float64 = 0.0)::Float64
    return (x === missing || x === nothing) ? default : Float64(x)
end

function _parse_gencost_csv(path::String)::GencostRaw{Float64}
    df = _read_csv(path)
    n = nrow(df)

    model = Vector{Int}(undef, n)
    startup = Vector{Float64}(undef, n)
    shutdown = Vector{Float64}(undef, n)
    ncost = Vector{Int}(undef, n)
    coeffs = Vector{Vector{Float64}}(undef, n)

    smodel = _sym("model")
    sstartup = _sym("startup")
    sshutdown = _sym("shutdown")
    sncost = _sym("ncost")

    @inbounds for i in 1:n
        model[i] = _safe_int(df[i, smodel], 0)
        startup[i] = _safe_f64(df[i, sstartup], 0.0)
        shutdown[i] = _safe_f64(df[i, sshutdown], 0.0)
        ncost[i] = _safe_int(df[i, sncost], 0)

        k = ncost[i]
        if k <= 0
            coeffs[i] = Float64[]
            continue
        end

        c = Vector{Float64}(undef, k)
        @inbounds for j in 1:k
            sj = _sym("c$(j)")
            if sj in propertynames(df)
                c[j] = _safe_f64(df[i, sj], 0.0)
            else
                c[j] = 0.0
            end
        end
        coeffs[i] = c
    end

    return GencostRaw{Float64}(model, startup, shutdown, ncost, coeffs)
end

function _parse_gencost_timeseries_csv(path::String)::GencostTimeseriesRaw{Float64}
    df = _read_csv(path)
    n = nrow(df)

    time = Vector{Int}(undef, n)
    gen_id = Vector{Int}(undef, n)
    model = Vector{Int}(undef, n)
    startup = Vector{Float64}(undef, n)
    shutdown = Vector{Float64}(undef, n)
    ncost = Vector{Int}(undef, n)
    coeffs = Vector{Vector{Float64}}(undef, n)
    tou = Vector{Float64}(undef, n)
    fit = Vector{Float64}(undef, n)

    stime = _sym("time")
    sgen_id = _sym("gen_id")
    smodel = _sym("model")
    sstartup = _sym("startup")
    sshutdown = _sym("shutdown")
    sncost = _sym("ncost")

    stou = _sym("tou")
    sfit = _sym("feed_in_tariff")
    has_tou = stou in propertynames(df)
    has_fit = sfit in propertynames(df)

    @inbounds for i in 1:n
        time[i] = _safe_int(df[i, stime], 1)
        gen_id[i] = _safe_int(df[i, sgen_id], 0)
        model[i] = _safe_int(df[i, smodel], 0)
        startup[i] = _safe_f64(df[i, sstartup], 0.0)
        shutdown[i] = _safe_f64(df[i, sshutdown], 0.0)
        ncost[i] = _safe_int(df[i, sncost], 0)

        k = ncost[i]
        if k <= 0
            coeffs[i] = Float64[]
        else
            c = Vector{Float64}(undef, k)
            @inbounds for j in 1:k
                sj = _sym("c$(j)")
                if sj in propertynames(df)
                    c[j] = _safe_f64(df[i, sj], 0.0)
                else
                    c[j] = 0.0
                end
            end
            coeffs[i] = c
        end

        tou[i] = has_tou ? _safe_f64(df[i, stou], 0.0) : 0.0
        fit[i] = has_fit ? _safe_f64(df[i, sfit], 0.0) : 0.0
    end

    return GencostTimeseriesRaw{Float64}(
        time, gen_id, model, startup, shutdown, ncost, coeffs, tou, fit,
    )
end

# -------------------------
# Prosumer timeseries: compat (p_pv/p_load/p_curtail) ou (generation/demand)
# -------------------------

function _prosumer_ts_from_tbl_csv(tbl::Dict{String, Any})
    n = length(tbl["time"])

    if haskey(tbl, "p_pv") && haskey(tbl, "p_load") && haskey(tbl, "p_curtail")
        return ProsumerTimeseriesRaw{Float64}(
            tbl["time"]::Vector{Int},
            tbl["bus"]::Vector{Int},
            tbl["p_pv"]::Vector{Float64},
            tbl["p_load"]::Vector{Float64},
            tbl["p_curtail"]::Vector{Float64},
        )
    end

    if haskey(tbl, "generation") && haskey(tbl, "demand")
        return ProsumerTimeseriesRaw{Float64}(
            tbl["time"]::Vector{Int},
            tbl["bus"]::Vector{Int},
            tbl["generation"]::Vector{Float64},
            tbl["demand"]::Vector{Float64},
            fill(0.0, n),
        )
    end

    @warn "prosumer_timeseries: missing expected columns. Filling zeros."
    return ProsumerTimeseriesRaw{Float64}(
        tbl["time"]::Vector{Int},
        tbl["bus"]::Vector{Int},
        fill(0.0, n),
        fill(0.0, n),
        fill(0.0, n),
    )
end

function _parse_case(dir::String, ::CsvFolder)::CaseRaw
    d = abspath(String(dir))

    case = CaseRaw(;
        name = splitdir(d)[end],
        source_type = "csv",
        source_version = "2",
        baseMVA = 100.0,
        T = 1,
    )

    base_path = joinpath(d, "baseMVA.csv")
    if isfile(base_path)
        df = _read_csv(base_path)
        if nrow(df) > 0 && _hascol(df, "baseMVA")
            x = df[1, _sym("baseMVA")]
            case.baseMVA = (x === missing) ? 100.0 : Float64(x)
        else
            @warn "baseMVA.csv empty or missing baseMVA. Using 100.0"
        end
    else
        @warn "baseMVA.csv not found. Using 100.0"
    end

    # -------------------------
    # core tables
    # -------------------------

    let p = joinpath(d, "bus.csv")
        if isfile(p)
            tbl = _parse_csv_table_dict(p, _mp_bus_columns_csv)
            case.bus = _bus_from_tbl(tbl)
        end
    end

    let p = joinpath(d, "gen.csv")
        if isfile(p)
            tbl = _parse_csv_table_dict(p, _mp_gen_columns_csv)
            n = length(tbl["gen_bus"])

            if !haskey(tbl, "plant_type")
                tbl["plant_type"] = fill("T", n)
            else
                pt = tbl["plant_type"]
                if !(pt isa Vector{String})
                    tbl["plant_type"] = String.(pt)
                end
                @inbounds for i in 1:n
                    if isempty(tbl["plant_type"][i])
                        tbl["plant_type"][i] = "T"
                    end
                end
            end

            if !haskey(tbl, "curtailment_cost")
                tbl["curtailment_cost"] = fill(0.0, n)
            end

            case.gen = _gen_from_tbl(tbl)
        end
    end

    let p = joinpath(d, "branch.csv")
        if isfile(p)
            tbl = _parse_csv_table_dict(p, _mp_branch_columns_csv)
            # NÃO converte angmin/angmax aqui.
            case.branch = _branch_from_tbl(tbl)
        end
    end

    let p = joinpath(d, "dcline.csv")
        if isfile(p)
            tbl = _parse_csv_table_dict(p, _mp_dcline_columns_csv)
            case.dcline = _dcline_from_tbl(tbl)
        end
    end

    let p = joinpath(d, "storage.csv")
        if isfile(p)
            tbl = _parse_csv_table_dict(p, _mp_storage_columns_csv)
            case.storage = _storage_from_tbl(tbl)
        end
    end

    let p = joinpath(d, "switch.csv")
        if isfile(p)
            tbl = _parse_csv_table_dict(p, _mp_switch_columns_csv)
            case.switch = _switch_from_tbl(tbl)
        end
    end

    let p = joinpath(d, "prosumer.csv")
        if isfile(p)
            tbl = _parse_csv_table_dict(p, _mp_prosumer_columns_csv)
            case.prosumer = _prosumer_from_tbl(tbl)
        end
    end

    let p = joinpath(d, "consumer.csv")
        if isfile(p)
            tbl = _parse_csv_table_dict(p, _mp_consumer_columns_csv)
            case.consumer = _consumer_from_tbl(tbl)
        end
    end

    let p = joinpath(d, "bus_name.csv")
        if isfile(p)
            tbl = _parse_csv_table_dict(p, _mp_bus_name_columns_csv)
            case.bus_name = _busname_from_tbl(tbl)
            if case.bus !== nothing && length(case.bus_name.name) != length(case.bus.bus_i)
                @error "Incorrect CSV case: number of bus names != number of buses."
            end
        end
    end

    let p = joinpath(d, "gencost.csv")
        if isfile(p)
            case.gencost = _parse_gencost_csv(p)
        end
    end

    # -------------------------
    # timeseries folder
    # -------------------------
    tdir = joinpath(d, "time")
    if isdir(tdir)
        let p = joinpath(tdir, "bus_timeseries.csv")
            if isfile(p)
                tbl = _parse_csv_table_dict(p, _mp_bus_timeseries_columns_csv)
                case.bus_timeseries = _bus_ts_from_tbl(tbl)
                if !isempty(case.bus_timeseries.time)
                    case.T = maximum(case.bus_timeseries.time)
                end
            end
        end

        let p = joinpath(tdir, "consumer_timeseries.csv")
            if isfile(p)
                tbl = _parse_csv_table_dict(p, _mp_consumer_timeseries_columns_csv)
                case.consumer_timeseries = _consumer_ts_from_tbl(tbl)
                if !isempty(case.consumer_timeseries.time)
                    case.T = max(case.T, maximum(case.consumer_timeseries.time))
                end
            end
        end

        let p = joinpath(tdir, "prosumer_timeseries.csv")
            if isfile(p)
                df = _read_csv(p)
                tbl = Dict{String, Any}()

                tbl["time"] = _col_or_default(df, "time", Int, 1)
                tbl["bus"] = _col_or_default(df, "bus", Int, 1)

                if _hascol(df, "p_pv") || _hascol(df, "p_load") || _hascol(df, "p_curtail")
                    tbl["p_pv"] = _col_or_default(df, "p_pv", Float64, 0.0)
                    tbl["p_load"] = _col_or_default(df, "p_load", Float64, 0.0)
                    tbl["p_curtail"] = _col_or_default(df, "p_curtail", Float64, 0.0)
                else
                    tbl["generation"] = _col_or_default(df, "generation", Float64, 0.0)
                    tbl["demand"] = _col_or_default(df, "demand", Float64, 0.0)
                end

                case.prosumer_timeseries = _prosumer_ts_from_tbl_csv(tbl)
                if !isempty(case.prosumer_timeseries.time)
                    case.T = max(case.T, maximum(case.prosumer_timeseries.time))
                end
            end
        end

        let p = joinpath(tdir, "gencost_timeseries.csv")
            if isfile(p)
                case.gencost_timeseries = _parse_gencost_timeseries_csv(p)
                if !isempty(case.gencost_timeseries.time)
                    case.T = max(case.T, maximum(case.gencost_timeseries.time))
                end
            end
        end
    end

    return case
end
