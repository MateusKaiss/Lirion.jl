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
