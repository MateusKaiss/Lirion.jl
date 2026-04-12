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

    # Cache columns once (avoids O(n) property lookups inside the loop).
    col_model = _getcol(df, "model")
    col_startup = _getcol(df, "startup")
    col_shutdown = _getcol(df, "shutdown")
    col_ncost = _getcol(df, "ncost")

    # Cache coefficient columns c1..cK once, growing as needed.
    max_k = 0
    @inbounds for i in 1:n
        k = _safe_int(col_ncost === nothing ? missing : col_ncost[i], 0)
        k > max_k && (max_k = k)
    end
    coeff_cols = Vector{Any}(undef, max_k)
    @inbounds for j in 1:max_k
        coeff_cols[j] = _getcol(df, "c$(j)")
    end

    @inbounds for i in 1:n
        model[i] = _safe_int(col_model === nothing ? missing : col_model[i], 0)
        startup[i] = _safe_f64(col_startup === nothing ? missing : col_startup[i], 0.0)
        shutdown[i] = _safe_f64(col_shutdown === nothing ? missing : col_shutdown[i], 0.0)
        ncost[i] = _safe_int(col_ncost === nothing ? missing : col_ncost[i], 0)

        k = ncost[i]
        if k <= 0
            coeffs[i] = Float64[]
            continue
        end

        c = Vector{Float64}(undef, k)
        @inbounds for j in 1:k
            col = coeff_cols[j]
            c[j] = col === nothing ? 0.0 : _safe_f64(col[i], 0.0)
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

    # Cache columns once.
    col_time = _getcol(df, "time")
    col_gen_id = _getcol(df, "gen_id")
    col_model = _getcol(df, "model")
    col_startup = _getcol(df, "startup")
    col_shutdown = _getcol(df, "shutdown")
    col_ncost = _getcol(df, "ncost")
    col_tou = _getcol(df, "tou")
    col_fit = _getcol(df, "feed_in_tariff")

    max_k = 0
    @inbounds for i in 1:n
        k = _safe_int(col_ncost === nothing ? missing : col_ncost[i], 0)
        k > max_k && (max_k = k)
    end
    coeff_cols = Vector{Any}(undef, max_k)
    @inbounds for j in 1:max_k
        coeff_cols[j] = _getcol(df, "c$(j)")
    end

    @inbounds for i in 1:n
        time[i] = _safe_int(col_time === nothing ? missing : col_time[i], 1)
        gen_id[i] = _safe_int(col_gen_id === nothing ? missing : col_gen_id[i], 0)
        model[i] = _safe_int(col_model === nothing ? missing : col_model[i], 0)
        startup[i] = _safe_f64(col_startup === nothing ? missing : col_startup[i], 0.0)
        shutdown[i] = _safe_f64(col_shutdown === nothing ? missing : col_shutdown[i], 0.0)
        ncost[i] = _safe_int(col_ncost === nothing ? missing : col_ncost[i], 0)

        k = ncost[i]
        if k <= 0
            coeffs[i] = Float64[]
        else
            c = Vector{Float64}(undef, k)
            @inbounds for j in 1:k
                col = coeff_cols[j]
                c[j] = col === nothing ? 0.0 : _safe_f64(col[i], 0.0)
            end
            coeffs[i] = c
        end

        tou[i] = col_tou === nothing ? 0.0 : _safe_f64(col_tou[i], 0.0)
        fit[i] = col_fit === nothing ? 0.0 : _safe_f64(col_fit[i], 0.0)
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
