# ---------- helpers ----------
_now_str() = Dates.format(Dates.now(), dateformat"yyyy-mm-dd HH:MM:SS")

# Tenta extrair n_buses do .m usando o loader do Lirion
function _n_buses_from_case(path::AbstractString)::Int
    try
        data = Lirion.load_input_data(path)  # usa teu loader
        # tenta vários formatos comuns
        if haskey(data, :bus)
            return size(data[:bus], 1)
        elseif haskey(data, "bus")
            return size(data["bus"], 1)
        elseif haskey(data, :n_buses)
            return Int(data[:n_buses])
        elseif haskey(data, "n_buses")
            return Int(data["n_buses"])
        else
            # fallback: 0 se não achar
            return 0
        end
    catch
        return 0
    end
end

function _list_cases(dir::AbstractString)
    isdir(dir) || error("Diretório não existe: $dir")
    return sort(glob("*.m", dir))
end

# Extrai (obj, iters, tempo_s) do retorno distribuído
function _triple_from_distributed(result_any)
    if result_any isa Tuple && length(result_any) >= 4
        obj    = Float64(result_any[1])
        iters  = Int(result_any[3])
        time_s = Float64(result_any[4])
        return obj, iters, time_s
    end
    error("Retorno distribuído inesperado; esperado tuple (obj, ..., iters, time_s).")
end

# Executa Lirion.main e captura erro/resultado
function _call_lirion_main(argsvec::Vector{String})
    try
        res = Lirion.main(argsvec)
        return (status = "ok", result = res, error = "")
    catch err
        return (status = "error", result = nothing, error = sprint(showerror, err))
    end
end

# ---------- só ADMM AC (bus) ----------
"""
run_bus_ac_only(dir; extra_args=String[], outpath="csv_bus_ac_angles.csv")

- Roda apenas o modelo distribuído ADMM **AC/bus** para todos os .m do diretório.
- Gera CSV com: [:timestamp, :case, :family, :variant, :n_buses, :objective, :iters, :time_s, :status, :error]
"""
function run_bus_ac_only(dir::AbstractString;
                         extra_args::Vector{String}=String[],
                         outpath::Union{Nothing,String}="csv_bus_ac_angles.csv")

    cases = _list_cases(dir)
    isempty(cases) && @warn "Nenhum .m encontrado em $dir"

    rows = Vector{NamedTuple}()
    for f in cases
        modelkey = "distributed-admm-ac"   # family=:ac + variant=:bus
        args = String[f, "--model", modelkey]
        append!(args, extra_args)

        nb = _n_buses_from_case(f)

        @info "[ADMM AC/BUS] $(basename(f)) | $modelkey | n_buses=$nb"
        call = _call_lirion_main(args)

        if call.status == "ok"
            try
                obj, iters, time_s = _triple_from_distributed(call.result)
                push!(rows, (
                    timestamp = _now_str(),
                    case = f, family = "ac", variant = "bus",
                    n_buses = nb,
                    objective = obj, iters = iters, time_s = time_s,
                    status = "ok", error = "",
                ))
            catch err
                push!(rows, (
                    timestamp = _now_str(),
                    case = f, family = "ac", variant = "bus",
                    n_buses = nb,
                    objective = NaN, iters = -1, time_s = NaN,
                    status = "error", error = sprint(showerror, err),
                ))
            end
        else
            push!(rows, (
                timestamp = _now_str(),
                case = f, family = "ac", variant = "bus",
                n_buses = nb,
                objective = NaN, iters = -1, time_s = NaN,
                status = "error", error = call.error,
            ))
        end
    end

    df = DataFrame(rows)
    outcsv = isnothing(outpath) ? joinpath(dir, "csv_bus_ac_angles.csv") : outpath
    CSV.write(outcsv, df)
    @info "CSV (AC/BUS) salvo em: $outcsv"
    show(df); println()
    return df
end

# AC / ADMM (variant=:bus) — só os rhos relevantes
extra_ac_bus = String[
    "--admm-max-iter","5000",
    "--admm-rho","1e5",          # fallback (usado se θ ou V não forem passados)
    "--admm-rho-theta","1e5",    # ρθ
    "--admm-rho-v","1e5",        # ρV (|V|)
    "--admm-eps-primal","1e-4",
    "--admm-eps-dual","1e-4",
    "--admm-adaptive","on",
    "--admm-mu-incr","3.0",
    "--admm-mu-decr","3.0",
    "--admm-tau-incr","1.1",
    "--admm-tau-decr","0.9",
]

# 57-bus
run_bus_ac_only("/home/mateus-kaiss/Projects/experiments/57-bus";
    extra_args=extra_ac_bus,
    outpath="csv_bus_ac_angles.csv")

# 200-bus
run_bus_ac_only("/home/mateus-kaiss/Projects/experiments/200-bus";
    extra_args=extra_ac_bus,
    outpath="csv_bus_ac_angles.csv")

# 300-bus
run_bus_ac_only("/home/mateus-kaiss/Projects/experiments/300-bus";
    extra_args=extra_ac_bus,
    outpath="csv_bus_ac_angles.csv")

# 1354-bus
run_bus_ac_only("/home/mateus-kaiss/Projects/experiments/1354-bus";
    extra_args=extra_ac_bus,
    outpath="csv_bus_ac_angles.csv")
