# ---------------------------------------------------------------------
# Expected columns from CSV files (same as Matpower)
# ---------------------------------------------------------------------
const _mp_data_names_csv = [
    "mpc.version", "mpc.baseMVA", "mpc.bus", "mpc.gen",
    "mpc.branch", "mpc.dcline", "mpc.gencost", "mpc.dclinecost",
    "mpc.bus_name", "mpc.storage", "mpc.switch", "mpc.prosumer",
    "mpc.bus_timeseries", "mpc.prosumer_timeseries",
    "mpc.consumer"
]

const _mp_bus_columns_csv = [
    ("bus_i", Int),
    ("bus_type", Int),
    ("pd", Float64), ("qd", Float64),
    ("gs", Float64), ("bs", Float64),
    ("area", Int),
    ("vm", Float64), ("va", Float64),
    ("base_kv", Float64),
    ("zone", Int),
    ("vmax", Float64), ("vmin", Float64),
    ("lam_p", Float64), ("lam_q", Float64),
    ("mu_vmax", Float64), ("mu_vmin", Float64),
]

const _mp_bus_name_columns_csv = [
    ("name", Union{String, SubString{String}})
]

const _mp_gen_columns_csv = [
    ("gen_bus", Int),
    ("pg", Float64), ("qg", Float64),
    ("qmax", Float64), ("qmin", Float64),
    ("vg", Float64),
    ("mbase", Float64),
    ("gen_status", Int),
    ("pmax", Float64), ("pmin", Float64),
    ("pc1", Float64),
    ("pc2", Float64),
    ("qc1min", Float64), ("qc1max", Float64),
    ("qc2min", Float64), ("qc2max", Float64),
    ("ramp_agc", Float64),
    ("ramp_10", Float64),
    ("ramp_30", Float64),
    ("ramp_q", Float64),
    ("apf", Float64),
    ("mu_pmax", Float64), ("mu_pmin", Float64),
    ("mu_qmax", Float64), ("mu_qmin", Float64),
    ("plant_type", Float64),
    ("curtailment_cost", Float64),
]

const _mp_branch_columns_csv = [
    ("f_bus", Int),
    ("t_bus", Int),
    ("br_r", Float64), ("br_x", Float64),
    ("br_b", Float64),
    ("rate_a", Float64),
    ("rate_b", Float64),
    ("rate_c", Float64),
    ("tap", Float64), ("shift", Float64),
    ("br_status", Int),
    ("angmin", Float64), ("angmax", Float64),
    ("pf", Float64), ("qf", Float64),
    ("pt", Float64), ("qt", Float64),
    ("mu_sf", Float64), ("mu_st", Float64),
    ("mu_angmin", Float64), ("mu_angmax", Float64),
]

const _mp_dcline_columns_csv = [
    ("f_bus", Int),
    ("t_bus", Int),
    ("br_status", Int),
    ("pf", Float64), ("pt", Float64),
    ("qf", Float64), ("qt", Float64),
    ("vf", Float64), ("vt", Float64),
    ("pmin", Float64), ("pmax", Float64),
    ("qminf", Float64), ("qmaxf", Float64),
    ("qmint", Float64), ("qmaxt", Float64),
    ("loss0", Float64),
    ("loss1", Float64),
    ("mu_pmin", Float64), ("mu_pmax", Float64),
    ("mu_qminf", Float64), ("mu_qmaxf", Float64),
    ("mu_qmint", Float64), ("mu_qmaxt", Float64),
]

const _mp_storage_columns_csv = [
    ("storage_bus", Int),
    ("ps", Float64), ("qs", Float64),
    ("energy", Float64), ("energy_rating", Float64),
    ("charge_rating", Float64), ("discharge_rating", Float64),
    ("charge_efficiency", Float64), ("discharge_efficiency", Float64),
    ("thermal_rating", Float64),
    ("qmin", Float64), ("qmax", Float64),
    ("r", Float64), ("x", Float64),
    ("p_loss", Float64), ("q_loss", Float64),
    ("status", Int),
]

const _mp_switch_columns_csv = [
    ("f_bus", Int), ("t_bus", Int),
    ("psw", Float64), ("qsw", Float64), ("state", Int),
    ("thermal_rating", Float64),
    ("status", Int),
]

const _mp_gencost_columns_csv = [
    ("model", Int),
    ("startup", Float64),
    ("shutdown", Float64),
    ("ncost", Int),
    ("cost_coefficients", Vector{Float64}),
]

const _mp_prosumer_columns_csv = [
    ("bus", Int),
    ("generation", Float64),
    ("demand", Float64),
    ("q_generation", Float64),
    ("q_demand", Float64),
    ("cost_injection", Float64),
    ("inverter_capacity", Float64),
]

const _mp_consumer_columns_csv = [
    ("bus", Int),
    ("demand", Float64),
    ("q_demand", Float64),
]

const _mp_consumer_timeseries_columns_csv = [
    ("time", Int),
    ("bus", Int),
    ("demand", Float64),
    ("q_demand", Float64),
]

const _mp_prosumer_timeseries_columns_csv = [
    ("time", Int),
    ("bus", Int),
    ("generation", Float64),
    ("demand", Float64),
    ("q_generation", Float64),
    ("q_demand", Float64),
]

const _mp_gencost_timeseries_columns_csv = [
    ("time", Int),
    ("gen_id", Int),
    ("model", Int),
    ("startup", Float64),
    ("shutdown", Float64),
    ("ncost", Int),
    ("cost_coefficients", Vector{Float64}),
    ("tou", Float64),
    ("feed_in_tariff", Float64),
]

const _mp_bus_timeseries_columns_csv = [
    ("time", Int),
    ("bus", Int),
    ("pd", Float64),
    ("qd", Float64),
]

function parse_csv_table(path::String, columns::Vector{Tuple{String, DataType}})
    df = CSV.read(path, DataFrame)
    table = Dict{String, Union{Vector, Nothing}}()

    for (col_name, col_type) in columns
        if hasproperty(df, col_name)
            table[col_name] = convert(Vector{col_type}, df[!, col_name])
        else
            @warn "Missing column '$col_name' in file $path. Setting as `nothing`."
            table[col_name] = Any[]
        end
    end

    return table
end

# ---------------------------------------------------------------------
# Specific function for gencost.csv
# ---------------------------------------------------------------------
function parse_gencost_table(path::String)
    df = CSV.read(path, DataFrame)
    table = Dict{String, Vector}()

    table["model"] = Vector{Int}()
    table["startup"] = Vector{Float64}()
    table["shutdown"] = Vector{Float64}()
    table["ncost"] = Vector{Int}()
    table["cost_coefficients"] = Vector{Vector{Float64}}()

    for row in eachrow(df)
        push!(table["model"], row.model)
        push!(table["startup"], row.startup)
        push!(table["shutdown"], row.shutdown)
        push!(table["ncost"], row.ncost)
        coeffs = []
        for i in 1:row.ncost
            push!(coeffs, row[Symbol("c$(i)")])
        end
        push!(table["cost_coefficients"], coeffs)
    end
    return table
end

# ---------------------------------------------------------------------
# Parse Gencost Timeseries CSV
# ---------------------------------------------------------------------
function parse_gencost_timeseries_csv(path::String)
    df = CSV.read(path, DataFrame)
    table = Dict{String, Vector}()

    table["time"] = Int[]
    table["gen_id"] = Int[]
    table["model"] = Int[]
    table["startup"] = Float64[]
    table["shutdown"] = Float64[]
    table["ncost"] = Int[]
    table["cost_coefficients"] = Vector{Vector{Float64}}()
    table["tou"] = Float64[]
    table["feed_in_tariff"] = Float64[]

    for row in eachrow(df)
        n = row.ncost
        coeffs = [row[Symbol("c$(i)")] for i in 1:n]
        push!(table["time"], row.time)
        push!(table["gen_id"], row.gen_id)
        push!(table["model"], row.model)
        push!(table["startup"], row.startup)
        push!(table["shutdown"], row.shutdown)
        push!(table["ncost"], n)
        push!(table["cost_coefficients"], coeffs)
        push!(table["tou"], row.tou)
        push!(table["feed_in_tariff"], row.feed_in_tariff)
    end
    return table
end

# ---------------------------------------------------------------------
# Main CSV Parser Function
# ---------------------------------------------------------------------
function parse_case_from_csv(dir::String)
    case = Dict{String, Any}()
    case["source_type"] = "csv"
    case["source_version"] = "2"
    case["name"] = splitdir(dir)[end]

    base_path = joinpath(dir, "baseMVA.csv")
    if isfile(base_path)
        base = CSV.read(base_path, DataFrame)
        if nrow(base) > 0 && hasproperty(base, :baseMVA)
            case["baseMVA"] = base[1, :baseMVA]
        else
            @warn "baseMVA.csv is empty or missing 'baseMVA'. Using 100.0"
            case["baseMVA"] = 100.0
        end
    else
        @warn "baseMVA.csv not found. Using default value 100.0"
        case["baseMVA"] = 100.0
    end

    for (name, cols) in [
        ("bus", _mp_bus_columns_csv),
        ("gen", _mp_gen_columns_csv),
        ("branch", _mp_branch_columns_csv),
        ("dcline", _mp_dcline_columns_csv),
        ("storage", _mp_storage_columns_csv),
        ("switch", _mp_switch_columns_csv),
        ("prosumer", _mp_prosumer_columns_csv),
        ("consumer", _mp_consumer_columns_csv),
        ("bus_name", _mp_bus_name_columns_csv),
    ]
        path = joinpath(dir, "$name.csv")
        if isfile(path)
            case[name] = parse_csv_table(path, cols)
        end
    end

    if haskey(case, "branch")
        if haskey(case["branch"], "angmin")
            case["branch"]["angmin"] = [deg2rad(x) for x in case["branch"]["angmin"]]
        end
        if haskey(case["branch"], "angmax")
            case["branch"]["angmax"] = [deg2rad(x) for x in case["branch"]["angmax"]]
        end
    end

    if isfile(joinpath(dir, "gencost.csv"))
        case["gencost"] = parse_gencost_table(joinpath(dir, "gencost.csv"))
    end

    if haskey(case, "gen")
        ngens = length(case["gen"]["gen_bus"])
        if !haskey(case["gen"], "plant_type") || isempty(case["gen"]["plant_type"])
            @warn "plant_type missing in gen.csv. Assuming 'T'."
            case["gen"]["plant_type"] = fill("T", ngens)
        end
        if !haskey(case["gen"], "curtailment_cost") ||
           isempty(case["gen"]["curtailment_cost"])
            @warn "curtailment_cost missing in gen.csv. Assuming 0.0."
            case["gen"]["curtailment_cost"] = fill(0.0, ngens)
        end
    end

    profile_dir = joinpath(dir, "time")
    if isdir(profile_dir)
        for file in readdir(profile_dir)
            fullpath = joinpath(profile_dir, file)
            key = splitext(file)[1]

            if key == "gencost_timeseries"
                case[key] = parse_gencost_timeseries_csv(fullpath)
            elseif key == "bus_timeseries"
                case[key] = parse_csv_table(fullpath, _mp_bus_timeseries_columns_csv)
            elseif key == "prosumer_timeseries"
                case[key] = parse_csv_table(fullpath, _mp_prosumer_timeseries_columns_csv)
            elseif key == "consumer_timeseries"
                case[key] = parse_csv_table(fullpath, _mp_consumer_timeseries_columns_csv)
            else
                @warn "Unknown timeseries file: $file. Reading as raw DataFrame."
                case[key] = CSV.read(fullpath, DataFrame)
            end
        end
    end

    return case
end
