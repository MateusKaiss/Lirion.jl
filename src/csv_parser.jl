# ---------------------------------------------------------------------
# Expected columns from CSV files (same as Matpower)
# ---------------------------------------------------------------------
const _mp_bus_columns = [
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

const _mp_gen_columns = [
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
    ("plant_type", String),
    ("curtailment_cost", Float64),
]

const _mp_branch_columns = [
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

const _mp_gencost_columns = [
    ("model", Int),
    ("startup", Float64),
    ("shutdown", Float64),
    ("ncost", Int),
    ("cost_coefficients", Vector{Float64}),
]

const _prosumer_columns = [
    ("bus", Int),
    ("generation", Float64),
    ("demand", Float64),
    ("q_generation", Float64),
    ("q_demand", Float64),
    ("cost_injection", Float64),
    ("inverter_capacity", Float64),
]

# ---------------------------------------------------------------------
# Function to read any common CSV table (except gencost)
# ---------------------------------------------------------------------
function parse_csv_table(path::String, columns::Vector{Tuple{String, DataType}})
    df = CSV.read(path, DataFrame)
    table = Dict{String, Vector}()

    for (col_name, col_type) in columns
        if hasproperty(df, col_name)
            table[col_name] = convert(Vector{col_type}, df[!, col_name])
        else
            @warn "Missing column '$col_name' in file $path. Filling with default values."
            default = col_type === String ? "" : zero(col_type)
            table[col_name] = fill(default, nrow(df))
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
# Main function: assembles the complete case from the folder with CSVs
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
            case["baseMVA"] = base[!, :baseMVA][1]
        else
            @warn "baseMVA.csv is empty or missing column 'baseMVA'. Using 100.0"
            case["baseMVA"] = 100.0
        end
    else
        @warn "baseMVA.csv not found. Using default value 100.0"
        case["baseMVA"] = 100.0
    end

    if isfile(joinpath(dir, "bus.csv"))
        case["bus"] = parse_csv_table(joinpath(dir, "bus.csv"), _mp_bus_columns)
    end

    if isfile(joinpath(dir, "gen.csv"))
        case["gen"] = parse_csv_table(joinpath(dir, "gen.csv"), _mp_gen_columns)
    end

    if isfile(joinpath(dir, "branch.csv"))
        case["branch"] = parse_csv_table(joinpath(dir, "branch.csv"), _mp_branch_columns)
    end

    if isfile(joinpath(dir, "gencost.csv"))
        case["gencost"] = parse_gencost_table(joinpath(dir, "gencost.csv"))
    end

    if isfile(joinpath(dir, "prosumer.csv"))
        case["prosumer"] = parse_csv_table(joinpath(dir, "prosumer.csv"), _prosumer_columns)
    end

    profile_dir = joinpath(dir, "time")
    if isdir(profile_dir)
        for file in readdir(profile_dir)
            fullpath = joinpath(profile_dir, file)
            case[file] = CSV.read(fullpath, DataFrame)
        end
    end

    return case
end