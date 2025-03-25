# ---------------------------------------------------------------------
# Specific data of Matpower
# ---------------------------------------------------------------------
const _mp_data_names = ["mpc.version", "mpc.baseMVA", "mpc.bus", "mpc.gen",
    "mpc.branch", "mpc.dcline", "mpc.gencost", "mpc.dclinecost",
    "mpc.bus_name", "mpc.storage", "mpc.switch",
]

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

const _mp_bus_name_columns = [
    ("name", Union{String, SubString{String}})
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
    ("plant_type", Float64),
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

const _mp_dcline_columns = [
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

const _mp_storage_columns = [
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

const _mp_switch_columns = [
    ("f_bus", Int), ("t_bus", Int),
    ("psw", Float64), ("qsw", Float64), ("state", Int),
    ("thermal_rating", Float64),
    ("status", Int),
]

const _mp_gencost_columns = [
    ("model", Int),
    ("startup", Float64),
    ("shutdown", Float64),
    ("ncost", Int),
    ("cost_coefficients", Vector{Float64}),
]

# ---------------------------------------------------------------------
# Matpower Parser
# ---------------------------------------------------------------------
function _parse_table(data::Array, columns::Vector{Tuple{String, DataType}})
    table = Dict{String, Vector}()
    for (column_name, _) in columns
        table[column_name] = Vector{Any}()
    end

    for row in data
        for (col_index, (column_name, _)) in enumerate(columns)
            if col_index <= length(row)
                push!(table[column_name], row[col_index])
            else
                @warn "Column '$column_name' does not exist in the data and will be ignored."
            end
        end
    end
    return table
end

function _parse_gencost_table(data::Array)
    table = Dict{String, Vector}()

    table["model"] = Vector{Int}()
    table["startup"] = Vector{Float64}()
    table["shutdown"] = Vector{Float64}()
    table["ncost"] = Vector{Int}()
    table["cost_coefficients"] = Vector{Vector{Float64}}()

    for row in data
        push!(table["model"], row[1])
        push!(table["startup"], row[2])
        push!(table["shutdown"], row[3])
        push!(table["ncost"], row[4])

        ncost = row[4]
        cost_coefficients = row[5:5+ncost-1]
        push!(table["cost_coefficients"], cost_coefficients)
    end

    return table
end

function _parse_matpower_string(data_string::String)
    matlab_data, func_name, colnames =
        InfrastructureModels.parse_matlab_string(data_string; extended = true)

    case = Dict{String, Any}()

    if func_name != nothing
        case["name"] = func_name
    else
        @warn(
            "No case name found in matpower file. The file seems to be missing \"function mpc = ...\""
        )
        case["name"] = "no_name_found"
    end

    case["source_type"] = "matpower"
    if haskey(matlab_data, "mpc.version")
        case["source_version"] = matlab_data["mpc.version"]
    else
        @warn(
            "no case version found in matpower file. The file seems to be missing \"mpc.version = ...\""
        )
        case["source_version"] = "0.0.0+"
    end

    if haskey(matlab_data, "mpc.baseMVA")
        case["baseMVA"] = matlab_data["mpc.baseMVA"]
    else
        @warn(
            "no baseMVA found in matpower file. The file seems to be missing \"mpc.baseMVA = ...\""
        )
        case["baseMVA"] = 1.0
    end

    if haskey(matlab_data, "mpc.bus")
        case["bus"] = _parse_table(matlab_data["mpc.bus"], _mp_bus_columns)
    else
        @error(
            "No bus table found in matpower file. The file seems to be missing \"mpc.bus = [...];\""
        )
    end

    if haskey(matlab_data, "mpc.gen")
        case["gen"] = _parse_table(matlab_data["mpc.gen"], _mp_gen_columns)

        ngens = length(case["gen"]["gen_bus"])
        if !haskey(case["gen"], "plant_type") || isempty(case["gen"]["plant_type"])
            @warn "Column 'plant_type' not found in mpc.gen. Assuming all plants are thermal ('T')."
            case["gen"]["plant_type"] = fill("T", ngens)
        end
        if !haskey(case["gen"], "curtailment_cost") ||
           isempty(case["gen"]["curtailment_cost"])
            @warn "Column 'curtailment_cost' not found in mpc.gen. Assuming zero curtailment cost."
            case["gen"]["curtailment_cost"] = fill(0.0, ngens)
        end

    else
        @error(
            "No gen table found in matpower file. The file seems to be missing \"mpc.gen = [...];\""
        )
    end

    if haskey(matlab_data, "mpc.gencost")
        case["gencost"] = _parse_gencost_table(matlab_data["mpc.gencost"])
    else
        @error(
            "No gencost table found in matpower file. The file seems to be missing \"mpc.gencost = [...];\""
        )
    end

    if haskey(matlab_data, "mpc.branch")
        case["branch"] = _parse_table(matlab_data["mpc.branch"], _mp_branch_columns)
    else
        @error(
            "No branch table found in matpower file. The file seems to be missing \"mpc.branch = [...];\""
        )
    end

    if haskey(matlab_data, "mpc.branch")
        branch_table = _parse_table(matlab_data["mpc.branch"], _mp_branch_columns)

        if haskey(branch_table, "angmin")
            branch_table["angmin"] = [x * π / 180 for x in branch_table["angmin"]]
        end
        if haskey(branch_table, "angmax")
            branch_table["angmax"] = [x * π / 180 for x in branch_table["angmax"]]
        end

        case["branch"] = branch_table
    else
        @error "No branch table found in matpower file."
    end

    if haskey(matlab_data, "mpc.dcline")
        case["dcline"] = _parse_table(matlab_data["mpc.dcline"], _mp_dcline_columns)
    end

    if haskey(matlab_data, "mpc.storage")
        case["storage"] = _parse_table(matlab_data["mpc.storage"], _mp_storage_columns)
    end

    if haskey(matlab_data, "mpc.switch")
        case["switch"] = _parse_table(matlab_data["mpc.switch"], _mp_switch_columns)
    end

    if haskey(matlab_data, "mpc.bus_name")
        case["bus_name"] = _parse_table(matlab_data["mpc.bus_name"], _mp_bus_name_columns)

        if length(case["bus_name"]["name"]) != length(case["bus"]["bus_i"])
            @error(
                "Incorrect Matpower file: the number of bus names does not match the number of buses."
            )
        end
    end

    return case
end
