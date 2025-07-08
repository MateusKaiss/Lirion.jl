# ===============================================================
# Parameter Extractors for PowerOptimization.jl
# All functions to preprocess and extract model-ready parameters
# ===============================================================

# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

"Gets a field from a dictionary with a fallback default"
getfieldsafe(data::Dict, field::String, default = missing) = get(data, field, default)

"Creates a bus index map from bus IDs"
function get_bus_index_map(bus_ids::AbstractVector)
    return Dict(b => i for (i, b) in enumerate(bus_ids))
end

# ---------------------------------------------------------------------
# Getting Bus Parameters
# ---------------------------------------------------------------------
function get_bus_parameters(bus_data)
    bus_ids = bus_data["bus_i"]
    reference_bus_idx = findfirst(==(3), bus_data["bus_type"])
    reference_bus = bus_ids[reference_bus_idx]
    vmin = Dict{Int, Float64}()
    vmax = Dict{Int, Float64}()
    vstart = Dict{Int, Union{Nothing, Float64}}()
    @inbounds for i in eachindex(bus_ids)
        b = bus_ids[i]
        vmin[b] = bus_data["vmin"][i]
        vmax[b] = bus_data["vmax"][i]
        vstart[b] = isnan(bus_data["vm"][i]) ? nothing : bus_data["vm"][i]
    end

    return (
        bus_ids = bus_ids,
        reference_bus = reference_bus,
        vmin = vmin,
        vmax = vmax,
        vstart = vstart,
    )
end

# ------------------------------------------------------------
# Getting Branch Parameters
# ------------------------------------------------------------
function get_branch_parameters(branch_data, matpower_data)
    resistance = branch_data["br_r"]
    reactance = branch_data["br_x"]
    tap_position_raw = branch_data["tap"]
    tap_angle_raw = branch_data["shift"]

    n = length(tap_position_raw)
    tap_position = ones(Float64, n)
    tap_angle = zeros(Float64, n)
    inv_tap = ones(Float64, n)

    @inbounds for i in 1:n
        t = tap_position_raw[i]
        a = tap_angle_raw[i]

        if isnan(t) || t == 0.0
            tap_position[i] = 1.0
            tap_angle[i] = 0.0
            inv_tap[i] = 1.0
        else
            tap_position[i] = t
            tap_angle[i] = deg2rad(a)
            inv_tap[i] = 1.0 / t
        end
    end

    active = Int[]
    @inbounds for i in eachindex(branch_data["f_bus"])
        if branch_data["br_status"][i] == 1
            push!(active, i)
        end
    end

    bus_ids = matpower_data["bus"]["bus_i"]

    den = resistance .^ 2 .+ reactance .^ 2
    g = resistance ./ den
    b = -reactance ./ den

    return (
        number_of_branches = length(branch_data["f_bus"]),
        bus_from = branch_data["f_bus"],
        bus_to = branch_data["t_bus"],
        bus_ids = bus_ids,
        bus_idx_map = get_bus_index_map(bus_ids),
        resistance = resistance,
        reactance = reactance,
        b_shunt = branch_data["br_b"],
        tap_position = tap_position,
        tap_angle = tap_angle,
        inv_tap = inv_tap,
        capacity = branch_data["rate_a"],
        min_ang = branch_data["angmin"],
        max_ang = branch_data["angmax"],
        existing = branch_data["br_status"],
        s_base = matpower_data["baseMVA"],
        initial_voltage_magnitude = matpower_data["bus"]["vm"],
        initial_voltage_angle = matpower_data["bus"]["va"],
        g = g,
        b = b,
        active = active,
    )
end

# ---------------------------------------------------------------------------
# Generator Parameters
# ---------------------------------------------------------------------------
function get_generator_parameters(gen_data)
    number_of_generators = length(gen_data["gen_bus"])
    start_pg = Dict{Int, Union{Nothing, Float64}}()
    start_qg = Dict{Int, Union{Nothing, Float64}}()
    for g in 1:number_of_generators
        start_pg[g] = isnan(gen_data["pg"][g]) ? nothing : gen_data["pg"][g]
        start_qg[g] = isnan(gen_data["qg"][g]) ? nothing : gen_data["qg"][g]
    end

    return (
        number_of_generators = number_of_generators,
        existing = gen_data["gen_status"] .== 1,
        pmn = gen_data["pmin"],
        pmx = gen_data["pmax"],
        qmn = gen_data["qmin"],
        qmx = gen_data["qmax"],
        initial_pg = start_pg,
        initial_qg = start_qg,
        plant_type = gen_data["plant_type"],
        curtailment_cost = gen_data["curtailment_cost"],
    )
end

# ---------------------------------------------------------------------
# Battery Parameters
# ---------------------------------------------------------------------
function get_battery_parameters(battery_data)
    return (
        number_of_batteries = length(battery_data["storage_bus"]),
        min_storage = battery_data["energy"],
        max_storage = battery_data["energy_rating"],
        charge_ramp = battery_data["charge_rating"],
        discharge_ramp = battery_data["discharge_rating"],
        initial_storage = battery_data["energy"],
        qmin = battery_data["qmin"],
        qmax = battery_data["qmax"],
        thermal_limit = battery_data["thermal_rating"],
        charge_efficiency = battery_data["charge_efficiency"],
        discharge_efficiency = battery_data["discharge_efficiency"],
        r = battery_data["r"],
        x = battery_data["x"],
        p_loss = battery_data["p_loss"],
        q_loss = battery_data["q_loss"],
        bus = battery_data["storage_bus"],
    )
end

# ---------------------------------------------------------------------
# Prosumer Parameters
# ---------------------------------------------------------------------
function get_prosumer_parameters(data)
    if haskey(data, "prosumer_timeseries")
        ts = data["prosumer_timeseries"]
        return (
            is_timeseries = true,
            number_of_prosumers = length(unique(ts["bus"])),
            time = ts["time"],
            bus_index = ts["bus"],
            generation = ts["generation"],
            demand = ts["demand"],
            q_generation = ts["q_generation"],
            q_demand = ts["q_demand"],
            cost_injection = haskey(data, "prosumer") ? data["prosumer"]["cost_injection"] : fill(0.0, length(unique(ts["bus"]))),
            inverter_capacity = haskey(data, "prosumer") ? data["prosumer"]["inverter_capacity"] : fill(Inf, length(unique(ts["bus"]))),
        )
    else
        s = data["prosumer"]
        if haskey(data, "bus_timeseries")
            times = unique(data["bus_timeseries"]["time"])
            ntimes = length(times)
            nb = length(s["bus"])
            return (
                is_timeseries = true,
                number_of_prosumers = nb,
                time = repeat(times, inner=nb),
                bus_index = repeat(s["bus"], outer=ntimes),
                generation = repeat(s["generation"], outer=ntimes),
                demand = repeat(s["demand"], outer=ntimes),
                q_generation = repeat(s["q_generation"], outer=ntimes),
                q_demand = repeat(s["q_demand"], outer=ntimes),
                cost_injection = repeat(s["cost_injection"], outer=ntimes),
                inverter_capacity = repeat(s["inverter_capacity"], outer=ntimes),
            )
        else
            return (
                is_timeseries = false,
                number_of_prosumers = length(s["bus"]),
                bus_index = s["bus"],
                generation = s["generation"],
                demand = s["demand"],
                q_generation = s["q_generation"],
                q_demand = s["q_demand"],
                cost_injection = s["cost_injection"],
                inverter_capacity = s["inverter_capacity"],
            )
        end
    end
end

# ---------------------------------------------------------------------
# Objective Parameters
# ---------------------------------------------------------------------
function get_objective_parameters(
    gen_data,
    battery_data,
    gencost_data,
    prosumer_data = Dict(),
    gencost_timeseries = nothing,
)
    has_timeseries = gencost_timeseries !== nothing

    if has_timeseries
        generator_costs = Dict{Int, Dict{Int, Vector{Float64}}}()
        gen_ids = gencost_timeseries["gen_id"]
        times = gencost_timeseries["time"]
        cost_coeffs = gencost_timeseries["cost_coefficients"]

        for i in eachindex(gen_ids)
            t = times[i]
            g = gen_ids[i]
            c = cost_coeffs[i]

            if !haskey(generator_costs, t)
                generator_costs[t] = Dict{Int, Vector{Float64}}()
            end

            generator_costs[t][g] = c
        end
    else
        generator_costs = Dict{Int, Dict{Int, Vector{Float64}}}()
        generator_costs[1] = Dict{Int, Vector{Float64}}()
        for (g, coefs) in enumerate(getfieldsafe(gencost_data, "cost_coefficients", []))
            generator_costs[1][g] = coefs
        end
    end

    return (
        number_of_generators = getfieldsafe(gen_data, "gen_bus", []) |> length,
        number_of_batteries = getfieldsafe(battery_data, "storage_bus", []) |> length,
        generator_costs = generator_costs,
        num_polynomial_coefs = getfieldsafe(gencost_data, "ncost", []),
        curtailment_cost = getfieldsafe(gen_data, "curtailment_cost", []),
        plant_type = getfieldsafe(gen_data, "plant_type", []),
        generator_status = getfieldsafe(gen_data, "gen_status", []),
        number_of_prosumers = getfieldsafe(prosumer_data, "bus", []) |> length,
        prosumer_costs = getfieldsafe(prosumer_data, "cost_injection", []),
        has_timeseries = has_timeseries,
    )
end

# ---------------------------------------------------------------------
# Power Balance Parameters (multi-temporal compatible)
# ---------------------------------------------------------------------
function get_power_balance_parameters(matpower_data)
    bus_ids = matpower_data["bus"]["bus_i"]
    return (
        bus_ids = bus_ids,
        bus_idx_map = get_bus_index_map(bus_ids),
        generator_bus_index = matpower_data["gen"]["gen_bus"],
        battery_bus_index = getfieldsafe(
            getfieldsafe(matpower_data, "storage", Dict()),
            "storage_bus",
            Int[],
        ),
        prosumer_bus_index = getfieldsafe(
            getfieldsafe(matpower_data, "prosumer", Dict()),
            "bus",
            Int[],
        ),
        demand_bus_index = bus_ids,
        branch_bus_from = matpower_data["branch"]["f_bus"],
        branch_bus_to = matpower_data["branch"]["t_bus"],
        shunt_conductance = matpower_data["bus"]["gs"],
        existing = matpower_data["gen"]["gen_status"],
        branch_status = matpower_data["branch"]["br_status"],
        has_timeseries = haskey(matpower_data, "bus_timeseries"),
        bus_to_generators = getfieldsafe(
            matpower_data,
            "bus_to_gens",
            Dict{Int, Vector{Int}}(),
        ),
        bus_to_batteries = getfieldsafe(
            matpower_data,
            "bus_to_batteries",
            Dict{Int, Vector{Int}}(),
        ),
        bus_to_prosumers = getfieldsafe(
            matpower_data,
            "bus_to_prosumers",
            Dict{Int, Vector{Int}}(),
        ),
        bus_to_branches_from = getfieldsafe(
            matpower_data,
            "bus_to_branches_from",
            Dict{Int, Vector{Int}}(),
        ),
        bus_to_branches_to = getfieldsafe(
            matpower_data,
            "bus_to_branches_to",
            Dict{Int, Vector{Int}}(),
        ),
    )
end

# ---------------------------------------------------------------------
# Reactive Power Balance Parameters
# ---------------------------------------------------------------------
function get_reactive_power_parameters(matpower_data)
    return (
        bus_ids = matpower_data["bus"]["bus_i"],
        shunt_susceptance = matpower_data["bus"]["bs"],
        generator_status = matpower_data["gen"]["gen_status"],
        has_timeseries = haskey(matpower_data, "bus_timeseries"),
        bus_idx_map = matpower_data["bus_idx_map"],
        branch_status = matpower_data["branch"]["br_status"],
        bus_to_generators = getfieldsafe(
            matpower_data,
            "bus_to_gens",
            Dict{Int, Vector{Int}}(),
        ),
        bus_to_batteries = getfieldsafe(
            matpower_data,
            "bus_to_batteries",
            Dict{Int, Vector{Int}}(),
        ),
        bus_to_prosumers = getfieldsafe(
            matpower_data,
            "bus_to_prosumers",
            Dict{Int, Vector{Int}}(),
        ),
        bus_to_branches_from = getfieldsafe(
            matpower_data,
            "bus_to_branches_from",
            Dict{Int, Vector{Int}}(),
        ),
        bus_to_branches_to = getfieldsafe(
            matpower_data,
            "bus_to_branches_to",
            Dict{Int, Vector{Int}}(),
        ),
    )
end
