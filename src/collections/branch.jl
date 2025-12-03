struct BranchData
    n::Int
    from::Vector{Int}
    to::Vector{Int}
    r::Vector{Float64}
    x::Vector{Float64}
    b_sh::Vector{Float64}
    rate_a::Vector{Float64}
    rate_b::Vector{Float64}
    rate_c::Vector{Float64}
    tap::Vector{Float64}
    shift::Vector{Float64}
    status::Vector{Int}
    angmin::Vector{Float64}
    angmax::Vector{Float64}


    global_index::Vector{Int}
    bus_ids::Vector{Int}
    bus_idx_map::Dict{Int,Int}
    tap_position::Vector{Float64}
    tap_angle::Vector{Float64}
    inv_tap::Vector{Float64}
    g::Vector{Float64}
    b::Vector{Float64}
    tr::Vector{Float64}
    ti::Vector{Float64} 
    active::Vector{Int}
    s_base::Float64
    initial_voltage_magnitude::Vector{Float64}
    initial_voltage_angle::Vector{Float64}
    zone::Vector{Int}
end

function get_branch_parameters(branch_data, matpower_data)
    resistance       = branch_data["br_r"]
    reactance        = branch_data["br_x"]
    tap_position_raw = branch_data["tap"]
    tap_angle_raw    = branch_data["shift"]

    n = length(tap_position_raw)
    local_indices = get(matpower_data, "branch_local_indices", nothing)
    global_index =
        local_indices === nothing ? collect(1:n) : Vector{Int}(local_indices)

    tap_position = ones(Float64, n)
    tap_angle    = zeros(Float64, n)
    inv_tap      = ones(Float64, n)
    tr           = zeros(Float64, n)
    ti           = zeros(Float64, n)

    @inbounds for i in 1:n
        t = tap_position_raw[i]
        a_deg = tap_angle_raw[i]

        if isnan(t) || t == 0.0
            tap_position[i] = 1.0
            tap_angle[i]    = 0.0
            inv_tap[i]      = 1.0
        else
            tap_position[i] = t
            tap_angle[i]    = deg2rad(a_deg)
            inv_tap[i]      = 1.0 / t
        end
        tr[i] = tap_position[i] * cos(tap_angle[i])
        ti[i] = tap_position[i] * sin(tap_angle[i])
    end

    active = Int[]
    @inbounds for i in eachindex(branch_data["f_bus"])
        if branch_data["br_status"][i] == 1
            push!(active, i)
        end
    end

    bus_ids = matpower_data["bus"]["bus_i"]

    den = resistance .^ 2 .+ reactance .^ 2
    g   = resistance ./ den
    b   = -reactance ./ den

    return (
        number_of_branches         = length(branch_data["f_bus"]),
        bus_from                   = branch_data["f_bus"],
        bus_to                     = branch_data["t_bus"],
        bus_ids                    = bus_ids,
        bus_idx_map                = get_bus_index_map(bus_ids),
        resistance                 = resistance,
        reactance                  = reactance,
        b_shunt                    = branch_data["br_b"],
        tap_position               = tap_position,
        tap_angle                  = tap_angle,
        inv_tap                    = inv_tap,
        tr                         = tr,
        ti                         = ti,
        capacity                   = branch_data["rate_a"],
        min_ang                    = branch_data["angmin"],
        max_ang                    = branch_data["angmax"],
        existing                   = branch_data["br_status"],
        s_base                     = matpower_data["baseMVA"],
        initial_voltage_magnitude  = matpower_data["bus"]["vm"],
        initial_voltage_angle      = matpower_data["bus"]["va"],
        g                          = g,
        b                          = b,
        active                     = active,
        zone                       = matpower_data["bus"]["zone"],
        global_index               = global_index,
    )
end

function build_branch_data(data::Dict{String,Any})::BranchData
    branch = data["branch"]
    n      = length(branch["f_bus"])

    params = get_branch_parameters(branch, data)

    return BranchData(
        n,
        branch["f_bus"],
        branch["t_bus"],
        branch["br_r"],
        branch["br_x"],
        branch["br_b"],
        branch["rate_a"],
        branch["rate_b"],
        branch["rate_c"],
        branch["tap"],
        branch["shift"],
        branch["br_status"],
        branch["angmin"],
        branch["angmax"],

        params.global_index,
        params.bus_ids,
        params.bus_idx_map,
        params.tap_position,
        params.tap_angle,
        params.inv_tap,
        params.g,
        params.b,
        params.tr,
        params.ti,
        params.active,
        params.s_base,
        params.initial_voltage_magnitude,
        params.initial_voltage_angle,
        params.zone,
    )
end