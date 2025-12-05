struct BranchData{T <: Real}
    n::Int
    s_base::T

    from::Vector{Int}
    to::Vector{Int}
    status::Vector{Int}
    global_index::Vector{Int}
    active::Vector{Int}

    r::Vector{T}
    x::Vector{T}
    b_sh::Vector{T}

    rate_a::Vector{T}
    rate_b::Vector{T}
    rate_c::Vector{T}

    tap::Vector{T}
    shift::Vector{T}
    angmin::Vector{T}
    angmax::Vector{T}

    tap_position::Vector{T}
    tap_angle::Vector{T}
    inv_tap::Vector{T}
    g::Vector{T}
    b::Vector{T}
    tr::Vector{T}
    ti::Vector{T}

    initial_voltage_magnitude::Vector{T}
    initial_voltage_angle::Vector{T}

    bus_ids::Vector{Int}
    bus_idx_map::Dict{Int, Int}
    zone::Vector{Int}
end

function get_branch_parameters(matpower_data::Dict{String, Any})
    branch_data = matpower_data["branch"]

    f_bus = branch_data["f_bus"]
    t_bus = branch_data["t_bus"]

    r = Float64.(branch_data["br_r"])
    x = Float64.(branch_data["br_x"])
    b_sh = Float64.(branch_data["br_b"])

    tap_raw = Float64.(branch_data["tap"])
    shift_deg = Float64.(branch_data["shift"])

    rate_a = Float64.(branch_data["rate_a"])
    rate_b = Float64.(branch_data["rate_b"])
    rate_c = Float64.(branch_data["rate_c"])

    angmin = Float64.(branch_data["angmin"])
    angmax = Float64.(branch_data["angmax"])

    status = branch_data["br_status"]

    n = length(f_bus)

    local_indices = get(matpower_data, "branch_local_indices", nothing)
    global_index =
        local_indices === nothing ? collect(1:n) : Vector{Int}(local_indices)

    tap_position = ones(Float64, n)
    tap_angle = zeros(Float64, n)
    inv_tap = ones(Float64, n)
    tr = zeros(Float64, n)
    ti = zeros(Float64, n)

    @inbounds for i in 1:n
        t = tap_raw[i]
        adeg = shift_deg[i]

        if isnan(t) || t == 0.0
            tap_position[i] = 1.0
            tap_angle[i] = 0.0
            inv_tap[i] = 1.0
        else
            tap_position[i] = t
            tap_angle[i] = deg2rad(adeg)
            inv_tap[i] = 1.0 / t
        end
        tr[i] = tap_position[i] * cos(tap_angle[i])
        ti[i] = tap_position[i] * sin(tap_angle[i])
    end

    active = Int[]
    @inbounds for i in eachindex(f_bus)
        if status[i] == 1
            push!(active, i)
        end
    end

    bus_ids = matpower_data["bus"]["bus_i"]
    bus_idx_map = get_bus_index_map(bus_ids)

    den = r .^ 2 .+ x .^ 2
    g = r ./ den
    b = -x ./ den

    return (
        n = n,
        from = f_bus,
        to = t_bus,
        status = status,
        global_index = global_index,
        active = active, r = r,
        x = x,
        b_sh = b_sh,
        rate_a = rate_a,
        rate_b = rate_b,
        rate_c = rate_c, tap = tap_raw,
        shift = shift_deg,
        angmin = angmin,
        angmax = angmax, tap_position = tap_position,
        tap_angle = tap_angle,
        inv_tap = inv_tap,
        g = g,
        b = b,
        tr = tr,
        ti = ti, s_base = Float64(matpower_data["baseMVA"]),
        initial_voltage_magnitude = Float64.(matpower_data["bus"]["vm"]),
        initial_voltage_angle = Float64.(matpower_data["bus"]["va"]), bus_ids = bus_ids,
        bus_idx_map = bus_idx_map,
        zone = matpower_data["bus"]["zone"],
    )
end

function build_branch_data(data::Dict{String, Any})::BranchData{Float64}
    p = get_branch_parameters(data)

    return BranchData{Float64}(
        p.n,
        p.s_base,
        p.from,
        p.to,
        p.status,
        p.global_index,
        p.active,
        p.r,
        p.x,
        p.b_sh,
        p.rate_a,
        p.rate_b,
        p.rate_c,
        p.tap,
        p.shift,
        p.angmin,
        p.angmax,
        p.tap_position,
        p.tap_angle,
        p.inv_tap,
        p.g,
        p.b,
        p.tr,
        p.ti,
        p.initial_voltage_magnitude,
        p.initial_voltage_angle,
        p.bus_ids,
        p.bus_idx_map,
        p.zone,
    )
end
