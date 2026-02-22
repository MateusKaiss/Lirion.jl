# ============================================================
# Branch data
# ============================================================

struct BranchData{T <: Real}
    n::Int
    s_base::T

    from::Vector{Int}
    to::Vector{Int}
    status::Vector{Int}
    global_index::Vector{Int}   # global branch ids (full case) for this local vector
    active::Vector{Int}         # local indices i where status[i] == 1

    r::Vector{T}
    x::Vector{T}
    b_sh::Vector{T}

    rate_a::Vector{T}
    rate_b::Vector{T}
    rate_c::Vector{T}

    tap::Vector{T}
    shift::Vector{T}            # radians (keep consistent with raw)
    angmin::Vector{T}
    angmax::Vector{T}

    tap_position::Vector{T}
    tap_angle::Vector{T}        # radians
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

# ------------------------------------------------------------------------------
# Tap + phase-shift helpers
# ------------------------------------------------------------------------------

@inline function _tap_params(tap_raw::T, shift_rad::T) where {T <: Real}
    if isnan(tap_raw) || tap_raw == zero(T)
        tap_pos = one(T)
        tap_ang = zero(T)
        inv_tap = one(T)
    else
        tap_pos = tap_raw
        tap_ang = shift_rad
        inv_tap = one(T) / tap_raw
    end

    tr = tap_pos * cos(tap_ang)
    ti = tap_pos * sin(tap_ang)
    return tap_pos, tap_ang, inv_tap, tr, ti
end

# ------------------------------------------------------------------------------
# Parameter extraction
# ------------------------------------------------------------------------------

function get_branch_parameters(raw::CaseRaw)
    @assert raw.branch !== nothing && raw.bus !== nothing "CaseRaw.branch/bus is nothing"
    br = raw.branch::BranchRaw{Float64}
    bus = raw.bus::BusRaw{Float64}

    from = br.f_bus
    to = br.t_bus

    r = br.br_r
    x = br.br_x
    b_sh = br.br_b

    tap_raw = br.tap
    shift_rad = br.shift

    rate_a = br.rate_a
    rate_b = br.rate_b
    rate_c = br.rate_c

    angmin = br.angmin
    angmax = br.angmax

    status = br.br_status
    n = length(from)

    local_indices = raw.dist.zone_branch_idx_in_full
    global_index = isempty(local_indices) ? collect(1:n) : Vector{Int}(local_indices)

    tap_position = ones(Float64, n)
    tap_angle = zeros(Float64, n)
    inv_tap = ones(Float64, n)
    tr = zeros(Float64, n)
    ti = zeros(Float64, n)

    @inbounds for i in 1:n
        tap_position[i], tap_angle[i], inv_tap[i], tr[i], ti[i] =
            _tap_params(tap_raw[i], shift_rad[i])
    end

    active = Int[]
    sizehint!(active, n)
    @inbounds for i in 1:n
        status[i] == 1 && push!(active, i)
    end

    bus_ids = bus.bus_i
    bus_idx_map = raw.maps.bus_idx_map

    den = r .^ 2 .+ x .^ 2
    g = r ./ den
    b = -x ./ den

    return (
        n = n,
        from = from,
        to = to,
        status = status,
        global_index = global_index,
        active = active, r = r,
        x = x,
        b_sh = b_sh, rate_a = rate_a,
        rate_b = rate_b,
        rate_c = rate_c, tap = tap_raw,
        shift = shift_rad,
        angmin = angmin,
        angmax = angmax, tap_position = tap_position,
        tap_angle = tap_angle,
        inv_tap = inv_tap,
        g = g,
        b = b,
        tr = tr,
        ti = ti, s_base = Float64(raw.baseMVA),
        initial_voltage_magnitude = Float64.(bus.vm),
        initial_voltage_angle = Float64.(bus.va), bus_ids = bus_ids,
        bus_idx_map = bus_idx_map,
        zone = bus.zone,
    )
end

# ------------------------------------------------------------------------------
# Builder
# ------------------------------------------------------------------------------

function build_branch_data(raw::CaseRaw)::BranchData{Float64}
    p = get_branch_parameters(raw)

    return BranchData{Float64}(
        p.n,
        p.s_base, p.from,
        p.to,
        p.status,
        p.global_index,
        p.active, p.r,
        p.x,
        p.b_sh, p.rate_a,
        p.rate_b,
        p.rate_c, p.tap,
        p.shift,
        p.angmin,
        p.angmax, p.tap_position,
        p.tap_angle,
        p.inv_tap,
        p.g,
        p.b,
        p.tr,
        p.ti, p.initial_voltage_magnitude,
        p.initial_voltage_angle, p.bus_ids,
        p.bus_idx_map,
        p.zone,
    )
end
