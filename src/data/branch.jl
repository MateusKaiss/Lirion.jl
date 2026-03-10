# ============================================================
# Branch data
# ============================================================

struct BranchData
    n::Int
    s_base::Float64

    from::Vector{Int}
    to::Vector{Int}
    status::Vector{Int}
    global_index::Vector{Int}
    active::Vector{Int}

    r::Vector{Float64}
    x::Vector{Float64}
    b_sh::Vector{Float64}

    rate_a::Vector{Float64}
    rate_b::Vector{Float64}
    rate_c::Vector{Float64}

    tap::Vector{Float64}
    shift::Vector{Float64}
    angmin::Vector{Float64}
    angmax::Vector{Float64}

    tap_position::Vector{Float64}
    tap_angle::Vector{Float64}
    inv_tap::Vector{Float64}
    g::Vector{Float64}
    b::Vector{Float64}
    tr::Vector{Float64}
    ti::Vector{Float64}

    initial_voltage_magnitude::Vector{Float64}
    initial_voltage_angle::Vector{Float64}

    bus_ids::Vector{Int}
    bus_idx_map::Dict{Int, Int}
    zone::Vector{Int}
end

# ---------------------------------------------------------------------
# Tap + phase-shift helpers
# ---------------------------------------------------------------------
@inline function _tap_params(tap_raw::Float64, shift_rad::Float64)
    if isnan(tap_raw) || tap_raw == 0.0
        tap_pos = 1.0
        tap_ang = 0.0
        inv_tap = 1.0
    else
        tap_pos = tap_raw
        tap_ang = shift_rad
        inv_tap = 1.0 / tap_raw
    end

    tr = tap_pos * cos(tap_ang)
    ti = tap_pos * sin(tap_ang)

    return tap_pos, tap_ang, inv_tap, tr, ti
end

# ---------------------------------------------------------------------
# Branch Parameters (CaseRaw)
# ---------------------------------------------------------------------
function get_branch_parameters(case_raw::CaseRaw)
    @assert case_raw.branch !== nothing && case_raw.bus !== nothing "CaseRaw.branch/bus is nothing"

    br = case_raw.branch::BranchRaw{Float64}
    bus = case_raw.bus::BusRaw{Float64}

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

    local_indices = case_raw.dist.zone_branch_idx_in_full
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

    den = r .^ 2 .+ x .^ 2
    g = r ./ den
    b = -x ./ den

    return (
        n = n,
        s_base = case_raw.baseMVA,
        from = from,
        to = to,
        status = status,
        global_index = global_index,
        active = active,
        r = r,
        x = x,
        b_sh = b_sh,
        rate_a = rate_a,
        rate_b = rate_b,
        rate_c = rate_c,
        tap = tap_raw,
        shift = shift_rad,
        angmin = angmin,
        angmax = angmax,
        tap_position = tap_position,
        tap_angle = tap_angle,
        inv_tap = inv_tap,
        g = g,
        b = b,
        tr = tr,
        ti = ti,
        initial_voltage_magnitude = bus.vm,
        initial_voltage_angle = bus.va,
        bus_ids = bus.bus_i,
        bus_idx_map = case_raw.maps.bus_idx_map,
        zone = bus.zone,
    )
end

function build_branch_data(case_raw::CaseRaw)::BranchData
    params = get_branch_parameters(case_raw)

    return BranchData(
        params.n,
        params.s_base,
        params.from,
        params.to,
        params.status,
        params.global_index,
        params.active,
        params.r,
        params.x,
        params.b_sh,
        params.rate_a,
        params.rate_b,
        params.rate_c,
        params.tap,
        params.shift,
        params.angmin,
        params.angmax,
        params.tap_position,
        params.tap_angle,
        params.inv_tap,
        params.g,
        params.b,
        params.tr,
        params.ti,
        params.initial_voltage_magnitude,
        params.initial_voltage_angle,
        params.bus_ids,
        params.bus_idx_map,
        params.zone,
    )
end
