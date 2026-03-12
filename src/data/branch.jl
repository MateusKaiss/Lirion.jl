# ============================================================
# Branch data
# ============================================================

struct BranchData
    n::Int
    s_base::Float64

    from::Vector{Int}
    to::Vector{Int}
    status::Vector{Int}
    active::Vector{Int}

    r::Vector{Float64}
    x::Vector{Float64}
    b_sh::Vector{Float64}
    rate_a::Vector{Float64}
    rate_b::Vector{Float64}
    rate_c::Vector{Float64}

    tap::Vector{Float64}
    shift::Vector{Float64}
    ang_min::Vector{Float64}
    ang_max::Vector{Float64}

    tap_position::Vector{Float64}
    tap_angle::Vector{Float64}
    inv_tap::Vector{Float64}

    g::Vector{Float64}
    b::Vector{Float64}

    tap_real::Vector{Float64}
    tap_imag::Vector{Float64}
end

@inline function _tap_params(tap_raw::Float64, shift_rad::Float64)
    if isnan(tap_raw) || tap_raw == 0.0
        tap_position = 1.0
        tap_angle = 0.0
        inv_tap = 1.0
    else
        tap_position = tap_raw
        tap_angle = shift_rad
        inv_tap = 1.0 / tap_raw
    end

    tap_real = tap_position * cos(tap_angle)
    tap_imag = tap_position * sin(tap_angle)

    return tap_position, tap_angle, inv_tap, tap_real, tap_imag
end

function get_branch_parameters(case_raw::CaseRaw)
    @assert case_raw.branch !== nothing && case_raw.bus !== nothing "CaseRaw.branch/bus is nothing"

    branch_raw = case_raw.branch::BranchRaw{Float64}

    from = branch_raw.f_bus
    to = branch_raw.t_bus
    status = branch_raw.br_status
    tap = branch_raw.tap
    shift = branch_raw.shift
    r = branch_raw.br_r
    x = branch_raw.br_x
    denominator = r .^ 2 .+ x .^ 2
    g = r ./ denominator
    b = -x ./ denominator

    n = length(from)

    tap_position = Vector{Float64}(undef, n)
    tap_angle = Vector{Float64}(undef, n)
    inv_tap = Vector{Float64}(undef, n)
    tap_real = Vector{Float64}(undef, n)
    tap_imag = Vector{Float64}(undef, n)

    @inbounds for i in 1:n
        tap_position[i], tap_angle[i], inv_tap[i], tap_real[i], tap_imag[i] =
            _tap_params(tap[i], shift[i])
    end

    active = Int[]
    sizehint!(active, n)
    @inbounds for i in 1:n
        status[i] == 1 && push!(active, i)
    end

    return (
        n = n,
        s_base = case_raw.baseMVA,
        from = from,
        to = to,
        status = status,
        active = active,
        r = r,
        x = x,
        b_sh = branch_raw.br_b,
        rate_a = branch_raw.rate_a,
        rate_b = branch_raw.rate_b,
        rate_c = branch_raw.rate_c,
        tap = tap,
        shift = shift,
        ang_min = branch_raw.angmin,
        ang_max = branch_raw.angmax,
        tap_position = tap_position,
        tap_angle = tap_angle,
        inv_tap = inv_tap,
        g = g,
        b = b,
        tap_real = tap_real,
        tap_imag = tap_imag,
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
        params.active,
        params.r,
        params.x,
        params.b_sh,
        params.rate_a,
        params.rate_b,
        params.rate_c,
        params.tap,
        params.shift,
        params.ang_min,
        params.ang_max,
        params.tap_position,
        params.tap_angle,
        params.inv_tap,
        params.g,
        params.b,
        params.tap_real,
        params.tap_imag,
    )
end
