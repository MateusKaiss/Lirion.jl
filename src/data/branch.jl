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

    # Pi-model shunts per side (standard matpower: g_fr = g_to = 0, b_fr = b_to = b_sh/2)
    g_fr::Vector{Float64}
    b_fr::Vector{Float64}
    g_to::Vector{Float64}
    b_to::Vector{Float64}

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
    n = length(from)

    g = Vector{Float64}(undef, n)
    b = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        z2 = r[i]^2 + x[i]^2
        if z2 == 0.0
            if status[i] == 1
                error(
                    "Branch $i (from=$(from[i]), to=$(to[i])) has zero series impedance " *
                    "(r=0 and x=0). Admittance y=1/(r+jx) is undefined. " *
                    "Either remove the branch, disable it (br_status=0), or fix the input data.",
                )
            end
            g[i] = 0.0
            b[i] = 0.0
        else
            g[i] = r[i] / z2
            b[i] = -x[i] / z2
        end
    end

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

    # Pi-model shunts: standard matpower has g_fr = g_to = 0, b_fr = b_to = b_sh/2.
    b_sh = branch_raw.br_b
    g_fr = zeros(Float64, n)
    b_fr = b_sh ./ 2
    g_to = zeros(Float64, n)
    b_to = b_sh ./ 2

    return (
        n = n,
        s_base = case_raw.baseMVA,
        from = from,
        to = to,
        status = status,
        active = active,
        r = r,
        x = x,
        b_sh = b_sh,
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
        g_fr = g_fr,
        b_fr = b_fr,
        g_to = g_to,
        b_to = b_to,
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
        params.g_fr,
        params.b_fr,
        params.g_to,
        params.b_to,
        params.tap_real,
        params.tap_imag,
    )
end
