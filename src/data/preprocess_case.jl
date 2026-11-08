# ============================================================
# Index maps
# ============================================================

function _build_bus_idx_map(bus_ids::AbstractVector{<:Integer})::Vector{Int}
    n = length(bus_ids)
    n == 0 && return Int[]
    max_id = 0
    @inbounds for i in 1:n
        id = Int(bus_ids[i])
        id > max_id && (max_id = id)
    end
    v = zeros(Int, max_id)
    @inbounds for i in 1:n
        v[Int(bus_ids[i])] = i
    end
    return v
end

function _group_items_by_bus(
    item_bus::AbstractVector{<:Integer},
    bus_idx_map::Vector{Int},
    nbus::Int,
)::Vector{Vector{Int}}
    result = [Int[] for _ in 1:nbus]
    @inbounds for i in eachindex(item_bus)
        b = Int(item_bus[i])
        idx = b > length(bus_idx_map) ? 0 : bus_idx_map[b]
        idx == 0 && continue
        push!(result[idx], i)
    end
    return result
end

function add_case_index_maps!(raw_case::CaseRaw)
    m = raw_case.maps

    nbus = raw_case.bus === nothing ? 0 : length(raw_case.bus.bus_i)

    if raw_case.bus !== nothing
        m.bus_idx_map = _build_bus_idx_map(raw_case.bus.bus_i)
    else
        m.bus_idx_map = Int[]
    end

    if raw_case.gen !== nothing
        gen_bus = raw_case.gen.gen_bus
        resize!(m.gen_idx_map, length(gen_bus))
        @inbounds for (i, b) in enumerate(gen_bus)
            m.gen_idx_map[i] = Int(b)
        end
        m.bus_to_gens = _group_items_by_bus(gen_bus, m.bus_idx_map, nbus)
    else
        empty!(m.gen_idx_map)
        m.bus_to_gens = [Int[] for _ in 1:nbus]
    end

    if raw_case.branch !== nothing
        m.bus_to_branches_from =
            _group_items_by_bus(raw_case.branch.f_bus, m.bus_idx_map, nbus)
        m.bus_to_branches_to =
            _group_items_by_bus(raw_case.branch.t_bus, m.bus_idx_map, nbus)
    else
        m.bus_to_branches_from = [Int[] for _ in 1:nbus]
        m.bus_to_branches_to = [Int[] for _ in 1:nbus]
    end

    m.bus_to_storages =
        raw_case.storage === nothing ? [Int[] for _ in 1:nbus] :
        _group_items_by_bus(raw_case.storage.storage_bus, m.bus_idx_map, nbus)

    return raw_case
end

# ============================================================
# Angle conversion
# ============================================================

function convert_angles_to_rad!(raw_case::CaseRaw)
    deg2rad = π / 180

    if raw_case.bus !== nothing
        va = raw_case.bus.va
        @inbounds for i in eachindex(va)
            va[i] *= deg2rad
        end
    end

    if raw_case.branch !== nothing
        br = raw_case.branch
        vshift = br.shift
        v_min = br.angmin
        v_max = br.angmax

        @inbounds for i in eachindex(vshift)
            vshift[i] *= deg2rad
        end
        @inbounds for i in eachindex(v_min)
            v_min[i] *= deg2rad
        end
        @inbounds for i in eachindex(v_max)
            v_max[i] *= deg2rad
        end
    end

    return raw_case
end

# ============================================================
# Angle-difference clamping
# ============================================================

function clamp_voltage_angle_differences!(
    raw_case::CaseRaw,
    default_pad::Float64 = 1.0472,
)
    raw_case.branch === nothing && return raw_case
    amin = raw_case.branch.angmin
    amax = raw_case.branch.angmax
    half_pi = π / 2

    @inbounds for i in eachindex(amin)
        if amin[i] <= -half_pi
            amin[i] = -default_pad
        end
        if amax[i] >= half_pi
            amax[i] = default_pad
        end
        if amin[i] == 0.0 && amax[i] == 0.0
            amin[i] = -default_pad
            amax[i] = default_pad
        end
    end
    return raw_case
end

# ============================================================
# Per-unit conversion
# ============================================================

function convert_to_pu!(raw_case::CaseRaw)
    Sbase = raw_case.baseMVA
    Sbase > 0 || throw(
        ArgumentError(
            "baseMVA must be positive (got $(Sbase)); cannot convert case to per-unit.",
        ),
    )

    if raw_case.bus !== nothing
        b = raw_case.bus
        b.pd ./= Sbase
        b.qd ./= Sbase
        b.gs ./= Sbase
        b.bs ./= Sbase
    end

    if raw_case.gen !== nothing
        g = raw_case.gen
        g.pg ./= Sbase
        g.qg ./= Sbase
        g.pmax ./= Sbase
        g.pmin ./= Sbase
        g.qmax ./= Sbase
        g.qmin ./= Sbase
    end

    if raw_case.branch !== nothing
        br = raw_case.branch
        br.rate_a ./= Sbase
        br.rate_b ./= Sbase
        br.rate_c ./= Sbase
    end

    if raw_case.storage !== nothing
        s = raw_case.storage
        s.ps ./= Sbase
        s.qs ./= Sbase
        s.energy ./= Sbase
        s.energy_rating ./= Sbase
        s.charge_rating ./= Sbase
        s.discharge_rating ./= Sbase
        s.thermal_rating ./= Sbase
        s.qmin ./= Sbase
        s.qmax ./= Sbase
        s.q_loss ./= Sbase
        s.p_loss ./= Sbase
    end

    return raw_case
end

function convert_timeseries_to_pu!(raw_case::CaseRaw)
    Sbase = raw_case.baseMVA
    Sbase > 0 || throw(
        ArgumentError(
            "baseMVA must be positive (got $(Sbase)); cannot convert time series to per-unit.",
        ),
    )

    if raw_case.bus_timeseries !== nothing
        ts = raw_case.bus_timeseries
        ts.pd ./= Sbase
        ts.qd ./= Sbase
        ts.gs ./= Sbase
        ts.bs ./= Sbase
    end

    if raw_case.gen_timeseries !== nothing
        ts = raw_case.gen_timeseries
        _rescale_finite!(ts.pmax, Sbase)
        _rescale_finite!(ts.pmin, Sbase)
        _rescale_finite!(ts.qmax, Sbase)
        _rescale_finite!(ts.qmin, Sbase)
    end

    return raw_case
end

@inline function _rescale_finite!(v::Vector{Float64}, s::Float64)
    @inbounds for i in eachindex(v)
        isnan(v[i]) || (v[i] /= s)
    end
    return v
end

# ============================================================
# Main entry point
# ============================================================

function preprocess_case_data!(
    raw_case::CaseRaw,
    ::Type{T} = OPFModel,
) where {T <: OPFModel}
    if is_centralized_model(T)
        all_time_steps = get_model_times(raw_case)
        time_horizon = isempty(all_time_steps) ? 1 : maximum(all_time_steps)
        raw_case.T = time_horizon
    end

    time_horizon = raw_case.T

    check_time_alignment(raw_case)

    ensure_bus_timeseries!(raw_case, time_horizon)
    ensure_gencost_timeseries!(raw_case, time_horizon)

    extend_bus_timeseries!(raw_case, time_horizon)
    extend_gencost_timeseries!(raw_case, time_horizon)
    extend_gen_timeseries!(raw_case, time_horizon)

    sort_timeseries_by_time_and_id!(raw_case)

    convert_angles_to_rad!(raw_case)
    clamp_voltage_angle_differences!(raw_case)
    add_case_index_maps!(raw_case)

    convert_to_pu!(raw_case)
    convert_timeseries_to_pu!(raw_case)

    return nothing
end
