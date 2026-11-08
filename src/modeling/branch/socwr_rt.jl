# ============================================================
# Per-bus-pair voltage products (canonical, mirrors PowerModels)
# ============================================================

# Group active branches by their ordered (from, to) bus pair.
# Returns (pair_keys, pair_of, reps):
#   pair_keys :: Vector{Tuple{Int,Int}}  distinct (from,to) keys, in first-seen order
#   pair_of   :: Dict{Int,Int}           branch id -> pair index (into pair_keys)
#   reps      :: Vector{Int}             pair index -> representative branch id
function _socwr_buspair_groups(br::BranchData)
    active = br.active
    pair_index = Dict{Tuple{Int,Int},Int}()
    pair_keys = Tuple{Int,Int}[]
    reps = Int[]
    pair_of = Dict{Int,Int}()

    @inbounds for l in active
        key = (br.from[l], br.to[l])
        k = get(pair_index, key, 0)
        if k == 0
            push!(pair_keys, key)
            push!(reps, l)
            k = length(pair_keys)
            pair_index[key] = k
        end
        pair_of[l] = k
    end

    return pair_keys, pair_of, reps
end

# ============================================================
# Per-pair R/T product bounds
# ============================================================

# Mirrors PowerModels `ref_calc_voltage_product_bounds` on aggregated pair angle limits.
function _rt_pair_product_bounds(
    bus::BusData, br::BranchData,
    pair_keys::Vector{Tuple{Int,Int}}, pair_of::Dict{Int,Int},
)
    K = length(pair_keys)

    # Aggregated angle limits per pair (tightest across parallel branches).
    ang_min = fill(-Inf, K)
    ang_max = fill(Inf, K)
    @inbounds for l in br.active
        k = pair_of[l]
        amn = br.ang_min[l]
        amx = br.ang_max[l]
        isnan(amn) || (ang_min[k] = max(ang_min[k], amn))
        isnan(amx) || (ang_max[k] = min(ang_max[k], amx))
    end

    R_min = fill(-Inf, K)
    R_max = fill(Inf, K)
    T_min = fill(-Inf, K)
    T_max = fill(Inf, K)

    @inbounds for k in 1:K
        i, j = pair_keys[k]
        idx_i = bus.id_to_idx[i]
        idx_j = bus.id_to_idx[j]

        vm_fr_min = bus.v_min[idx_i]
        vm_fr_max = bus.v_max[idx_i]
        vm_to_min = bus.v_min[idx_j]
        vm_to_max = bus.v_max[idx_j]

        amin = ang_min[k]
        amax = ang_max[k]
        (isfinite(amin) && isfinite(amax)) || continue

        if amin >= 0
            R_max[k] = vm_fr_max * vm_to_max * cos(amin)
            R_min[k] = vm_fr_min * vm_to_min * cos(amax)
            T_max[k] = vm_fr_max * vm_to_max * sin(amax)
            T_min[k] = vm_fr_min * vm_to_min * sin(amin)
        elseif amax <= 0
            R_max[k] = vm_fr_max * vm_to_max * cos(amax)
            R_min[k] = vm_fr_min * vm_to_min * cos(amin)
            T_max[k] = vm_fr_min * vm_to_min * sin(amax)
            T_min[k] = vm_fr_max * vm_to_max * sin(amin)
        else # amin < 0 < amax
            R_max[k] = vm_fr_max * vm_to_max * 1.0
            R_min[k] = vm_fr_min * vm_to_min * min(cos(amin), cos(amax))
            T_max[k] = vm_fr_max * vm_to_max * sin(amax)
            T_min[k] = vm_fr_max * vm_to_max * sin(amin)
        end
    end

    return R_min, R_max, T_min, T_max
end

# ============================================================
# R/T variables with bounds (per-bus-pair, branch-indexed views)
# ============================================================

function add_rt_variables_with_bounds!(
    model::JuMP.Model,
    bus::BusData,
    br::BranchData;
    bounded::Bool = true,
)
    T = get(model.ext, :time_periods, 1)
    active = br.active

    pair_keys, pair_of, reps = _socwr_buspair_groups(br)
    K = length(pair_keys)

    # One product variable per bus pair (canonical), anonymous so the object
    # dictionary entries :R / :T_var below are free for the branch-indexed views.
    Rp = @variable(model, [k in 1:K, t in 1:T], base_name = "R_pair")
    Tp = @variable(model, [k in 1:K, t in 1:T], base_name = "T_pair")

    # Start values from the representative branch of each pair.
    vm_start = bus.vm_start
    va_start = bus.va_start
    shf = br.tap_angle
    id_to_idx = bus.id_to_idx

    @inbounds for k in 1:K, t in 1:T
        l = reps[k]
        idx_i = id_to_idx[br.from[l]]
        idx_j = id_to_idx[br.to[l]]

        Vm_i = vm_start[idx_i]
        Vm_j = vm_start[idx_j]
        Va_i = va_start[idx_i]
        Va_j = va_start[idx_j]

        if isfinite(Vm_i) && isfinite(Vm_j) &&
           isfinite(Va_i) && isfinite(Va_j) && isfinite(shf[l])
            set_start_value(Rp[k, t], Vm_i * Vm_j * cos(Va_i - Va_j - shf[l]))
            set_start_value(Tp[k, t], Vm_i * Vm_j * sin(Va_i - Va_j - shf[l]))
        end
    end

    if bounded
        R_min, R_max, T_min, T_max =
            _rt_pair_product_bounds(bus, br, pair_keys, pair_of)

        @inbounds for k in 1:K, t in 1:T
            isfinite(R_min[k]) && set_lower_bound(Rp[k, t], R_min[k])
            isfinite(R_max[k]) && set_upper_bound(Rp[k, t], R_max[k])
            isfinite(T_min[k]) && set_lower_bound(Tp[k, t], T_min[k])
            isfinite(T_max[k]) && set_upper_bound(Tp[k, t], T_max[k])
        end
    end

    # branch-indexed views aliasing the shared pair variables
    R_data = Dict{Tuple{Int,Int},JuMP.VariableRef}()
    T_data = Dict{Tuple{Int,Int},JuMP.VariableRef}()
    @inbounds for l in active, t in 1:T
        k = pair_of[l]
        R_data[(l, t)] = Rp[k, t]
        T_data[(l, t)] = Tp[k, t]
    end

    model[:R] = JuMP.Containers.SparseAxisArray(R_data)
    model[:T_var] = JuMP.Containers.SparseAxisArray(T_data)

    # Representative branch per pair, so callers can emit ONE Jabr cone per pair.
    model.ext[:socwr_pair_reps] = reps

    return nothing
end

# ============================================================
# Tangent angle cuts
# ============================================================

@inline function _tan_clip(x; eps = 1e-6)
    b = (pi / 2) - eps
    x = min(max(x, -b), b)
    return tan(x)
end

function add_rt_tangent_angle_cuts!(model::JuMP.Model, br::BranchData)
    T = get(model.ext, :time_periods, 1)
    active = br.active

    R = model[:R]
    T_var = model[:T_var]

    @constraint(model, [l in active, t in 1:T; !isnan(br.ang_max[l])],
        T_var[l, t] <= _tan_clip(br.ang_max[l]) * R[l, t]
    )

    @constraint(model, [l in active, t in 1:T; !isnan(br.ang_min[l])],
        T_var[l, t] >= _tan_clip(br.ang_min[l]) * R[l, t]
    )

    return nothing
end
