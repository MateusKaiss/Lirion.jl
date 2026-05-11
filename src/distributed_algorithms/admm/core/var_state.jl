# ============================================================
# Struct
# ============================================================

mutable struct VarState
    consensus::Vector{Float64}
    consensus_prev::Vector{Float64}
    duals::Matrix{Float64}
    rho::Base.RefValue{Float64}
    prev_rho::Base.RefValue{Float64}
end

# ============================================================
# Constructor
# ============================================================

function VarState(
    n::Int,
    n_zones::Int,
    rho_init::Float64;
    init_value::Float64 = 0.0,
)
    consensus = fill(init_value, n)
    # `prev_rho` starts as NaN so refresh_rho! always rebuilds the
    # cached quadratics on the first call (NaN != anything, including
    # itself). Subsequent calls compare equal and short-circuit when ρ
    # is unchanged.
    return VarState(
        consensus,
        fill(init_value, n),
        zeros(n, n_zones),
        Ref(rho_init),
        Ref(NaN),
    )
end

# Twin VarState that shares its rho refs with another (used so p_dummy
# follows the same rho schedule as the tie-line p_flow it complements).
function VarState(
    n::Int,
    n_zones::Int,
    rho::Base.RefValue{Float64},
    prev_rho::Base.RefValue{Float64};
    init_value::Float64 = 0.0,
)
    return VarState(
        fill(init_value, n),
        fill(init_value, n),
        zeros(n, n_zones),
        rho,
        prev_rho,
    )
end
