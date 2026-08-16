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
    return VarState(
        consensus,
        fill(init_value, n),
        zeros(n, n_zones),
        Ref(rho_init),
        Ref(NaN),
    )
end

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
