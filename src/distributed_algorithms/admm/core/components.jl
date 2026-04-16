mutable struct VarState
    consensus::Vector{Float64}
    consensus_prev::Vector{Float64}
    duals::Matrix{Float64}
    rho::Base.RefValue{Float64}
    prev_rho::Base.RefValue{Float64}
end

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
        Ref(rho_init),
    )
end

# ============================================================
# BranchPair — for DummyBus (Method 2) where zone A's branch
# (i → d) must have the same p_flow (and q_flow / R / T) as
# zone B's branch (d → j). Unlike interzonal branches, the
# two branches have DIFFERENT global indices.
# ============================================================

struct BranchPair
    branch_a_global::Int   # global branch idx (in zone A's and system's layout)
    zone_a::Int            # 1-based position in partitions.parts
    branch_b_global::Int   # global branch idx for the B-side half
    zone_b::Int            # 1-based position in partitions.parts
end
