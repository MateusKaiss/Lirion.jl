# ============================================================
# ALADIN-SOCWR — workspaces (State / Power / Hybrid)

# ------------------------------------------------------------
# Pair descriptors
# ------------------------------------------------------------

struct SharedBusPairSOCWR
    consensus_idx::Int
    bus_id::Int
    zone_i::Int
    v2_col_i::Int
    zone_j::Int
    v2_col_j::Int
end

struct SharedEdgePairSOCWR
    consensus_idx::Int
    edge_global::Int
    zone_i::Int
    edge_local_i::Int   # position in zone_i.edge_r_cols / edge_t_cols
    zone_j::Int
    edge_local_j::Int
end

# ------------------------------------------------------------
# STATE workspace
# ------------------------------------------------------------

Base.@kwdef mutable struct SOCWRStateWorkspaceALADIN <: AbstractALADINWorkspace
    partitions::Partitions{PartitionSOCWR}
    id::Vector{Int}
    interzonal::Vector{Int}

    zones::Vector{AladinZoneDataSOCWR}

    # State consensus on v² at shared buses and on R/T at interzonal
    # branches.
    bus_pairs::Vector{SharedBusPairSOCWR}
    edge_pairs::Vector{SharedEdgePairSOCWR}

    # Duals keyed on canonical variable symbol.
    vars::NamedTuple{(:v2, :r, :t), Tuple{LambdaState, LambdaState, LambdaState}}

    qp_optimizer::Any

    rho::Base.RefValue{Float64}
    mu::Base.RefValue{Float64}

    alpha1::Float64 = 1.0
    alpha2::Float64 = 1.0
    alpha3::Float64 = 1.0

    # Initial (cap) value of α₂ for the Armijo line search each iter.
    # Set from cfg.alpha2 at workspace build time; the effective α₂
    # used by primal_update! is stored in `alpha2` and refreshed each
    # iteration by `_line_search_step!`.
    alpha2_init::Float64 = 1.0
    line_search::Bool = true

    # Σ_i proximal weights per variable class. The shared state
    # (v², R, T) gets a stronger pull than free interior variables,
    # matching the structural sparsity of the Phase-2 consensus rows.
    sigma_v2::Float64 = 10.0
    sigma_r::Float64 = 5.0
    sigma_t::Float64 = 5.0
    sigma_theta::Float64 = 1.0
    sigma_pg::Float64 = 1.0
    sigma_qg::Float64 = 1.0
    sigma_other::Float64 = 1.0

    r_primal::Float64 = 0.0
    r_progress::Float64 = 0.0
    lambda_step::Float64 = 0.0
    eps_primal::Float64 = 0.0
    eps_progress::Float64 = 0.0
    iter::Int = 0
    hist::ALADINHistory = ALADINHistory()
end

# ------------------------------------------------------------
# POWER workspace
# ------------------------------------------------------------

Base.@kwdef mutable struct SOCWRPowerWorkspaceALADIN <: AbstractALADINWorkspace
    partitions::Partitions{PartitionSOCWR}
    id::Vector{Int}
    interzonal::Vector{Int}

    zones::Vector{AladinZoneDataSOCWR}

    edge_pairs::Vector{SharedEdgePairSOCWR}

    vars::NamedTuple{(:p_flow, :q_flow), Tuple{LambdaState, LambdaState}}

    qp_optimizer::Any

    rho::Base.RefValue{Float64}
    mu::Base.RefValue{Float64}

    alpha1::Float64 = 1.0
    alpha2::Float64 = 1.0
    alpha3::Float64 = 1.0

    # Initial (cap) value of α₂ for the Armijo line search each iter.
    # Set from cfg.alpha2 at workspace build time; the effective α₂
    # used by primal_update! is stored in `alpha2` and refreshed each
    # iteration by `_line_search_step!`.
    alpha2_init::Float64 = 1.0
    line_search::Bool = true

    sigma_v2::Float64 = 5.0
    sigma_r::Float64 = 5.0
    sigma_t::Float64 = 5.0
    sigma_theta::Float64 = 1.0
    sigma_pg::Float64 = 1.0
    sigma_qg::Float64 = 1.0
    sigma_other::Float64 = 1.0

    r_primal::Float64 = 0.0
    r_progress::Float64 = 0.0
    lambda_step::Float64 = 0.0
    eps_primal::Float64 = 0.0
    eps_progress::Float64 = 0.0
    iter::Int = 0
    hist::ALADINHistory = ALADINHistory()
end

# ------------------------------------------------------------
# HYBRID workspace
# ------------------------------------------------------------

Base.@kwdef mutable struct SOCWRHybridWorkspaceALADIN <: AbstractALADINWorkspace
    partitions::Partitions{PartitionSOCWR}
    id::Vector{Int}
    interzonal::Vector{Int}

    zones::Vector{AladinZoneDataSOCWR}

    bus_pairs::Vector{SharedBusPairSOCWR}
    edge_pairs::Vector{SharedEdgePairSOCWR}

    vars::NamedTuple{
        (:v2, :p_flow, :q_flow, :r, :t),
        Tuple{LambdaState, LambdaState, LambdaState, LambdaState, LambdaState},
    }

    qp_optimizer::Any

    rho::Base.RefValue{Float64}
    mu::Base.RefValue{Float64}

    alpha1::Float64 = 1.0
    alpha2::Float64 = 1.0
    alpha3::Float64 = 1.0

    # Initial (cap) value of α₂ for the Armijo line search each iter.
    # Set from cfg.alpha2 at workspace build time; the effective α₂
    # used by primal_update! is stored in `alpha2` and refreshed each
    # iteration by `_line_search_step!`.
    alpha2_init::Float64 = 1.0
    line_search::Bool = true

    sigma_v2::Float64 = 10.0
    sigma_r::Float64 = 5.0
    sigma_t::Float64 = 5.0
    sigma_theta::Float64 = 1.0
    sigma_pg::Float64 = 1.0
    sigma_qg::Float64 = 1.0
    sigma_other::Float64 = 1.0

    r_primal::Float64 = 0.0
    r_progress::Float64 = 0.0
    lambda_step::Float64 = 0.0
    eps_primal::Float64 = 0.0
    eps_progress::Float64 = 0.0
    iter::Int = 0
    hist::ALADINHistory = ALADINHistory()
end

# ------------------------------------------------------------
# Rho helpers
# ------------------------------------------------------------

_get_rhos(ws::SOCWRStateWorkspaceALADIN) = (rho = ws.rho[],)
_get_rhos(ws::SOCWRPowerWorkspaceALADIN) = (rho = ws.rho[],)
_get_rhos(ws::SOCWRHybridWorkspaceALADIN) = (rho = ws.rho[],)
