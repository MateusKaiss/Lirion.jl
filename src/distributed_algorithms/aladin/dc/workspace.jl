struct SharedBusPair
    consensus_idx::Int   # index into sys.full_shared_bus_ids (global)
    bus_id::Int
    zone_i::Int          # position of zone i in partitions.parts
    col_i::Int           # column of θ[bus_id] in zone i's flat variable layout
    zone_j::Int
    col_j::Int
end

struct SharedEdgePair
    consensus_idx::Int   # index into sys.interzonal_branches
    edge_global::Int     # global branch index
    zone_i::Int
    edge_local_i::Int    # position of this edge in zone_i's zone_edges vector
    zone_j::Int
    edge_local_j::Int
end

# ------------------------------------------------------------
# DC State workspace
# ------------------------------------------------------------

Base.@kwdef mutable struct DCStateWorkspaceALADIN <: AbstractALADINWorkspace
    partitions::Partitions{PartitionDC}
    id::Vector{Int}

    # Per-zone canonical ALADIN data (variable layout, sensitivities, etc.)
    zones::Vector{AladinZoneData}

    # Pairwise bus consensus constraints (one per shared bus, for K=2).
    bus_pairs::Vector{SharedBusPair}

    # Duals keyed on canonical variable symbol. DC State couples only θ,
    # so `vars = (theta = LambdaState(n_pairs),)`. Legacy aliases
    # `ws.vars.theta.lambda` / `ws.vars.theta.lambda_prev` forward to
    # `vars.theta.lambda` / `vars.theta.lambda_prev` via the getproperty
    # layer in `aladin/core/components.jl`.
    vars::NamedTuple{(:theta,), Tuple{LambdaState}}

    # Phase-2 QP optimizer factory.
    qp_optimizer::Any

    # Parameters (held as Refs so adaptive schemes can update them).
    rho::Base.RefValue{Float64}
    mu::Base.RefValue{Float64}

    # Line-search parameters (default 1.0 → full step).
    alpha1::Float64 = 1.0
    alpha2::Float64 = 1.0
    alpha3::Float64 = 1.0

    # HFD Σ_i diagonal weights (per variable class).
    sigma_theta_shared::Float64 = 10.0
    sigma_theta_other::Float64 = 10.0
    sigma_pg::Float64 = 1.0
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
# DC Power workspace
# ------------------------------------------------------------

Base.@kwdef mutable struct DCPowerWorkspaceALADIN <: AbstractALADINWorkspace
    partitions::Partitions{PartitionDC}
    id::Vector{Int}
    interzonal::Vector{Int}

    zones::Vector{AladinZoneData}
    edge_pairs::Vector{SharedEdgePair}

    # Duals: DC Power couples only branch p_flow. Legacy alias
    # `ws.vars.p_flow.lambda` → `vars.p_flow.lambda`.
    vars::NamedTuple{(:p_flow,), Tuple{LambdaState}}

    qp_optimizer::Any

    rho::Base.RefValue{Float64}
    mu::Base.RefValue{Float64}

    alpha1::Float64 = 1.0
    alpha2::Float64 = 1.0
    alpha3::Float64 = 1.0

    sigma_theta_shared::Float64 = 10.0
    sigma_theta_other::Float64 = 10.0
    sigma_pg::Float64 = 1.0
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
# DC Hybrid workspace — both bus and edge couplings simultaneously
# ------------------------------------------------------------

Base.@kwdef mutable struct DCHybridWorkspaceALADIN <: AbstractALADINWorkspace
    partitions::Partitions{PartitionDC}
    id::Vector{Int}
    interzonal::Vector{Int}

    zones::Vector{AladinZoneData}
    bus_pairs::Vector{SharedBusPair}
    edge_pairs::Vector{SharedEdgePair}

    # Duals: DC Hybrid couples both bus θ and branch p_flow. Legacy
    # aliases `lambda_bus` → vars.theta, `lambda_edge` → vars.p_flow.
    vars::NamedTuple{(:theta, :p_flow), Tuple{LambdaState, LambdaState}}

    qp_optimizer::Any

    rho_bus::Base.RefValue{Float64}
    rho_edge::Base.RefValue{Float64}
    mu::Base.RefValue{Float64}

    alpha1::Float64 = 1.0
    alpha2::Float64 = 1.0
    alpha3::Float64 = 1.0

    sigma_theta_shared::Float64 = 10.0
    sigma_theta_other::Float64 = 10.0
    sigma_pg::Float64 = 1.0
    sigma_other::Float64 = 1.0

    # Per-component diagnostics (not part of the stopping rule)
    r_primal_bus::Float64 = 0.0
    r_primal_edge::Float64 = 0.0
    r_progress_bus::Float64 = 0.0
    r_progress_edge::Float64 = 0.0

    r_primal::Float64 = 0.0
    r_progress::Float64 = 0.0
    lambda_step::Float64 = 0.0
    eps_primal::Float64 = 0.0
    eps_progress::Float64 = 0.0
    iter::Int = 0
    hist::ALADINHistory = ALADINHistory()
end

# ------------------------------------------------------------
# Rho helpers (mirror ADMM's _get_rhos for history/logging)
# ------------------------------------------------------------

_get_rhos(ws::DCStateWorkspaceALADIN) = (rho = ws.rho[],)
_get_rhos(ws::DCPowerWorkspaceALADIN) = (rho = ws.rho[],)
_get_rhos(ws::DCHybridWorkspaceALADIN) =
    (rho_bus = ws.rho_bus[], rho_edge = ws.rho_edge[])
