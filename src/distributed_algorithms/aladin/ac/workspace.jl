# ------------------------------------------------------------
# Pair descriptors
# ------------------------------------------------------------

struct SharedBusPairAC
    consensus_idx::Int
    bus_id::Int
    zone_i::Int
    theta_col_i::Int
    vm_col_i::Int
    zone_j::Int
    theta_col_j::Int
    vm_col_j::Int
end

struct SharedEdgePairAC
    consensus_idx::Int
    edge_global::Int
    zone_i::Int
    edge_local_i::Int
    zone_j::Int
    edge_local_j::Int
end

# ------------------------------------------------------------
# STATE workspace
# ------------------------------------------------------------

Base.@kwdef mutable struct ACStateWorkspaceALADIN <: AbstractALADINWorkspace
    partitions::Partitions{PartitionAC}
    id::Vector{Int}
    interzonal::Vector{Int}

    zones::Vector{AladinZoneDataAC}

    bus_pairs::Vector{SharedBusPairAC}

    # Duals keyed on canonical variable symbol. Legacy accessors
    # `ws.vars.theta.lambda`, `ws.vars.v_mag.lambda` forward to `vars.theta.lambda`
    # and `vars.v_mag.lambda`.
    vars::NamedTuple{(:theta, :v_mag), Tuple{LambdaState, LambdaState}}

    qp_optimizer::Any

    rho::Base.RefValue{Float64}
    mu::Base.RefValue{Float64}

    alpha1::Float64 = 1.0
    alpha2::Float64 = 1.0
    alpha3::Float64 = 1.0

    # Σ_i proximal weights per variable class. Defaults are tuned so
    # that tightly-coupled shared state (θ, v_m) feels a stronger
    # proximal pull than free interior variables — a standard
    # engineering choice that matches the structural sparsity of the
    # Phase-2 QP consensus rows.
    sigma_theta::Float64 = 10.0
    sigma_vm::Float64 = 5.0
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

Base.@kwdef mutable struct ACPowerWorkspaceALADIN <: AbstractALADINWorkspace
    partitions::Partitions{PartitionAC}
    id::Vector{Int}
    interzonal::Vector{Int}

    zones::Vector{AladinZoneDataAC}

    edge_pairs::Vector{SharedEdgePairAC}

    # Duals: branch p_flow + q_flow. Legacy aliases lambda_p/lambda_q.
    vars::NamedTuple{(:p_flow, :q_flow), Tuple{LambdaState, LambdaState}}

    qp_optimizer::Any

    rho::Base.RefValue{Float64}
    mu::Base.RefValue{Float64}

    alpha1::Float64 = 1.0
    alpha2::Float64 = 1.0
    alpha3::Float64 = 1.0

    sigma_theta::Float64 = 10.0
    sigma_vm::Float64 = 5.0
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

Base.@kwdef mutable struct ACHybridWorkspaceALADIN <: AbstractALADINWorkspace
    partitions::Partitions{PartitionAC}
    id::Vector{Int}
    interzonal::Vector{Int}

    zones::Vector{AladinZoneDataAC}

    bus_pairs::Vector{SharedBusPairAC}
    edge_pairs::Vector{SharedEdgePairAC}

    # Duals: bus (θ, |V|) + branch (p, q). Legacy accessors
    # `ws.vars.theta.lambda`, `ws.vars.v_mag.lambda`, `ws.vars.p_flow.lambda`, `ws.vars.q_flow.lambda`.
    vars::NamedTuple{
        (:theta, :v_mag, :p_flow, :q_flow),
        Tuple{LambdaState, LambdaState, LambdaState, LambdaState},
    }

    qp_optimizer::Any

    rho::Base.RefValue{Float64}
    mu::Base.RefValue{Float64}

    alpha1::Float64 = 1.0
    alpha2::Float64 = 1.0
    alpha3::Float64 = 1.0

    sigma_theta::Float64 = 10.0
    sigma_vm::Float64 = 5.0
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

_get_rhos(ws::ACStateWorkspaceALADIN) = (rho = ws.rho[],)
_get_rhos(ws::ACPowerWorkspaceALADIN) = (rho = ws.rho[],)
_get_rhos(ws::ACHybridWorkspaceALADIN) = (rho = ws.rho[],)
