# ======================================================================
# Formulation / Coupling Traits
# ======================================================================
#
# These traits are the single source of truth for *which primal
# variables participate in the consensus* for each (formulation,
# coupling) combination. All other code (dual updates, consensus
# averaging, residual computation, workspace field layout) consumes
# these traits instead of hard-coding per-variant dispatch.
#
# Canonical coupling vocabulary:
#
#   DC      ∋ :theta (bus angle), :p_flow (branch active flow)
#   AC      ∋ :theta, :v_mag, :p_flow, :q_flow
#   SOCWR   ∋ :v2 (|V|² relaxation), :p_flow, :q_flow, :r, :t
#
# The tuple order below is *load-bearing*: it dictates the storage
# order in DualState / ConsensusState and the iteration order in the
# unified dual_update! loop. Do not reorder without updating both.
# ======================================================================

# ----------------------------------------------------------------------
# consensus_vars(F, C) :: NTuple{N, Symbol}
# ----------------------------------------------------------------------
# The primal quantities shared across the partition boundary for each
# (formulation, coupling) pair. Each symbol identifies a *kind* of
# variable, not a storage field — the storage-name mapping lives in
# `admm_dual_field` / `aladin_dual_field` below.

@inline consensus_vars(::Type{DC}, ::Type{State}) = (:theta,)
@inline consensus_vars(::Type{DC}, ::Type{Power}) = (:p_flow,)
@inline consensus_vars(::Type{DC}, ::Type{Hybrid}) = (:theta, :p_flow)

@inline consensus_vars(::Type{AC}, ::Type{State}) = (:theta, :v_mag)
@inline consensus_vars(::Type{AC}, ::Type{Power}) = (:p_flow, :q_flow)
@inline consensus_vars(::Type{AC}, ::Type{Hybrid}) = (:theta, :v_mag, :p_flow, :q_flow)

@inline consensus_vars(::Type{SOCWR}, ::Type{State}) = (:v2, :r, :t)
@inline consensus_vars(::Type{SOCWR}, ::Type{Power}) = (:p_flow, :q_flow)
@inline consensus_vars(::Type{SOCWR}, ::Type{Hybrid}) = (:v2, :p_flow, :q_flow, :r, :t)

# Sugar: lift from a tag value.
@inline consensus_vars(tag::DistributedAlgorithm) =
    consensus_vars(formulation(tag), coupling(tag))

# ----------------------------------------------------------------------
# is_bus_var(sym) / is_branch_var(sym)
# ----------------------------------------------------------------------
# Classify a consensus variable by whether it lives on shared buses
# (state coupling) or on interzonal branches (power coupling). This is
# the single place to change the classification rule; consensus!,
# dual_update! and residuals! all dispatch off it.

@inline is_bus_var(sym::Symbol) = sym === :theta || sym === :v_mag || sym === :v2
@inline is_branch_var(sym::Symbol) =
    sym === :p_flow || sym === :q_flow || sym === :r || sym === :t

# ----------------------------------------------------------------------
# ADMM storage-name mapping
# ----------------------------------------------------------------------
# Maps a canonical variable symbol to the field names used by ADMM
# workspaces today. This is the only place that encodes the historical
# naming (`duals_v_mag` vs `duals_v2`, etc.). Future workspace
# parametrization (Phase 4) can drop this and store arrays in a
# homogeneous `NamedTuple` keyed on the canonical symbols.

@inline function admm_dual_field(sym::Symbol)
    sym === :theta && return :duals_theta
    sym === :v_mag && return :duals_v_mag
    sym === :v2 && return :duals_v2
    sym === :p_flow && return :duals_p_flow
    sym === :q_flow && return :duals_q_flow
    sym === :r && return :duals_r
    sym === :t && return :duals_t
    return error("admm_dual_field: unknown consensus variable $(sym)")
end

@inline function admm_consensus_field(sym::Symbol)
    sym === :theta && return :theta_consensus
    sym === :v_mag && return :v_mag_consensus
    sym === :v2 && return :v2_consensus
    sym === :p_flow && return :p_flow_consensus
    sym === :q_flow && return :q_flow_consensus
    sym === :r && return :r_consensus
    sym === :t && return :t_consensus
    return error("admm_consensus_field: unknown consensus variable $(sym)")
end

@inline function admm_rho_field(sym::Symbol)
    sym === :theta && return :rho_theta
    sym === :v_mag && return :rho_v
    sym === :v2 && return :rho_v2
    sym === :p_flow && return :rho_p
    sym === :q_flow && return :rho_q
    sym === :r && return :rho_r
    sym === :t && return :rho_t
    return error("admm_rho_field: unknown consensus variable $(sym)")
end

# ----------------------------------------------------------------------
# ALADIN storage-name mapping
# ----------------------------------------------------------------------
# ALADIN uses `lambda_*` instead of `duals_*` and `vm` instead of
# `v_mag` (following HFD 2016 notation).

@inline function aladin_dual_field(sym::Symbol)
    sym === :theta && return :lambda_theta
    sym === :v_mag && return :lambda_vm
    sym === :v2 && return :lambda_v2
    sym === :p_flow && return :lambda_p
    sym === :q_flow && return :lambda_q
    sym === :r && return :lambda_r
    sym === :t && return :lambda_t
    return error("aladin_dual_field: unknown consensus variable $(sym)")
end

@inline function aladin_dual_prev_field(sym::Symbol)
    sym === :theta && return :lambda_theta_prev
    sym === :v_mag && return :lambda_vm_prev
    sym === :v2 && return :lambda_v2_prev
    sym === :p_flow && return :lambda_p_prev
    sym === :q_flow && return :lambda_q_prev
    sym === :r && return :lambda_r_prev
    sym === :t && return :lambda_t_prev
    return error("aladin_dual_prev_field: unknown consensus variable $(sym)")
end

@inline function aladin_nu_field(sym::Symbol)
    sym === :theta && return :nu_theta
    sym === :v_mag && return :nu_vm
    sym === :v2 && return :nu_v2
    sym === :p_flow && return :nu_p
    sym === :q_flow && return :nu_q
    sym === :r && return :nu_r
    sym === :t && return :nu_t
    return error("aladin_nu_field: unknown consensus variable $(sym)")
end
