# ======================================================================
# Distributed Workspace Contract
# ======================================================================
#
# `AbstractADMMWorkspace` is declared in src/types.jl and
# `AbstractALADINWorkspace` in aladin/core/aladin_support.jl. Phase 1
# keeps both in place as-is; the unification into a single
# `AbstractDistributedWorkspace` root is deferred to the workspace
# parametrization phase so we avoid cascading changes this early in
# the refactor.
#
# This file exists purely to document, in one place, the method
# contract every concrete workspace must satisfy. Think of it as a
# README-as-code: the comments here are authoritative.
#
# ----------------------------------------------------------------------
# ADMM WORKSPACE CONTRACT
# ----------------------------------------------------------------------
# For any concrete `ws <: AbstractADMMWorkspace`:
#
#   local_solve!(ws)
#   save_prev!(ws)
#   consensus!(ws)
#   dual_update!(ws)
#   compute_residuals!(ws)
#   compute_tolerances!(ws, eps_abs::Float64, eps_rel::Float64)
#   stop_if_converged(ws) :: Bool
#   rho_update!(ws, params::ADMMParams)
#   refresh_rho!(ws)
#   _record_iter!(ws)
#
# ----------------------------------------------------------------------
# ALADIN WORKSPACE CONTRACT
# ----------------------------------------------------------------------
# For any concrete `ws <: AbstractALADINWorkspace`:
#
#   aladin_local_solve!(ws)
#   solve_coordinator_qp!(ws) :: NamedTuple      # returns `nu`
#   primal_update!(ws)
#   dual_update!(ws, nu::NamedTuple)
#   compute_residuals!(ws)
#   compute_tolerances!(ws, eps_abs::Float64, eps_rel::Float64)
#   stop_if_converged(ws) :: Bool
#   adapt_rho_mu!(ws, params::ALADINParams)
#   _record_iter!(ws)
#
# Both families additionally need a `formulation(ws)` and
# `coupling(ws)` accessor once Phase 4/5 parametrize the structs.
# ======================================================================
