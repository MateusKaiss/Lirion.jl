# ============================================================
# Flow expression accessors
# ------------------------------------------------------------
# Helpers that extract branch flow expressions (or their
# evaluated values) from the JuMP model. Used by consensus,
# dual updates, residuals, and objective builders.
# ============================================================

# ------------------------------------------------------------
# Active power flow (DC + generic)
# ------------------------------------------------------------

@inline _p_flow_from_exprs(model::JuMP.Model, l::Int, t::Int = 1) =
    model[:p_flow_branches][l, t]

@inline _p_flow_from_exprs(part::PartitionDC, l::Int; t::Int = 1) =
    _p_flow_from_exprs(part.model, l, t)

# ---------------------------------------------------------------------
# NOTE / TODO (DC flow helpers / indexing mismatch — investigate)
#
# Try to make DC follow the same pattern as AC/SOCP (use
# model[:p_flow_branches] via _p_flow_from_exprs), but keep
# _dc_flow_val_model for now to avoid wrong indexing/results to the paper.
# ---------------------------------------------------------------------

# Safe scalar evaluation of an AffExpr flow:
# - if model has values: value(p)
# - otherwise: fall back to `default` (does NOT use start_value on AffExpr)
@inline function _p_flow_value_or_zero(model::JuMP.Model, p; default::Float64 = 0.0)
    if JuMP.has_values(model)
        v = try
            JuMP.value(p)
        catch
            NaN
        end
        return isfinite(v) ? Float64(v) : default
    end
    return default
end

@inline function _dc_flow_val_model(model::JuMP.Model, part::PartitionDC, l::Int)
    p = _p_flow_from_exprs(part, l; t = 1)
    return _p_flow_value_or_zero(model, p; default = 0.0)
end

# ------------------------------------------------------------
# Active + reactive power flow (AC / SOCP)
# ------------------------------------------------------------

@inline _pq_flow_from_exprs(model::JuMP.Model, l::Int, t::Int = 1) =
    (model[:p_flow_branches][l, t], model[:q_flow_branches][l, t])

@inline _pq_flow_from_exprs(part::PartitionAC, l::Int; t::Int = 1) =
    _pq_flow_from_exprs(part.model, l, t)

@inline _pq_flow_from_exprs(part::PartitionSOCP, l::Int; t::Int = 1) =
    _pq_flow_from_exprs(part.model, l, t)

# ------------------------------------------------------------
# Second-order cone auxiliary variables (SOCP)
# ------------------------------------------------------------

@inline _rt_from_exprs(model::JuMP.Model, l::Int, t::Int = 1) =
    (model[:R][l, t], model[:T_var][l, t])

# ------------------------------------------------------------
# Evaluated flow values (used by consensus, dual updates, residuals)
# ------------------------------------------------------------

@inline _dc_flow_value(part::PartitionDC, l::Int) =
    _dc_flow_val_model(part.model, part, l)

@inline function _socp_pq_values(part::PartitionSOCP, l::Int)
    p, q = _pq_flow_from_exprs(part.model, l)
    return JuMP.value(p), JuMP.value(q)
end

@inline function _socp_rt_values(part::PartitionSOCP, l::Int)
    r, t = _rt_from_exprs(part.model, l)
    return JuMP.value(r), JuMP.value(t)
end
