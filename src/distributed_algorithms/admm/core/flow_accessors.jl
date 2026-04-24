# ============================================================
# Expressions
# ============================================================

# ------------------------------------------------------------
# DC
# ------------------------------------------------------------

@inline _p_flow_from_exprs(model::JuMP.Model, l::Int, t::Int = 1) =
    model[:p_flow_branches][l, t]

@inline _p_flow_from_exprs(part::PartitionDC, l::Int; t::Int = 1) =
    _p_flow_from_exprs(part.model, l, t)

# ------------------------------------------------------------
# AC
# ------------------------------------------------------------

@inline _pq_flow_from_exprs(model::JuMP.Model, l::Int, t::Int = 1) =
    (model[:p_flow_branches][l, t], model[:q_flow_branches][l, t])

@inline _pq_flow_from_exprs(part::PartitionAC, l::Int; t::Int = 1) =
    _pq_flow_from_exprs(part.model, l, t)

# ------------------------------------------------------------
# SOCWR
# ------------------------------------------------------------

@inline _pq_flow_from_exprs(part::PartitionSOCWR, l::Int; t::Int = 1) =
    _pq_flow_from_exprs(part.model, l, t)

@inline _rt_from_exprs(model::JuMP.Model, l::Int, t::Int = 1) =
    (model[:R][l, t], model[:T_var][l, t])

# ============================================================
# Values
# ============================================================

# ------------------------------------------------------------
# DC
# ------------------------------------------------------------

@inline function _dc_flow_value(part::PartitionDC, l::Int)
    p = _p_flow_from_exprs(part, l)
    return JuMP.value(p)
end

# ------------------------------------------------------------
# SOCWR
# ------------------------------------------------------------

@inline function _socwr_pq_values(part::PartitionSOCWR, l::Int)
    p, q = _pq_flow_from_exprs(part.model, l)
    return JuMP.value(p), JuMP.value(q)
end

@inline function _socwr_p_value(part::PartitionSOCWR, l::Int)
    p, _ = _pq_flow_from_exprs(part.model, l)
    return JuMP.value(p)
end

@inline function _socwr_q_value(part::PartitionSOCWR, l::Int)
    _, q = _pq_flow_from_exprs(part.model, l)
    return JuMP.value(q)
end

@inline function _socwr_rt_values(part::PartitionSOCWR, l::Int)
    r, t = _rt_from_exprs(part.model, l)
    return JuMP.value(r), JuMP.value(t)
end
