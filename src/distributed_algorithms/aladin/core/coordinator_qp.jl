# ============================================================
# ALADIN — Phase 2 coordinator QP (shared infrastructure)
# ============================================================

# ============================================================
# QP skeleton
# ============================================================

function _build_qp_skeleton(ws::AbstractALADINWorkspace)
    model = JuMP.Model(ws.qp_optimizer)
    n_zones = length(ws.zones)

    dy_all = Vector{Vector{JuMP.VariableRef}}(undef, n_zones)
    @inbounds for i in 1:n_zones
        dy_all[i] =
            JuMP.@variable(model, [1:ws.zones[i].n_vars], base_name = "dy$(i)")
    end

    qp_obj = zero(JuMP.QuadExpr)
    @inbounds for i in 1:n_zones
        _add_zone_cost!(qp_obj, ws.zones[i], dy_all[i])
    end

    @inbounds for i in 1:n_zones
        _add_active_constraints!(model, ws.zones[i], dy_all[i])
        _add_active_bound_constraints!(model, ws.zones[i], dy_all[i])
    end

    return model, dy_all, qp_obj
end

# ============================================================
# Δy extraction
# ============================================================

function _extract_delta!(ws::AbstractALADINWorkspace, dy_all)
    @inbounds for i in eachindex(ws.zones)
        delta = ws.zones[i].delta
        v = dy_all[i]
        for k in 1:ws.zones[i].n_vars
            delta[k] = JuMP.value(v[k])
        end
    end
    return nothing
end

# ============================================================
# Zone cost — diagonal Hessian default
# ============================================================

@inline function _add_zone_cost!(
    qp_obj::JuMP.QuadExpr,
    zone::AbstractAladinZoneData,
    dy::Vector{JuMP.VariableRef},
)
    @inbounds for k in 1:zone.n_vars
        gk = zone.g[k]
        if gk != 0.0
            JuMP.add_to_expression!(qp_obj, gk, dy[k])
        end
        hk = zone.hess_diag[k]
        if hk != 0.0
            JuMP.add_to_expression!(qp_obj, 0.5 * hk, dy[k], dy[k])
        end
    end
    return nothing
end

# ============================================================
# Active set — equality rows + box bounds
# ============================================================

# ------------------------------------------------------------
# Active equality rows
# ------------------------------------------------------------

function _add_active_constraints!(
    model::JuMP.Model,
    zone,
    dy::Vector{JuMP.VariableRef},
)
    C = build_active_jacobian(zone)
    nrows = size(C, 1)
    if nrows == 0
        return nothing
    end
    rows = SparseArrays.rowvals(C)
    vals = SparseArrays.nonzeros(C)

    row_exprs = [zero(JuMP.AffExpr) for _ in 1:nrows]
    @inbounds for col in 1:size(C, 2)
        for idx in SparseArrays.nzrange(C, col)
            r = rows[idx]
            JuMP.add_to_expression!(row_exprs[r], vals[idx], dy[col])
        end
    end
    @inbounds for r in 1:nrows
        JuMP.@constraint(model, row_exprs[r] == 0)
    end
    return nothing
end

# ------------------------------------------------------------
# Active box bounds
# ------------------------------------------------------------

function _add_active_bound_constraints!(
    model::JuMP.Model,
    zone,
    dy::Vector{JuMP.VariableRef},
)
    @inbounds for k in zone.active_bound_cols
        JuMP.@constraint(model, dy[k] == 0)
    end
    return nothing
end

# ============================================================
# Termination guard
# ============================================================

function _check_qp_status!(model::JuMP.Model, tag::AbstractString)
    ts = JuMP.termination_status(model)
    acceptable =
        ts == MOI.OPTIMAL ||
        ts == MOI.LOCALLY_SOLVED ||
        ts == MOI.ALMOST_OPTIMAL ||
        ts == MOI.ALMOST_LOCALLY_SOLVED
    if !acceptable
        raw = try
            JuMP.raw_status(model)
        catch
            "raw_status unavailable"
        end
        error("""
        $tag failed
          termination_status = $ts
          raw_status         = $raw
        """)
    end
    return nothing
end
