# ============================================================
# ALADIN-DC — Zone-level data layout & builder
# ============================================================

# ============================================================
# Struct
# ============================================================

mutable struct AladinZoneData <: AbstractAladinZoneData
    # Local subproblem (LP in DC, NLP in AC, SOCP in SOCWR) solved by Ipopt at
    # each ALADIN iteration. Persistent across iterations — only its objective
    # is updated by update_objective!. Distinct from the coordinator QP
    # (created on-the-fly in _build_qp_skeleton and solved by OSQP).
    local_model::JuMP.Model

    # Flat variable layout (all JuMP variables of this zone's local model).
    all_vars::Vector{JuMP.VariableRef}
    n_vars::Int
    var_to_col::Dict{JuMP.VariableRef, Int}

    pg_cols::Vector{Int}
    pg_c2_sbase2::Vector{Float64}   # 2·c2·sbase² (diagonal of Hessian on pg)
    pg_c1_sbase::Vector{Float64}    # c1·sbase    (linear term in gradient)

    # Shared-bus θ columns, aligned with part.zone_shared_bus_ids.
    shared_theta_cols::Vector{Int}

    edge_from_cols::Vector{Int}
    edge_to_cols::Vector{Int}
    edge_b::Vector{Float64}

    # Per-variable box bounds (NaN if absent).
    lb::Vector{Float64}
    ub::Vector{Float64}

    # Previous iterate (proximal center) and current solve result.
    x_prev::Vector{Float64}
    y::Vector{Float64}

    g::Vector{Float64}

    hess_diag::Vector{Float64}

    # Equality-constraint Jacobian (constant for linear DC equalities).
    eq_jac::SparseArrays.SparseMatrixCSC{Float64, Int}

    # Coordinator QP output (Δy correction in full variable space).
    delta::Vector{Float64}

    # Active-set detection buffer: columns at box bounds at current y.
    active_bound_cols::Vector{Int}
    quad_static::JuMP.QuadExpr
end

# ============================================================
# Variable / cost / topology collectors
# ============================================================

# ------------------------------------------------------------
# Generator cost layout
# ------------------------------------------------------------

function _collect_pg_layout(
    model::JuMP.Model,
    var_to_col::Dict{JuMP.VariableRef, Int},
    obj::ObjectiveData,
)
    pg_cols = Int[]
    pg_c2_sbase2 = Float64[]
    pg_c1_sbase = Float64[]

    haskey(model, :p_generation) || return pg_cols, pg_c2_sbase2, pg_c1_sbase

    pg_var = model[:p_generation]
    sbase = obj.sbase

    @inbounds for g in 1:obj.number_of_generators
        obj.generator_status[g] == 1 || continue
        v = pg_var[g, 1]
        col = get(var_to_col, v, 0)
        col == 0 && continue
        push!(pg_cols, col)

        coefs = obj.generator_costs[1][g]
        @assert length(coefs) <= 3 "ALADIN: generator cost must be at most quadratic (got degree $(length(coefs)-1))"
        c2 = (length(coefs) >= 3) ? coefs[end-2] : 0.0
        c1 = (length(coefs) >= 2) ? coefs[end-1] : 0.0

        push!(pg_c2_sbase2, 2.0 * c2 * sbase^2)
        push!(pg_c1_sbase, c1 * sbase)
    end

    return pg_cols, pg_c2_sbase2, pg_c1_sbase
end

# ------------------------------------------------------------
# Shared θ columns
# ------------------------------------------------------------

function _collect_shared_theta_cols(
    model::JuMP.Model,
    var_to_col::Dict{JuMP.VariableRef, Int},
    shared_bus_ids::Vector{Int},
)
    out = Vector{Int}(undef, length(shared_bus_ids))
    theta_var = model[:theta]
    @inbounds for i in eachindex(shared_bus_ids)
        v = theta_var[shared_bus_ids[i], 1]
        out[i] = var_to_col[v]
    end
    return out
end

# ------------------------------------------------------------
# Interzonal edge layout
# ------------------------------------------------------------

function _collect_edge_layout(
    model::JuMP.Model,
    var_to_col::Dict{JuMP.VariableRef, Int},
    part::PartitionDC;
    include_edges::Bool,
)
    include_edges || return Int[], Int[], Float64[]

    zone_edges = part.zone_edges
    n = length(zone_edges)

    from_cols = Vector{Int}(undef, n)
    to_cols = Vector{Int}(undef, n)
    b_vals = Vector{Float64}(undef, n)

    branch_map = part.branch_idx_full_to_zone
    branch = part.branch
    theta_var = model[:theta]

    @inbounds for (k, e) in enumerate(zone_edges)
        l = branch_map[e]
        from_b = branch.from[l]
        to_b = branch.to[l]
        from_cols[k] = var_to_col[theta_var[from_b, 1]]
        to_cols[k] = var_to_col[theta_var[to_b, 1]]
        b_vals[k] = branch.b[l]
    end

    return from_cols, to_cols, b_vals
end

# ============================================================
# Builder
# ============================================================

function build_aladin_zone_data(
    part::PartitionDC;
    include_edges::Bool = false,
)::AladinZoneData
    obj = part.objective_params
    model = part.local_model::JuMP.Model

    all_vars, n_vars, var_to_col = _build_var_layout(model)

    pg_cols, pg_c2_sbase2, pg_c1_sbase = _collect_pg_layout(model, var_to_col, obj)
    shared_theta_cols = _collect_shared_theta_cols(model, var_to_col, part.zone_shared_bus_ids)
    edge_from_cols, edge_to_cols, edge_b = _collect_edge_layout(model, var_to_col, part; include_edges)
    lb, ub = _collect_bounds(all_vars)
    
    hess_diag = fill(1e-6, n_vars)
    @inbounds for (k, col) in enumerate(pg_cols)
        hess_diag[col] += pg_c2_sbase2[k]
    end
    @assert all(isfinite, hess_diag) "DC Hessian diagonal has non-finite entries"
    eq_jac = _extract_equality_jacobian(model, var_to_col, n_vars)

    return AladinZoneData(
        model,
        all_vars,
        n_vars,
        var_to_col,
        pg_cols,
        pg_c2_sbase2,
        pg_c1_sbase,
        shared_theta_cols,
        edge_from_cols,
        edge_to_cols,
        edge_b,
        lb,
        ub,
        _init_from_start_values(all_vars),  # x_prev
        zeros(Float64, n_vars),  # y
        zeros(Float64, n_vars),  # g
        hess_diag,
        eq_jac,
        zeros(Float64, n_vars),  # delta
        Int[],                   # active_bound_cols
        zero(JuMP.QuadExpr),     # quad_static (filled by _build_quad_static_dc!)
    )
end
