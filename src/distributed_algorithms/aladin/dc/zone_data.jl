# ============================================================
# ALADIN-DC — Zone-level data layout & builder
# ============================================================
#
# Defines `AladinZoneData` (the per-zone state container used by all
# Phase-1 / Phase-2 kernels) and the static initialisation that runs
# once at workspace build time:
#   - flat variable layout (`_build_var_layout`),
#   - generator cost coefficients (`_collect_pg_layout`),
#   - shared-bus θ columns (`_collect_shared_theta_cols`),
#   - interzonal edge layout (`_collect_edge_layout`),
#   - per-variable box bounds (`_collect_bounds`),
#   - linear equality Jacobian (`_extract_equality_jacobian`).
#
# Per-iteration kernels live in `sensitivities.jl`; active-set logic
# lives in `active_set.jl`.

# ============================================================
# Struct
# ============================================================

mutable struct AladinZoneData
    model::JuMP.Model

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
    C_eq::SparseArrays.SparseMatrixCSC{Float64, Int}

    # Coordinator QP output (Δy correction in full variable space).
    delta::Vector{Float64}

    # Active-set detection buffer: columns at box bounds at current y.
    active_bound_cols::Vector{Int}

    # Cached quadratic part of the Phase-1 proximal term:
    #   (ρ/2)·Σ σ_k·v_k²
    # Constant in (ρ, σ); rebuilt only when ρ changes (`adapt_rho_mu!`).
    # Mirrors `obj_static` in ADMM. The linear x_prev component is
    # rebuilt every iteration in `update_objective!`.
    quad_static::JuMP.QuadExpr
end

# ============================================================
# Variable / cost / topology collectors
# ============================================================

function _build_var_layout(model::JuMP.Model)
    all_vars = JuMP.all_variables(model)
    n = length(all_vars)
    var_to_col = Dict{JuMP.VariableRef, Int}()
    sizehint!(var_to_col, n)
    @inbounds for i in 1:n
        var_to_col[all_vars[i]] = i
    end
    return all_vars, n, var_to_col
end

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
        # The variable may not exist for filtered indices; skip gracefully.
        v = try
            pg_var[g, 1]
        catch
            continue
        end
        col = get(var_to_col, v, 0)
        col == 0 && continue
        push!(pg_cols, col)

        coefs =
            (g <= length(obj.generator_costs[1])) ? obj.generator_costs[1][g] : Float64[]
        # Matpower convention: cost_coefficients stored as [c_{n}, c_{n-1}, ..., c_0]
        # descending degrees. For quadratic (length=3): [c2, c1, c0].
        c2 = (length(coefs) >= 3) ? coefs[end-2] : 0.0
        c1 = (length(coefs) >= 2) ? coefs[end-1] : 0.0

        push!(pg_c2_sbase2, 2.0 * c2 * sbase^2)
        push!(pg_c1_sbase, c1 * sbase)
    end

    return pg_cols, pg_c2_sbase2, pg_c1_sbase
end

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

function _collect_edge_layout(
    model::JuMP.Model,
    var_to_col::Dict{JuMP.VariableRef, Int},
    part::PartitionDC,
)
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

function _collect_bounds(all_vars::Vector{JuMP.VariableRef})
    n = length(all_vars)
    lb = fill(NaN, n)
    ub = fill(NaN, n)
    @inbounds for i in 1:n
        v = all_vars[i]
        if JuMP.has_lower_bound(v)
            lb[i] = JuMP.lower_bound(v)
        end
        if JuMP.has_upper_bound(v)
            ub[i] = JuMP.upper_bound(v)
        end
    end
    return lb, ub
end

# ============================================================
# Equality Jacobian (constant — DC equalities are affine)
# ============================================================

function _extract_equality_jacobian(
    model::JuMP.Model,
    var_to_col::Dict{JuMP.VariableRef, Int},
    n_vars::Int,
)
    rows = Int[]
    cols = Int[]
    vals = Float64[]
    nrows = 0

    for (F, S) in JuMP.list_of_constraint_types(model)
        F <: JuMP.VariableRef && continue
        S === MOI.EqualTo{Float64} || continue

        for ci in JuMP.all_constraints(model, F, S)
            cobj = JuMP.constraint_object(ci)
            func = cobj.func
            nrows += 1
            _push_affine_row!(rows, cols, vals, nrows, func, var_to_col)
        end
    end

    return SparseArrays.sparse(rows, cols, vals, nrows, n_vars)
end

@inline function _push_affine_row!(
    rows::Vector{Int},
    cols::Vector{Int},
    vals::Vector{Float64},
    row::Int,
    func::JuMP.AffExpr,
    var_to_col::Dict{JuMP.VariableRef, Int},
)
    for (v, coef) in func.terms
        col = get(var_to_col, v, 0)
        col == 0 && continue
        push!(rows, row)
        push!(cols, col)
        push!(vals, coef)
    end
    return nothing
end

# Non-affine function — shouldn't happen in DC OPF, but guard anyway.
@inline function _push_affine_row!(
    ::Vector{Int},
    ::Vector{Int},
    ::Vector{Float64},
    ::Int,
    func,
    ::Dict,
)
    return error(
        "ALADIN-DC sensitivity extraction: expected affine equality constraint, " *
        "got $(typeof(func)). Canonical ALADIN-DC only supports linear local models.",
    )
end

# ============================================================
# Builder
# ============================================================

function build_aladin_zone_data(
    part::PartitionDC;
    include_edges::Bool = false,
)::AladinZoneData
    obj = part.objective_params
    model = part.model::JuMP.Model

    all_vars, n_vars, var_to_col = _build_var_layout(model)

    pg_cols, pg_c2_sbase2, pg_c1_sbase = _collect_pg_layout(model, var_to_col, obj)

    shared_theta_cols =
        _collect_shared_theta_cols(model, var_to_col, part.zone_shared_bus_ids)

    if include_edges
        edge_from_cols, edge_to_cols, edge_b =
            _collect_edge_layout(model, var_to_col, part)
    else
        edge_from_cols = Int[]
        edge_to_cols = Int[]
        edge_b = Float64[]
    end

    lb, ub = _collect_bounds(all_vars)
    hess_diag = fill(1e-6, n_vars)
    @inbounds for (k, col) in enumerate(pg_cols)
        hess_diag[col] += pg_c2_sbase2[k]
    end
    @assert all(isfinite, hess_diag) "DC Hessian diagonal has non-finite entries"
    C_eq = _extract_equality_jacobian(model, var_to_col, n_vars)

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
        C_eq,
        zeros(Float64, n_vars),  # delta
        Int[],                   # active_bound_cols
        zero(JuMP.QuadExpr),     # quad_static (filled by _build_quad_static_dc!)
    )
end
