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
end

# ------------------------------------------------------------
# Variable layout & gradient/Hessian data
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# Build AladinZoneData from a PartitionDC (post-build_local_model!)
# ------------------------------------------------------------

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

    # DC QP Hessian: cost (2·c2·sbase² on pg columns) plus a uniform
    # Levenberg-Marquardt floor δ on every column. The floor satisfies the
    # H_i + ρΣ_i ≻ 0 hypothesis of HFD Theorem 4 cheaply and matches the
    # convention used in SOCWR (`_refresh_hess_diag!`) and AC
    # (`refresh_sensitivities_ac!`). The Phase-1 NLP still applies the full
    # (ρ/2)·σ shift over primal coordinates via `_prox_all_quadexpr!`; the
    # δ here only regularises the QP geometry without folding the proximal
    # weights.
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
    )
end

# ------------------------------------------------------------
# Per-iteration kernels: extract y, compute g, detect active bounds
# ------------------------------------------------------------

function extract_solve_values!(zone::AladinZoneData, model::JuMP.Model)
    y = zone.y
    @inbounds for i in 1:zone.n_vars
        y[i] = JuMP.value(zone.all_vars[i])
    end
    return nothing
end

function compute_cost_gradient!(zone::AladinZoneData)
    g = zone.g
    fill!(g, 0.0)
    y = zone.y
    @inbounds for k in eachindex(zone.pg_cols)
        col = zone.pg_cols[k]
        g[col] = zone.pg_c2_sbase2[k] * y[col] + zone.pg_c1_sbase[k]
    end
    return nothing
end

function detect_active_bounds!(zone::AladinZoneData; tol::Float64 = 1e-6)
    active = zone.active_bound_cols
    empty!(active)
    y = zone.y
    @inbounds for i in 1:zone.n_vars
        lo = zone.lb[i]
        hi = zone.ub[i]
        v = y[i]
        if !isnan(lo) && v - lo < tol
            push!(active, i)
        elseif !isnan(hi) && hi - v < tol
            push!(active, i)
        end
    end
    return nothing
end

function build_active_jacobian(
    zone::AladinZoneData;
    dual_tol::Float64 = 1e-6,
    dep_tol::Float64 = 1e-10,
)::SparseArrays.SparseMatrixCSC{Float64, Int}
    model = zone.model
    var_to_col = zone.var_to_col
    n_vars = zone.n_vars

    # Start with the precomputed equality rows.
    C_eq = zone.C_eq
    eq_nrows = size(C_eq, 1)

    # Collect active inequality rows by walking constraint types.
    ineq_rows = Int[]
    ineq_cols = Int[]
    ineq_vals = Vector{Float64}()
    row_offset = 0

    has_duals = JuMP.has_duals(model)

    for (F, S) in JuMP.list_of_constraint_types(model)
        F <: JuMP.VariableRef && continue
        S === MOI.EqualTo{Float64} && continue  # already in C_eq

        for ci in JuMP.all_constraints(model, F, S)
            if has_duals
                dval = try
                    JuMP.dual(ci)
                catch
                    0.0
                end
                abs(dval) > dual_tol || continue
            else
                continue
            end

            cobj = JuMP.constraint_object(ci)
            row_offset += 1
            _push_affine_row_offset!(
                ineq_rows, ineq_cols, ineq_vals,
                row_offset, cobj.func, var_to_col,
            )
        end
    end

    C_stacked = if row_offset == 0
        C_eq
    else
        C_ineq =
            SparseArrays.sparse(ineq_rows, ineq_cols, ineq_vals, row_offset, n_vars)
        vcat(C_eq, C_ineq)
    end

    return _delete_dependent_rows(C_stacked; tol = dep_tol)
end

function _delete_dependent_rows(
    C::SparseArrays.SparseMatrixCSC{Float64, Int};
    tol::Float64 = 1e-10,
)::SparseArrays.SparseMatrixCSC{Float64, Int}
    n_rows = size(C, 1)
    n_rows <= 1 && return C

    # Dense transpose: columns of A correspond to rows of C.
    # Sparse QR with column pivoting is not available in Julia stdlib
    # (`LinearAlgebra.qr(::SparseMatrixCSC, ::ColumnNorm)` is unsupported);
    # for robust LICQ rank detection we therefore densify here. This is
    # O(n_rows · n_vars) memory and O(n_rows²) time per call. Acceptable
    # for case ≤ 200 buses; for larger cases this is a known performance
    # ceiling and could be replaced by a sparse rank-revealing factorisation
    # (SuiteSparse SPQR or a streaming QR) as future work.
    A = Matrix(C')
    if size(A, 1) == 0 || size(A, 2) == 0
        return C
    end

    F = LinearAlgebra.qr(A, LinearAlgebra.ColumnNorm())

    # Count diagonal entries of R above tolerance → rank.
    r_diag_len = min(size(F.R, 1), size(F.R, 2))
    rnk = 0
    @inbounds for i in 1:r_diag_len
        if abs(F.R[i, i]) > tol
            rnk += 1
        end
    end

    rnk == n_rows && return C
    keep = sort(F.p[1:rnk])
    return C[keep, :]
end

@inline function _push_affine_row_offset!(
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

@inline _push_affine_row_offset!(
    ::Vector{Int},
    ::Vector{Int},
    ::Vector{Float64},
    ::Int,
    ::Any,
    ::Dict,
) = nothing

@inline refresh_lagrangian_hessian!(::AladinZoneData, ::AbstractALADINWorkspace) =
    nothing
