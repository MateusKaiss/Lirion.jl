# ============================================================
# ALADIN-SOCWR — sensitivities (canonical Houska-Frasch-Diehl 2016)

mutable struct AladinZoneDataSOCWR
    # JuMP model reference — used every iteration to walk constraints
    # and extract active-inequality duals.
    model::JuMP.Model

    # Flat variable layout.
    all_vars::Vector{JuMP.VariableRef}
    n_vars::Int
    var_to_col::Dict{JuMP.VariableRef, Int}

    # Column classes for Σ_i sigma weights in the proximal term.
    theta_cols::Vector{Int}
    v2_cols::Vector{Int}
    r_cols::Vector{Int}
    t_cols::Vector{Int}
    pg_cols::Vector{Int}
    qg_cols::Vector{Int}
    curt_cols::Vector{Int}

    # Cost-support data on pg (quadratic cost in p_generation).
    pg_c2_sbase2::Vector{Float64}     # 2·c2·sbase²
    pg_c1_sbase::Vector{Float64}      # c1·sbase

    # Shared-bus v² column mapping (one entry per part.zone_shared_bus_ids).
    shared_v2_cols::Vector{Int}
    edge_r_cols::Vector{Int}
    edge_t_cols::Vector{Int}
    edge_flow_cols::Vector{Vector{Int}}
    edge_pflow_coefs::Vector{Vector{Float64}}
    edge_qflow_coefs::Vector{Vector{Float64}}
    edge_pflow_const::Vector{Float64}
    edge_qflow_const::Vector{Float64}

    # Per-variable box bounds (NaN when absent).
    lb::Vector{Float64}
    ub::Vector{Float64}

    # Previous iterate (proximal centre) and current Phase-1 solution.
    x_prev::Vector{Float64}
    y::Vector{Float64}

    # Cost gradient ∇f(y) at the current y (refreshed every iteration).
    g::Vector{Float64}

    # Lagrangian Hessian diagonal. Structurally constant (2·c2·sbase²
    # on pg cols) plus an iteration-dependent proximal shift ρ·σ added
    # by `_refresh_hess_diag!`.
    hess_diag::Vector{Float64}
    base_hess_diag::Vector{Float64}   # constant part, cached at build time

    # Constant equality Jacobian (cached once; SOCWR equalities are
    # all affine, so they never change across iterations).
    C_eq::SparseArrays.SparseMatrixCSC{Float64, Int}

    # Phase-2 QP correction Δy_i.
    delta::Vector{Float64}

    # Active box bounds at current y.
    active_bound_cols::Vector{Int}
end

# ------------------------------------------------------------
# Column-class detection
# ------------------------------------------------------------

function _collect_socwr_cols(
    model::JuMP.Model,
    var_to_col::Dict{JuMP.VariableRef, Int},
    obj::ObjectiveData,
)
    theta_cols = Int[]
    v2_cols = Int[]
    r_cols = Int[]
    t_cols = Int[]
    pg_cols = Int[]
    qg_cols = Int[]
    curt_cols = Int[]

    pg_c2_sbase2 = Float64[]
    pg_c1_sbase = Float64[]

    sbase = obj.sbase

    if haskey(model, :theta)
        for v in values(model[:theta].data)
            col = get(var_to_col, v, 0)
            col == 0 && continue
            push!(theta_cols, col)
        end
    end

    if haskey(model, :v_magnitude_sq)
        for v in values(model[:v_magnitude_sq].data)
            col = get(var_to_col, v, 0)
            col == 0 && continue
            push!(v2_cols, col)
        end
    end

    if haskey(model, :R)
        for v in values(model[:R].data)
            col = get(var_to_col, v, 0)
            col == 0 && continue
            push!(r_cols, col)
        end
    end

    if haskey(model, :T_var)
        for v in values(model[:T_var].data)
            col = get(var_to_col, v, 0)
            col == 0 && continue
            push!(t_cols, col)
        end
    end

    if haskey(model, :p_generation)
        pg_var = model[:p_generation]
        for g in 1:obj.number_of_generators
            obj.generator_status[g] == 1 || continue
            v = try
                pg_var[g, 1]
            catch
                continue
            end
            col = get(var_to_col, v, 0)
            col == 0 && continue
            push!(pg_cols, col)

            coefs =
                (g <= length(obj.generator_costs[1])) ? obj.generator_costs[1][g] :
                Float64[]
            # Matpower convention: descending degrees. For quadratic (len 3)
            # coefs = [c2, c1, c0].
            c2 = (length(coefs) >= 3) ? coefs[end-2] : 0.0
            c1 = (length(coefs) >= 2) ? coefs[end-1] : 0.0
            push!(pg_c2_sbase2, 2.0 * c2 * sbase^2)
            push!(pg_c1_sbase, c1 * sbase)
        end
    end

    if haskey(model, :q_generation)
        qg_var = model[:q_generation]
        for g in 1:obj.number_of_generators
            obj.generator_status[g] == 1 || continue
            v = try
                qg_var[g, 1]
            catch
                continue
            end
            col = get(var_to_col, v, 0)
            col == 0 && continue
            push!(qg_cols, col)
        end
    end

    if haskey(model, :curtailment)
        for v in values(model[:curtailment].data)
            col = get(var_to_col, v, 0)
            col == 0 && continue
            push!(curt_cols, col)
        end
    end

    return (
        theta_cols = theta_cols,
        v2_cols = v2_cols,
        r_cols = r_cols,
        t_cols = t_cols,
        pg_cols = pg_cols,
        qg_cols = qg_cols,
        curt_cols = curt_cols,
        pg_c2_sbase2 = pg_c2_sbase2,
        pg_c1_sbase = pg_c1_sbase,
    )
end

function _collect_shared_v2_cols(
    model::JuMP.Model,
    var_to_col::Dict{JuMP.VariableRef, Int},
    shared_bus_ids::Vector{Int},
)
    out = Vector{Int}(undef, length(shared_bus_ids))
    v2_var = model[:v_magnitude_sq]
    @inbounds for i in eachindex(shared_bus_ids)
        out[i] = var_to_col[v2_var[shared_bus_ids[i], 1]]
    end
    return out
end

function _collect_socwr_edge_rt_cols(
    model::JuMP.Model,
    var_to_col::Dict{JuMP.VariableRef, Int},
    part::PartitionSOCWR,
)
    zone_edges = part.zone_edges
    n = length(zone_edges)

    r_cols = Vector{Int}(undef, n)
    t_cols = Vector{Int}(undef, n)

    R = model[:R]
    Tv = model[:T_var]
    branch_map = part.branch_idx_full_to_zone

    @inbounds for (k, e) in enumerate(zone_edges)
        l = branch_map[e]
        r_cols[k] = var_to_col[R[l, 1]]
        t_cols[k] = var_to_col[Tv[l, 1]]
    end

    return r_cols, t_cols
end

# ------------------------------------------------------------
# Affine expansion of p_flow_branches[l,1] and q_flow_branches[l,1]
# ------------------------------------------------------------
# SOCWR branch flows (see src/modeling/branch/socwr.jl) are pure
# affine combinations of (v²[from], R[l], T[l]). Walking the
# @expression's .terms yields the exact coefficient vector of the
# flow w.r.t. the zone's flat variable layout, which is everything
# we need for the Power/Hybrid coordinator QP consensus rows.

function _collect_socwr_edge_flow_expansion(
    model::JuMP.Model,
    var_to_col::Dict{JuMP.VariableRef, Int},
    part::PartitionSOCWR,
)
    zone_edges = part.zone_edges
    n = length(zone_edges)

    cols = Vector{Vector{Int}}(undef, n)
    coefs_p = Vector{Vector{Float64}}(undef, n)
    coefs_q = Vector{Vector{Float64}}(undef, n)
    const_p = Vector{Float64}(undef, n)
    const_q = Vector{Float64}(undef, n)

    pflow = model[:p_flow_branches]
    qflow = model[:q_flow_branches]
    branch_map = part.branch_idx_full_to_zone

    @inbounds for (k, e) in enumerate(zone_edges)
        l = branch_map[e]
        p_expr = pflow[l, 1]
        q_expr = qflow[l, 1]

        # p_flow and q_flow reference the same set of variables
        # (v²[from], R[l], T[l]) — walk p's terms to fix the col order
        # and look up the matching coefficient in q.
        p_terms = p_expr.terms
        q_terms = q_expr.terms
        k_cols = Vector{Int}(undef, length(p_terms))
        k_coefs_p = Vector{Float64}(undef, length(p_terms))
        k_coefs_q = Vector{Float64}(undef, length(p_terms))

        i = 0
        for (v, cp) in p_terms
            col = get(var_to_col, v, 0)
            col == 0 && continue
            i += 1
            k_cols[i] = col
            k_coefs_p[i] = cp
            k_coefs_q[i] = get(q_terms, v, 0.0)
        end
        resize!(k_cols, i)
        resize!(k_coefs_p, i)
        resize!(k_coefs_q, i)

        # Any q variable that was not in p's terms (unlikely for
        # SOCWR, but safe) is appended here.
        for (v, cq) in q_terms
            col = get(var_to_col, v, 0)
            col == 0 && continue
            if !(col in k_cols)
                push!(k_cols, col)
                push!(k_coefs_p, 0.0)
                push!(k_coefs_q, cq)
            end
        end

        cols[k] = k_cols
        coefs_p[k] = k_coefs_p
        coefs_q[k] = k_coefs_q
        const_p[k] = p_expr.constant
        const_q[k] = q_expr.constant
    end

    return cols, coefs_p, coefs_q, const_p, const_q
end

# ------------------------------------------------------------
# Build AladinZoneDataSOCWR from a PartitionSOCWR
# ------------------------------------------------------------

function build_aladin_zone_data_socwr(
    part::PartitionSOCWR;
    include_edges::Bool = true,
)::AladinZoneDataSOCWR
    obj = part.objective_params
    model = part.model::JuMP.Model

    all_vars, n_vars, var_to_col = _build_var_layout(model)

    cols = _collect_socwr_cols(model, var_to_col, obj)
    shared_v2_cols =
        _collect_shared_v2_cols(model, var_to_col, part.zone_shared_bus_ids)

    if include_edges
        edge_r_cols, edge_t_cols =
            _collect_socwr_edge_rt_cols(model, var_to_col, part)
        edge_flow_cols, edge_pflow_coefs, edge_qflow_coefs,
        edge_pflow_const, edge_qflow_const =
            _collect_socwr_edge_flow_expansion(model, var_to_col, part)
    else
        edge_r_cols = Int[]
        edge_t_cols = Int[]
        edge_flow_cols = Vector{Vector{Int}}()
        edge_pflow_coefs = Vector{Vector{Float64}}()
        edge_qflow_coefs = Vector{Vector{Float64}}()
        edge_pflow_const = Float64[]
        edge_qflow_const = Float64[]
    end

    lb, ub = _collect_bounds(all_vars)

    # Constant Hessian diagonal (2·c2·sbase² on pg cols).
    base_hess_diag = zeros(Float64, n_vars)
    @inbounds for (k, col) in enumerate(cols.pg_cols)
        base_hess_diag[col] = cols.pg_c2_sbase2[k]
    end

    C_eq = _extract_equality_jacobian(model, var_to_col, n_vars)

    return AladinZoneDataSOCWR(
        model,
        all_vars,
        n_vars,
        var_to_col,
        cols.theta_cols,
        cols.v2_cols,
        cols.r_cols,
        cols.t_cols,
        cols.pg_cols,
        cols.qg_cols,
        cols.curt_cols,
        cols.pg_c2_sbase2,
        cols.pg_c1_sbase,
        shared_v2_cols,
        edge_r_cols,
        edge_t_cols,
        edge_flow_cols,
        edge_pflow_coefs,
        edge_qflow_coefs,
        edge_pflow_const,
        edge_qflow_const,
        lb,
        ub,
        _init_from_start_values(all_vars),  # x_prev
        zeros(Float64, n_vars),     # y
        zeros(Float64, n_vars),     # g
        copy(base_hess_diag),       # hess_diag (will be refreshed each iter)
        base_hess_diag,
        C_eq,
        zeros(Float64, n_vars),     # delta
        Int[],                      # active_bound_cols
    )
end

# ------------------------------------------------------------
# Per-iteration kernels
# ------------------------------------------------------------

function extract_solve_values!(zone::AladinZoneDataSOCWR, model::JuMP.Model)
    y = zone.y
    @inbounds for i in 1:zone.n_vars
        y[i] = JuMP.value(zone.all_vars[i])
    end
    return nothing
end

function compute_cost_gradient!(zone::AladinZoneDataSOCWR)
    g = zone.g
    fill!(g, 0.0)
    y = zone.y
    @inbounds for k in eachindex(zone.pg_cols)
        col = zone.pg_cols[k]
        g[col] = zone.pg_c2_sbase2[k] * y[col] + zone.pg_c1_sbase[k]
    end
    return nothing
end

function detect_active_bounds!(zone::AladinZoneDataSOCWR; tol::Float64 = 1e-6)
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

# ------------------------------------------------------------
# Proximal Hessian refresh
# ------------------------------------------------------------
#
# HFD 2016 Theorem 4.1 allows H_i to be any positive-definite
# approximation of the local Lagrangian Hessian. The cost term alone
# is PSD but only on pg columns, which leaves the remaining classes
# (θ, v², R, T, qg, curtailment) with zero curvature. We recover
# strict positive definiteness on the whole space by adding the
# curvature of the proximal term itself — exactly as the Phase-1
# subproblem optimises. This is the same construction used by
# ALADIN-DC and is the standard SQP PSD fix.

function _refresh_hess_diag!(
    zone::AladinZoneDataSOCWR,
    rho::Float64;
    sigma_theta::Float64,
    sigma_v2::Float64,
    sigma_r::Float64,
    sigma_t::Float64,
    sigma_pg::Float64,
    sigma_qg::Float64,
    sigma_other::Float64,
    delta::Float64 = 1e-6,
)
    copyto!(zone.hess_diag, zone.base_hess_diag)

    if rho > 0.0
        add_rho! = (col_set, σ) -> begin
            s = rho * σ
            s > 0.0 || return
            @inbounds for c in col_set
                zone.hess_diag[c] += s
            end
        end
        add_rho!(zone.theta_cols, sigma_theta)
        add_rho!(zone.v2_cols, sigma_v2)
        add_rho!(zone.r_cols, sigma_r)
        add_rho!(zone.t_cols, sigma_t)
        add_rho!(zone.pg_cols, sigma_pg)
        add_rho!(zone.qg_cols, sigma_qg)
        add_rho!(zone.curt_cols, sigma_other)
    end

    # Ensure a strictly positive floor on every column (a column with
    # zero σ would still leave a zero on the diagonal, making the QP
    # indefinite on that direction).
    @inbounds for c in 1:zone.n_vars
        if zone.hess_diag[c] < delta
            zone.hess_diag[c] = delta
        end
    end
    return nothing
end

# ------------------------------------------------------------
# Cone active-face normals (canonical HFD Alg. 1 for conic problems)
# ------------------------------------------------------------
# At an active (primal on the cone boundary) conic constraint
# `func(y) ∈ K`, HFD 2016 Algorithm 1 requires preserving the
# active-face hyperplane in the Phase-2 coordinator QP via a row
# `nᵀ · ∇func · Δy = 0`, where n is the outward normal of K at
# `func(y*)` and ∇func is the Jacobian of the affine constraint
# function w.r.t. the zone variables.
#
# We derive n from the squared-slack representation of each cone
# (so both n and its scaling are unambiguous):
#
#   Rotated SOC  K_rot = {(u, v, t₁, t₂, ...) : 2uv ≥ Σ tₖ², u,v ≥ 0}
#     slack g = Σ tₖ² - 2uv  ≤ 0
#     ∇g = (-2v, -2u, 2t₁, 2t₂, ...)
#
#   Standard SOC K_soc = {(t, x₁, x₂, ...) : t² ≥ Σ xₖ², t ≥ 0}
#     slack g = Σ xₖ² - t²   ≤ 0
#     ∇g = (-2t, 2x₁, 2x₂, ...)
#
# For affine functions in cones the row is constant-direction in
# variable space but its location in the cone depends on y*, so it
# must be refreshed every iteration (unlike the equality Jacobian
# C_eq which is structurally constant).

@inline function _cone_active_normal(
    ::MOI.RotatedSecondOrderCone,
    v::AbstractVector{Float64},
)
    k = length(v)
    n = zeros(Float64, k)
    n[1] = -2.0 * v[2]
    n[2] = -2.0 * v[1]
    @inbounds for i in 3:k
        n[i] = 2.0 * v[i]
    end
    return n
end

@inline function _cone_active_normal(
    ::MOI.SecondOrderCone,
    v::AbstractVector{Float64},
)
    k = length(v)
    n = zeros(Float64, k)
    n[1] = -2.0 * v[1]
    @inbounds for i in 2:k
        n[i] = 2.0 * v[i]
    end
    return n
end

@inline function _eval_aff_at_y(
    f::JuMP.AffExpr,
    y::Vector{Float64},
    var_to_col::Dict{JuMP.VariableRef, Int},
)
    acc = f.constant
    for (v, coef) in f.terms
        col = get(var_to_col, v, 0)
        col == 0 && continue
        acc += coef * y[col]
    end
    return acc
end

# ------------------------------------------------------------
# Active-constraint Jacobian
# ------------------------------------------------------------
# C_active = [ C_eq ; C_ineq_active ; C_cone_active ]
#
# SOCWR constraint categories:
#   * affine equalities (power balance, curtailment links) → C_eq
#     (precomputed at build time, structurally constant)
#   * affine scalar inequalities (branch angle bounds)     → walked
#     via scalar dual threshold
#   * vector-valued cone constraints (Jabr RSOC, thermal SOC) → walked
#     via cone dual norm threshold; at each active cone we compute the
#     outward normal in the cone's image space and pull it back to the
#     zone's variable layout via the affine function's Jacobian.
#
# This is the full HFD 2016 Algorithm 1 active-set treatment for conic
# problems. Inactive cones contribute zero to the Lagrangian Hessian
# (the affine functions have zero Hessian) and zero rows here; active
# cones contribute one row per constraint (the tangent hyperplane).

function build_active_jacobian(
    zone::AladinZoneDataSOCWR;
    dual_tol::Float64 = 1e-6,
    cone_tol::Float64 = 1e-6,
    dep_tol::Float64 = 1e-10,
)::SparseArrays.SparseMatrixCSC{Float64, Int}
    model = zone.model
    var_to_col = zone.var_to_col
    n_vars = zone.n_vars
    y = zone.y

    C_eq = zone.C_eq

    ineq_rows = Int[]
    ineq_cols = Int[]
    ineq_vals = Float64[]
    row_offset = 0

    has_duals = JuMP.has_duals(model)

    for (F, S) in JuMP.list_of_constraint_types(model)
        F <: JuMP.VariableRef && continue

        # --- Scalar affine inequality (angle bounds, etc.) ---
        if F <: JuMP.AffExpr
            S === MOI.EqualTo{Float64} && continue
            for ci in JuMP.all_constraints(model, F, S)
                has_duals || continue
                dval = try
                    JuMP.dual(ci)
                catch
                    0.0
                end
                (dval isa Real) || continue
                abs(dval) > dual_tol || continue

                cobj = JuMP.constraint_object(ci)
                row_offset += 1
                _push_affine_row_offset!(
                    ineq_rows, ineq_cols, ineq_vals,
                    row_offset, cobj.func, var_to_col,
                )
            end
            continue
        end

        # --- Vector cone constraints (RSOC, SOC) ---
        if F <: AbstractVector
            is_rsoc = S <: MOI.RotatedSecondOrderCone
            is_soc = S <: MOI.SecondOrderCone
            (is_rsoc || is_soc) || continue

            for ci in JuMP.all_constraints(model, F, S)
                has_duals || continue
                dval_vec = try
                    JuMP.dual(ci)
                catch
                    nothing
                end
                (dval_vec isa AbstractVector) || continue
                isempty(dval_vec) && continue
                LinearAlgebra.norm(dval_vec) > cone_tol || continue

                cobj = JuMP.constraint_object(ci)
                func = cobj.func
                cone_set = cobj.set

                m = length(func)
                v_star = Vector{Float64}(undef, m)
                ok = true
                @inbounds for k in 1:m
                    fk = func[k]
                    if fk isa JuMP.AffExpr
                        v_star[k] = _eval_aff_at_y(fk, y, var_to_col)
                    else
                        ok = false
                        break
                    end
                end
                ok || continue

                nrm = _cone_active_normal(cone_set, v_star)

                row_offset += 1
                @inbounds for k in 1:m
                    nk = nrm[k]
                    nk == 0.0 && continue
                    fk = func[k]
                    fk isa JuMP.AffExpr || continue
                    for (v, coef) in fk.terms
                        col = get(var_to_col, v, 0)
                        col == 0 && continue
                        push!(ineq_rows, row_offset)
                        push!(ineq_cols, col)
                        push!(ineq_vals, nk * coef)
                    end
                end
            end
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

# ------------------------------------------------------------
# Sensitivity refresh hook (called from aladin/core/local_solve.jl)
# ------------------------------------------------------------

@inline function refresh_lagrangian_hessian!(
    zone::AladinZoneDataSOCWR,
    ws::AbstractALADINWorkspace,
)
    _refresh_hess_diag!(
        zone,
        ws.rho[];
        sigma_theta = ws.sigma_theta,
        sigma_v2 = ws.sigma_v2,
        sigma_r = ws.sigma_r,
        sigma_t = ws.sigma_t,
        sigma_pg = ws.sigma_pg,
        sigma_qg = ws.sigma_qg,
        sigma_other = ws.sigma_other,
    )
    return nothing
end
