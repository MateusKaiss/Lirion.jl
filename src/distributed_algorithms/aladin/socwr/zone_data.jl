# ============================================================
# ALADIN-SOCWR — Zone-level data layout & builder
# ============================================================

# ============================================================
# Struct
# ============================================================

mutable struct AladinZoneDataSOCWR <: AbstractAladinZoneData
    local_model::JuMP.Model

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
    eq_jac::SparseArrays.SparseMatrixCSC{Float64, Int}

    # Phase-2 QP correction Δy_i.
    delta::Vector{Float64}

    # Active box bounds at current y.
    active_bound_cols::Vector{Int}

    # Cached quadratic part of Phase-1 proximal: (ρ/2)·Σ σ_k·v_k²
    # — see header note in `dc/zone_data.jl`. Mirrors `obj_static`
    # in ADMM. Refreshed only when ρ changes via `adapt_rho_mu!`.
    quad_static::JuMP.QuadExpr
end

# ============================================================
# Variable / cost / topology collectors
# ============================================================

# ------------------------------------------------------------
# Column classes (θ / v² / R / T / pg / qg / curt) + cost coefficients
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

    # `obj.generator_status[g] == 1` ⇔ JuMP sparse axis `gen.existing[g]`
    # (both derived from `gen_raw.gen_status`), so `pg_var[g, 1]` / `qg_var[g, 1]`
    # are guaranteed to exist whenever the status filter passes.
    if haskey(model, :p_generation)
        pg_var = model[:p_generation]
        for g in 1:obj.number_of_generators
            obj.generator_status[g] == 1 || continue
            v = pg_var[g, 1]
            col = get(var_to_col, v, 0)
            col == 0 && continue
            push!(pg_cols, col)

            # `obj.generator_costs[1]` is pre-allocated with `obj.number_of_generators`
            # entries in `build_objective_data`, so direct indexing is safe.
            # Matpower convention: descending degrees. For quadratic (len 3)
            # coefs = [c2, c1, c0].
            coefs = obj.generator_costs[1][g]
            @assert length(coefs) <= 3 "ALADIN: generator cost must be at most quadratic (got degree $(length(coefs)-1))"
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
            v = qg_var[g, 1]
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

# ------------------------------------------------------------
# Shared v² columns
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# Interzonal edge layout — R, T columns
# ------------------------------------------------------------

function _collect_socwr_edge_rt_cols(
    model::JuMP.Model,
    var_to_col::Dict{JuMP.VariableRef, Int},
    part::PartitionSOCWR;
    include_edges::Bool,
)
    include_edges || return Int[], Int[]

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
# Interzonal edge layout — affine flow expansion
# ------------------------------------------------------------

function _collect_socwr_edge_flow_expansion(
    model::JuMP.Model,
    var_to_col::Dict{JuMP.VariableRef, Int},
    part::PartitionSOCWR;
    include_edges::Bool,
)
    include_edges || return (
        Vector{Vector{Int}}(),
        Vector{Vector{Float64}}(),
        Vector{Vector{Float64}}(),
        Float64[],
        Float64[],
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

# ============================================================
# Builder
# ============================================================

function build_aladin_zone_data(
    part::PartitionSOCWR;
    include_edges::Bool = true,
)::AladinZoneDataSOCWR
    obj = part.objective_params
    model = part.local_model::JuMP.Model

    all_vars, n_vars, var_to_col = _build_var_layout(model)

    cols = _collect_socwr_cols(model, var_to_col, obj)
    shared_v2_cols =
        _collect_shared_v2_cols(model, var_to_col, part.zone_shared_bus_ids)

    edge_r_cols, edge_t_cols =
        _collect_socwr_edge_rt_cols(model, var_to_col, part; include_edges)
    edge_flow_cols, edge_pflow_coefs, edge_qflow_coefs,
    edge_pflow_const, edge_qflow_const =
        _collect_socwr_edge_flow_expansion(model, var_to_col, part; include_edges)

    lb, ub = _collect_bounds(all_vars)

    # Constant Hessian diagonal (2·c2·sbase² on pg cols).
    base_hess_diag = zeros(Float64, n_vars)
    @inbounds for (k, col) in enumerate(cols.pg_cols)
        base_hess_diag[col] = cols.pg_c2_sbase2[k]
    end

    eq_jac = _extract_equality_jacobian(model, var_to_col, n_vars)

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
        eq_jac,
        zeros(Float64, n_vars),     # delta
        Int[],                      # active_bound_cols
        zero(JuMP.QuadExpr),        # quad_static
    )
end
