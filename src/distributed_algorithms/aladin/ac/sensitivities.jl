mutable struct AladinZoneDataAC
    # JuMP model reference — used every iteration to walk constraints
    # and evaluate their current duals.
    model::JuMP.Model

    # Flat variable layout.
    all_vars::Vector{JuMP.VariableRef}
    n_vars::Int
    var_to_col::Dict{JuMP.VariableRef, Int}

    # Column classes for Σ_i sigma weights in the proximal term.
    theta_cols::Vector{Int}
    vm_cols::Vector{Int}
    pg_cols::Vector{Int}
    qg_cols::Vector{Int}
    curt_cols::Vector{Int}

    # Cost-support data on pg (quadratic cost in p_generation).
    pg_c2_sbase2::Vector{Float64}     # 2·c2·sbase²
    pg_c1_sbase::Vector{Float64}      # c1·sbase

    # Shared-bus column mappings (used by the coordinator QP when
    # assembling State/Hybrid consensus). Aligned with
    # `part.zone_shared_bus_ids`.
    shared_theta_cols::Vector{Int}
    shared_vm_cols::Vector{Int}

    # Per-interzonal-edge layout for Power/Hybrid couplings.
    # Each entry k corresponds to one entry in `part.zone_edges`:
    #   edge_from_theta_col[k] = col of θ[from_bus]
    #   edge_to_theta_col[k]   = col of θ[to_bus]
    #   edge_from_vm_col[k]    = col of v_m[from_bus]
    #   edge_to_vm_col[k]      = col of v_m[to_bus]
    # And the branch parameters needed to compute ∇p_flow / ∇q_flow
    # analytically at the current y:
    edge_from_theta_col::Vector{Int}
    edge_to_theta_col::Vector{Int}
    edge_from_vm_col::Vector{Int}
    edge_to_vm_col::Vector{Int}
    edge_g::Vector{Float64}
    edge_b::Vector{Float64}
    edge_b_sh::Vector{Float64}
    edge_inv_tap::Vector{Float64}
    edge_tap_angle::Vector{Float64}

    # Per-variable box bounds (NaN if absent).
    lb::Vector{Float64}
    ub::Vector{Float64}

    # Previous iterate (proximal centre) and Phase-1 solve result.
    x_prev::Vector{Float64}
    y::Vector{Float64}

    # g_i and H_i refreshed every iteration by `refresh_sensitivities_ac!`.
    g::Vector{Float64}
    hess::SparseArrays.SparseMatrixCSC{Float64, Int}

    # Active-constraint Jacobian refreshed every iteration.
    jac_cache::SparseArrays.SparseMatrixCSC{Float64, Int}

    # QP output.
    delta::Vector{Float64}

    # Active box bounds at current y.
    active_bound_cols::Vector{Int}
end

# ------------------------------------------------------------
# Column-class detection
# ------------------------------------------------------------

function _collect_ac_cols(
    model::JuMP.Model,
    var_to_col::Dict{JuMP.VariableRef, Int},
    obj::ObjectiveData,
)
    theta_cols = Int[]
    vm_cols = Int[]
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

    if haskey(model, :v_magnitude)
        for v in values(model[:v_magnitude].data)
            col = get(var_to_col, v, 0)
            col == 0 && continue
            push!(vm_cols, col)
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
        vm_cols = vm_cols,
        pg_cols = pg_cols,
        qg_cols = qg_cols,
        curt_cols = curt_cols,
        pg_c2_sbase2 = pg_c2_sbase2,
        pg_c1_sbase = pg_c1_sbase,
    )
end

function _collect_shared_ac_cols(
    model::JuMP.Model,
    var_to_col::Dict{JuMP.VariableRef, Int},
    shared_bus_ids::Vector{Int},
)
    theta_var = model[:theta]
    vm_var = model[:v_magnitude]
    theta_cols = Vector{Int}(undef, length(shared_bus_ids))
    vm_cols = Vector{Int}(undef, length(shared_bus_ids))
    @inbounds for i in eachindex(shared_bus_ids)
        b = shared_bus_ids[i]
        theta_cols[i] = var_to_col[theta_var[b, 1]]
        vm_cols[i] = var_to_col[vm_var[b, 1]]
    end
    return theta_cols, vm_cols
end

function _collect_ac_edge_layout(
    model::JuMP.Model,
    var_to_col::Dict{JuMP.VariableRef, Int},
    part::PartitionAC,
)
    zone_edges = part.zone_edges
    n = length(zone_edges)

    from_theta = Vector{Int}(undef, n)
    to_theta = Vector{Int}(undef, n)
    from_vm = Vector{Int}(undef, n)
    to_vm = Vector{Int}(undef, n)
    g_vec = Vector{Float64}(undef, n)
    b_vec = Vector{Float64}(undef, n)
    b_sh_vec = Vector{Float64}(undef, n)
    inv_tap_vec = Vector{Float64}(undef, n)
    tap_angle_vec = Vector{Float64}(undef, n)

    theta_var = model[:theta]
    vm_var = model[:v_magnitude]
    branch = part.branch
    branch_map = part.branch_idx_full_to_zone

    @inbounds for (k, e) in enumerate(zone_edges)
        l = branch_map[e]
        from_b = branch.from[l]
        to_b = branch.to[l]

        from_theta[k] = var_to_col[theta_var[from_b, 1]]
        to_theta[k] = var_to_col[theta_var[to_b, 1]]
        from_vm[k] = var_to_col[vm_var[from_b, 1]]
        to_vm[k] = var_to_col[vm_var[to_b, 1]]

        g_vec[k] = branch.g[l]
        b_vec[k] = branch.b[l]
        b_sh_vec[k] = branch.b_sh[l]
        inv_tap_vec[k] = branch.inv_tap[l]
        tap_angle_vec[k] = branch.tap_angle[l]
    end

    return (
        edge_from_theta_col = from_theta,
        edge_to_theta_col = to_theta,
        edge_from_vm_col = from_vm,
        edge_to_vm_col = to_vm,
        edge_g = g_vec,
        edge_b = b_vec,
        edge_b_sh = b_sh_vec,
        edge_inv_tap = inv_tap_vec,
        edge_tap_angle = tap_angle_vec,
    )
end

# ------------------------------------------------------------
# Build AladinZoneDataAC
# ------------------------------------------------------------

function build_aladin_zone_data_ac(
    part::PartitionAC;
    include_edges::Bool = false,
)::AladinZoneDataAC
    obj = part.objective_params
    model = part.model::JuMP.Model

    all_vars, n_vars, var_to_col = _build_var_layout(model)

    cols = _collect_ac_cols(model, var_to_col, obj)
    shared_theta_cols, shared_vm_cols =
        _collect_shared_ac_cols(model, var_to_col, part.zone_shared_bus_ids)

    edata =
        include_edges ?
        _collect_ac_edge_layout(model, var_to_col, part) :
        (
            edge_from_theta_col = Int[],
            edge_to_theta_col = Int[],
            edge_from_vm_col = Int[],
            edge_to_vm_col = Int[],
            edge_g = Float64[],
            edge_b = Float64[],
            edge_b_sh = Float64[],
            edge_inv_tap = Float64[],
            edge_tap_angle = Float64[],
        )

    lb, ub = _collect_bounds(all_vars)

    hess_initial = SparseArrays.spzeros(Float64, n_vars, n_vars)
    jac_initial = SparseArrays.spzeros(Float64, 0, n_vars)

    return AladinZoneDataAC(
        model,
        all_vars,
        n_vars,
        var_to_col,
        cols.theta_cols,
        cols.vm_cols,
        cols.pg_cols,
        cols.qg_cols,
        cols.curt_cols,
        cols.pg_c2_sbase2,
        cols.pg_c1_sbase,
        shared_theta_cols,
        shared_vm_cols,
        edata.edge_from_theta_col,
        edata.edge_to_theta_col,
        edata.edge_from_vm_col,
        edata.edge_to_vm_col,
        edata.edge_g,
        edata.edge_b,
        edata.edge_b_sh,
        edata.edge_inv_tap,
        edata.edge_tap_angle,
        lb,
        ub,
        zeros(Float64, n_vars),
        zeros(Float64, n_vars),
        zeros(Float64, n_vars),
        hess_initial,
        jac_initial,
        zeros(Float64, n_vars),
        Int[],
    )
end

# ------------------------------------------------------------
# Per-iteration kernels
# ------------------------------------------------------------

function extract_solve_values!(zone::AladinZoneDataAC, model::JuMP.Model)
    y = zone.y
    @inbounds for i in 1:zone.n_vars
        y[i] = JuMP.value(zone.all_vars[i])
    end
    return nothing
end

function compute_cost_gradient!(zone::AladinZoneDataAC)
    g = zone.g
    fill!(g, 0.0)
    y = zone.y
    @inbounds for k in eachindex(zone.pg_cols)
        col = zone.pg_cols[k]
        g[col] = zone.pg_c2_sbase2[k] * y[col] + zone.pg_c1_sbase[k]
    end
    return nothing
end

function detect_active_bounds!(zone::AladinZoneDataAC; tol::Float64 = 1e-6)
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
# Sensitivity extraction via an auxiliary MOI.Nonlinear.Model
# ------------------------------------------------------------
# We build a single aux NLP per zone per outer iteration. It contains:
#
#   * every active constraint of the zone's local JuMP model (active
#     in the strong-complementarity sense — equalities always,
#     inequalities when |κ_i| > tol),
#
#   * the generator quadratic cost as the NLP objective.
#
# This is exactly the `(f, g_i)` data HFD 2016 writes in Algorithm 1
# for the evaluation of g_i, H_i, C_i at y_i:
#
#       g_i = ∇f(y_i)
#       H_i = ∇²f(y_i) + Σ_{j ∈ active} κ_j ∇²g_j(y_i)
#       C_i = [∇g_j(y_i)ᵀ]_{j ∈ active}
#
# MOI's reverse-mode evaluator computes all three via autodiff:
#     MOI.eval_constraint_jacobian → C_i rows
#     MOI.eval_hessian_lagrangian  → H_i (lower triangle; we
#                                          symmetrise below)
#
# Inactive constraints are legitimately omitted because their
# multipliers are zero at a smooth KKT point and therefore contribute
# nothing to H_i. This is algorithmically equivalent to passing the
# full constraint list with zero multipliers, while avoiding the
# autodiff cost of unused derivatives.
# ============================================================

function refresh_sensitivities_ac!(
    zone::AladinZoneDataAC;
    rho::Float64 = 0.0,
    sigma_theta::Float64 = 0.0,
    sigma_vm::Float64 = 0.0,
    sigma_pg::Float64 = 0.0,
    sigma_qg::Float64 = 0.0,
    sigma_other::Float64 = 0.0,
    dual_tol::Float64 = 1e-3,
    dep_tol::Float64 = 1e-10,
)
    model = zone.model
    n_vars = zone.n_vars
    has_duals = JuMP.has_duals(model)

    nlp = MOI.Nonlinear.Model()
    rows = JuMP.ConstraintRef[]   # real constraints corresponding to aux NLP rows

    # ---- Add active constraints to the aux NLP ----
    for (F, S) in JuMP.list_of_constraint_types(model)
        F <: JuMP.VariableRef && continue

        for ci in JuMP.all_constraints(model, F, S)
            keep = if S === MOI.EqualTo{Float64}
                true
            elseif has_duals
                dval = try
                    JuMP.dual(ci)
                catch
                    0.0
                end
                abs(dval) > dual_tol
            else
                false
            end
            keep || continue

            cobj = JuMP.constraint_object(ci)
            MOI.Nonlinear.add_constraint(nlp, cobj.func, cobj.set)
            push!(rows, ci)
        end
    end

    # ---- Set objective = generation cost ----
    obj_expr = _build_ac_cost_expr(zone)
    if obj_expr !== nothing
        MOI.Nonlinear.set_objective(nlp, obj_expr)
    end

    # ---- Build the evaluator once, initialise for Jac + Hess ----
    ordered_vars = [JuMP.index(v) for v in zone.all_vars]
    evaluator = MOI.Nonlinear.Evaluator(
        nlp,
        MOI.Nonlinear.SparseReverseMode(),
        ordered_vars,
    )
    MOI.initialize(evaluator, [:Jac, :Hess])

    y_primal = zone.y

    # ---- Jacobian ----
    J_struct = MOI.jacobian_structure(evaluator)
    V_J = zeros(length(J_struct))
    MOI.eval_constraint_jacobian(evaluator, V_J, y_primal)
    n_rows = length(rows)
    if n_rows == 0
        zone.jac_cache = SparseArrays.spzeros(Float64, 0, n_vars)
    else
        I_J = [t[1] for t in J_struct]
        C_J = [t[2] for t in J_struct]
        C_full = SparseArrays.sparse(I_J, C_J, V_J, n_rows, n_vars)
        zone.jac_cache = _delete_dependent_rows(C_full; tol = dep_tol)
    end

    # ---- Hessian of Lagrangian ----
    H_struct = MOI.hessian_lagrangian_structure(evaluator)
    V_H = zeros(length(H_struct))
    y_duals = Float64[-JuMP.dual(ci) for ci in rows]
    if !isempty(H_struct)
        MOI.eval_hessian_lagrangian(evaluator, V_H, y_primal, 1.0, y_duals)
    end

    I_H = [t[1] for t in H_struct]
    C_H = [t[2] for t in H_struct]
    H_lower = SparseArrays.sparse(I_H, C_H, V_H, n_vars, n_vars)
    # MOI returns only the lower triangle — symmetrise explicitly.
    H_full =
        H_lower + SparseArrays.transpose(H_lower) -
        SparseArrays.spdiagm(0 => LinearAlgebra.diag(H_lower))

    # Augmented-Lagrangian curvature (ρ·diag(Σ)). HFD 2016 Theorem 4.1
    # specifies H_i ≈ ∇²L_i with the explicit "≈" allowing any positive
    # definite approximation of the local Lagrangian Hessian. For
    # non-convex AC OPF the raw Hessian ∇²f + Σ κ_j ∇²g_j computed
    # above is generally indefinite (active power-balance multipliers
    # can drive whole eigen-blocks negative), so we add the curvature
    # of the proximal term that the Phase-1 subproblem itself
    # optimises. The resulting block is exactly the Hessian of the
    # PENALISED local Lagrangian, which is positive-definite by
    # construction whenever ρ > 0 and Σ ≻ 0. This is the standard SQP
    # construction for non-convex NLPs and recovers the canonical
    # Newton step in the limit ρ → 0 on convex problems.
    if rho > 0.0
        add! = (col_set, σ) -> begin
            s = rho * σ
            s > 0.0 || return
            for c in col_set
                H_full[c, c] += s
            end
        end
        add!(zone.theta_cols, sigma_theta)
        add!(zone.vm_cols, sigma_vm)
        add!(zone.pg_cols, sigma_pg)
        add!(zone.qg_cols, sigma_qg)
        add!(zone.curt_cols, sigma_other)
    end

    # PSD projection by minimum-eigenvalue shift.
    #
    # For non-convex AC OPF the Lagrangian Hessian ∇²f + Σ κ_j ∇²g_j
    # can be indefinite even after the proximal regularisation. The
    # standard SQP remedy (used by Ipopt's interior-point algorithm
    # and by the original ALADIN convergence proofs in HFD 2016
    # Section 4) is to shift the diagonal of H_i so that
    #     H_i ≽ δ · I   for some small δ > 0.
    # We compute the smallest eigenvalue of the dense H_i and add
    #     τ = max(0, δ − λ_min(H_i))
    # to every diagonal entry. This is the minimum perturbation that
    # makes H_i positive definite, preserves the eigen-directions on
    # which the original Hessian was already curved, and degenerates
    # to the identity transformation on convex problems.
    if n_vars > 0
        H_dense = Matrix(H_full)
        H_dense = 0.5 * (H_dense + H_dense')   # numerical symmetry
        λ_min = LinearAlgebra.eigvals(LinearAlgebra.Symmetric(H_dense))[1]
        δ = 1e-6
        if λ_min < δ
            shift = δ - λ_min
            @inbounds for c in 1:n_vars
                H_full[c, c] += shift
            end
        end
    end

    zone.hess = H_full
    return nothing
end

# Build the objective `f(x) = Σ (c2·pg² + c1·pg)` as a Julia Expr tree
# referring to the real model's MOI.VariableIndex literals.
function _build_ac_cost_expr(zone::AladinZoneDataAC)
    parts = Any[]
    @inbounds for k in eachindex(zone.pg_cols)
        col = zone.pg_cols[k]
        vi = JuMP.index(zone.all_vars[col])
        c2 = zone.pg_c2_sbase2[k] / 2.0
        c1 = zone.pg_c1_sbase[k]
        if c2 != 0.0
            push!(parts, :($c2 * $vi * $vi))
        end
        if c1 != 0.0
            push!(parts, :($c1 * $vi))
        end
    end
    isempty(parts) && return nothing
    length(parts) == 1 && return parts[1]
    return Expr(:call, :+, parts...)
end

# Hot-path return used by coordinator_qp's `_add_active_constraints!`.
@inline build_active_jacobian(
    zone::AladinZoneDataAC;
    kwargs...,
)::SparseArrays.SparseMatrixCSC{Float64, Int} = zone.jac_cache
