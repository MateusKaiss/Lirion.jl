# ============================================================
# ALADIN-AC — Per-iteration sensitivities (MOI.Nonlinear.Evaluator)
# ============================================================

# ============================================================
# AC cost expression
# ============================================================

# Shared by AC and the SOCWRNL exact-Hessian path (explicit Union keeps DC zones out).
function _build_ac_cost_expr(zone::Union{AladinZoneDataAC, AladinZoneDataSOCWR})
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

# ============================================================
# Active-constraint NLP + evaluator
# ============================================================

@inline function _ac_keep_constraint(
    ci,
    ::Type{S},
    has_duals::Bool,
    dual_tol::Float64,
) where {S}
    S === MOI.EqualTo{Float64} && return true
    has_duals || return false
    dval = try
        JuMP.dual(ci)
    catch
        0.0
    end
    return abs(dval) > dual_tol
end

# Modified gradient (HFD Algorithm 2, step 3): gᵢ = ∇fᵢ(yᵢ) + (Cᵢ* − Cᵢ)ᵀκᵢ.
const _MODIFIED_GRADIENT = Ref(false)

# Weak-row floor: rows with |dual| in (floor, dual_tol] enter the evaluator (not C).
const _MODGRAD_DUAL_FLOOR = 1e-9

# Shared with the SOCWRNL exact-Hessian path (see `_build_ac_cost_expr`).
function _build_ac_evaluator(
    zone::Union{AladinZoneDataAC, AladinZoneDataSOCWR},
    dual_tol::Float64,
)
    model = zone.local_model
    has_duals = JuMP.has_duals(model)
    modgrad = _MODIFIED_GRADIENT[] && has_duals

    nlp = MOI.Nonlinear.Model()
    rows = JuMP.ConstraintRef[]
    strong = Bool[]   # true -> C̃ row (CΔ=0); false -> gradient-correction row

    for (F, S) in JuMP.list_of_constraint_types(model)
        F <: JuMP.VariableRef && continue
        for ci in JuMP.all_constraints(model, F, S)
            is_strong = _ac_keep_constraint(ci, S, has_duals, dual_tol)
            keep = is_strong
            if !keep && modgrad
                dval = try
                    JuMP.dual(ci)
                catch
                    0.0
                end
                keep = abs(dval) > _MODGRAD_DUAL_FLOOR
            end
            keep || continue
            cobj = JuMP.constraint_object(ci)
            MOI.Nonlinear.add_constraint(nlp, cobj.func, cobj.set)
            push!(rows, ci)
            push!(strong, is_strong)
        end
    end

    obj_expr = _build_ac_cost_expr(zone)
    obj_expr === nothing || MOI.Nonlinear.set_objective(nlp, obj_expr)

    ordered_vars = [JuMP.index(v) for v in zone.all_vars]
    evaluator = MOI.Nonlinear.Evaluator(
        nlp,
        MOI.Nonlinear.SparseReverseMode(),
        ordered_vars,
    )
    MOI.initialize(evaluator, [:Jac, :Hess])
    return evaluator, rows, strong
end

# ============================================================
# Active-constraint Jacobian
# ============================================================

# Returns (C̃ from strong rows, gradient correction Σ_weak κⱼ∇hⱼ from weak rows).
function _extract_active_jacobian_ac(
    evaluator,
    rows::Vector{JuMP.ConstraintRef},
    strong::Vector{Bool},
    y::Vector{Float64},
    n_vars::Int,
    dep_tol::Float64,
)
    n_rows = length(rows)
    g_corr = zeros(n_vars)
    n_rows == 0 &&
        return SparseArrays.spzeros(Float64, 0, n_vars), g_corr

    J_struct = MOI.jacobian_structure(evaluator)
    V_J = zeros(length(J_struct))
    MOI.eval_constraint_jacobian(evaluator, V_J, y)

    I_J = [t[1] for t in J_struct]
    C_J = [t[2] for t in J_struct]
    C_all = SparseArrays.sparse(I_J, C_J, V_J, n_rows, n_vars)

    if all(strong)
        return _delete_dependent_rows(C_all; tol = dep_tol), g_corr
    end

    sidx = findall(strong)
    widx = findall(!, strong)
    # weak rows: g_corr = Σ (-dualⱼ)·∇hⱼ (internal Lagrangian sign)
    Ct = SparseArrays.sparse(C_all')          # n_vars × n_rows, column = row grad
    rv = SparseArrays.rowvals(Ct)
    nz = SparseArrays.nonzeros(Ct)
    for r in widx
        kappa = -(try; JuMP.dual(rows[r]); catch; 0.0; end)
        kappa == 0.0 && continue
        @inbounds for k in SparseArrays.nzrange(Ct, r)
            g_corr[rv[k]] += kappa * nz[k]
        end
    end
    C_strong = C_all[sidx, :]
    return _delete_dependent_rows(C_strong; tol = dep_tol), g_corr
end

# ============================================================
# Lagrangian Hessian (symmetrised from lower triangle)
# ============================================================

function _extract_lagrangian_hessian_ac(
    evaluator,
    rows::Vector{JuMP.ConstraintRef},
    y::Vector{Float64},
    n_vars::Int,
)::SparseArrays.SparseMatrixCSC{Float64, Int}
    H_struct = MOI.hessian_lagrangian_structure(evaluator)
    # the symmetrisation below assumes lower-triangular structure entries
    @assert all(t -> t[1] >= t[2], H_struct) "MOI Hessian structure has upper-triangular entries; the lower-triangle symmetrisation in _extract_lagrangian_hessian_ac would double-count them"
    V_H = zeros(length(H_struct))
    y_duals = Float64[-JuMP.dual(ci) for ci in rows]
    if !isempty(H_struct)
        MOI.eval_hessian_lagrangian(evaluator, V_H, y, 1.0, y_duals)
    end

    I_H = [t[1] for t in H_struct]
    C_H = [t[2] for t in H_struct]
    H_lower = SparseArrays.sparse(I_H, C_H, V_H, n_vars, n_vars)

    return H_lower + SparseArrays.transpose(H_lower) -
           SparseArrays.spdiagm(0 => LinearAlgebra.diag(H_lower))
end

# ============================================================
# Diagonal regularisation: proximal ρ·σ + Tikhonov
# ============================================================

function _add_prox_curvature_ac!(
    H::SparseArrays.SparseMatrixCSC{Float64, Int},
    zone::AladinZoneDataAC,
    rho::Float64;
    sigma_theta::Float64,
    sigma_vm::Float64,
    sigma_pg::Float64,
    sigma_qg::Float64,
    sigma_other::Float64,
)
    _accumulate_prox_diagonal!(
        rho,
        (
            (zone.theta_cols, sigma_theta),
            (zone.vm_cols, sigma_vm),
            (zone.pg_cols, sigma_pg),
            (zone.qg_cols, sigma_qg),
            (zone.curt_cols, sigma_other),
        ),
    ) do c, s
        return H[c, c] += s
    end
    return nothing
end

@inline function _add_tikhonov_ac!(
    H::SparseArrays.SparseMatrixCSC{Float64, Int},
    n_vars::Int;
    delta::Float64 = 1e-6,
)
    @inbounds for c in 1:n_vars
        H[c, c] += delta
    end
    @assert all(isfinite, H.nzval) "AC Hessian has non-finite entries"
    return nothing
end

# ============================================================
# Orchestrator
# ============================================================

# Canonical-QP toggle: B_i ≈ ∇²L_i only (HFD Alg. 2), without the proximal ρΣ term.
const _AC_QP_NO_PROX = Ref(false)

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
    y = zone.y
    n_vars = zone.n_vars

    evaluator, rows, strong = _build_ac_evaluator(zone, dual_tol)

    zone.active_jac, g_corr =
        _extract_active_jacobian_ac(evaluator, rows, strong, y, n_vars, dep_tol)
    if _MODIFIED_GRADIENT[]
        g = zone.g
        @inbounds for c in 1:n_vars
            g[c] += g_corr[c]
        end
        if (try; JuMP.has_duals(zone.local_model); catch; false; end)
            act = zone.active_bound_cols
            for c in 1:n_vars
                c in act && continue
                vr = zone.all_vars[c]
                du = 0.0
                try
                    JuMP.has_lower_bound(vr) &&
                        (du += JuMP.dual(JuMP.LowerBoundRef(vr)))
                    JuMP.has_upper_bound(vr) &&
                        (du += JuMP.dual(JuMP.UpperBoundRef(vr)))
                catch
                end
                # internal Lagrangian sign, matching y_duals = -dual
                du != 0.0 && (g[c] += -du)
            end
        end
    end
    H = _extract_lagrangian_hessian_ac(evaluator, rows, y, n_vars)

    if !_AC_QP_NO_PROX[]
        _add_prox_curvature_ac!(
            H, zone, rho;
            sigma_theta = sigma_theta,
            sigma_vm = sigma_vm,
            sigma_pg = sigma_pg,
            sigma_qg = sigma_qg,
            sigma_other = sigma_other,
        )
    end
    n_vars > 0 && _add_tikhonov_ac!(H, n_vars)
    zone.hess = H
    return nothing
end

# ============================================================
# :pflow — Gauss-Newton flow curvature + PSD safeguard (AC only)
# ============================================================

# Accumulate ρ·σ·(g)(g)ᵀ over the 4 coupling columns of one edge into (I,J,V).
@inline function _push_flow_outer!(
    I::Vector{Int}, J::Vector{Int}, V::Vector{Float64},
    cols::NTuple{4, Int}, grad::NTuple{4, Float64}, coef::Float64,
)
    coef == 0.0 && return nothing
    @inbounds for a in 1:4
        ga = grad[a]
        ga == 0.0 && continue
        ca = cols[a]
        for b in 1:4
            v = coef * ga * grad[b]
            v == 0.0 && continue
            push!(I, ca)
            push!(J, cols[b])
            push!(V, v)
        end
    end
    return nothing
end

function _add_flow_gn_curvature_ac(
    H::SparseArrays.SparseMatrixCSC{Float64, Int},
    zone::AladinZoneDataAC,
    ws,
    zone_pos::Int,
    rho::Float64,
    sigma_pflow::Float64,
    sigma_qflow::Float64,
)::SparseArrays.SparseMatrixCSC{Float64, Int}
    n = zone.n_vars
    I = Int[]
    J = Int[]
    V = Float64[]
    rp = rho * sigma_pflow
    rq = rho * sigma_qflow

    @inbounds for pair in ws.interzonal_line_pairs
        le = pair.zone_i == zone_pos ? pair.line_local_i :
             (pair.zone_j == zone_pos ? pair.line_local_j : 0)
        le == 0 && continue
        lin = _ac_flow_linearisation(zone, le)
        cols = (
            zone.edge_from_theta_col[le], zone.edge_to_theta_col[le],
            zone.edge_from_vm_col[le], zone.edge_to_vm_col[le],
        )
        _push_flow_outer!(I, J, V, cols,
            (lin.dpf_dtheta_f, lin.dpf_dtheta_t, lin.dpf_dvf, lin.dpf_dvt), rp)
        _push_flow_outer!(I, J, V, cols,
            (lin.dqf_dtheta_f, lin.dqf_dtheta_t, lin.dqf_dvf, lin.dqf_dvt), rq)
    end

    @inbounds for pair in ws.dummy_line_pairs
        if pair.zone_i == zone_pos
            # half-A: the flow proximal is on the TO-side flow; use the same linearisation
            lin = _dummy_ac_flow_linearisation_to(
                zone.y, pair.from_theta_col_i, pair.to_theta_col_i,
                pair.from_vm_col_i, pair.to_vm_col_i,
                pair.g_i, pair.b_i, pair.b_sh_i, pair.inv_tap_i, pair.tap_angle_i)
            cols = (pair.from_theta_col_i, pair.to_theta_col_i,
                pair.from_vm_col_i, pair.to_vm_col_i)
        elseif pair.zone_j == zone_pos
            lin = _dummy_ac_flow_linearisation(
                zone.y, pair.from_theta_col_j, pair.to_theta_col_j,
                pair.from_vm_col_j, pair.to_vm_col_j,
                pair.g_j, pair.b_j, pair.b_sh_j, pair.inv_tap_j, pair.tap_angle_j)
            cols = (pair.from_theta_col_j, pair.to_theta_col_j,
                pair.from_vm_col_j, pair.to_vm_col_j)
        else
            continue
        end
        _push_flow_outer!(I, J, V, cols,
            (lin.dpf_dtheta_f, lin.dpf_dtheta_t, lin.dpf_dvf, lin.dpf_dvt), rp)
        _push_flow_outer!(I, J, V, cols,
            (lin.dqf_dtheta_f, lin.dqf_dtheta_t, lin.dqf_dvf, lin.dqf_dvt), rq)
    end

    isempty(V) && return H
    return H + SparseArrays.sparse(I, J, V, n, n)
end

# ------------------------------------------------------------
# PSD-safeguard mode (LIRION_AC_PSD_MODE ∈ {shift, flip, clip})
# ------------------------------------------------------------
function _resolve_ac_psd_mode()
    raw = get(ENV, "LIRION_AC_PSD_MODE", "shift")
    s = lowercase(strip(raw))
    s == "shift" && return :shift
    s == "flip" && return :flip
    s == "clip" && return :clip
    error("LIRION_AC_PSD_MODE=$(repr(raw)) is invalid; use one of shift|flip|clip")
end
const _PSD_MODE = Ref{Symbol}(_resolve_ac_psd_mode())

function _ensure_psd_ac!(
    H::SparseArrays.SparseMatrixCSC{Float64, Int};
    margin::Float64 = 1e-4,
    max_tries::Int = 64,
)
    n = size(H, 1)
    n == 0 && return H
    prox = _QP_HESS_REG[]
    if prox != 0.0
        @inbounds for c in 1:n
            H[c, c] += prox
        end
    end
    # Cheap path: already PD -> no shift (avoid over-regularizing the common case).
    if (try; LinearAlgebra.cholesky(LinearAlgebra.Symmetric(H)); true; catch; false; end)
        return H
    end
    # Spectral safeguard (:flip / :clip), dense-eigen branch only; returns a new matrix.
    mode = _PSD_MODE[]
    if mode !== :shift && n <= 800
        F = LinearAlgebra.eigen(LinearAlgebra.Symmetric(Matrix(H)))
        lam = F.values
        eps_floor = margin * max(1.0, maximum(abs, lam))   # same 1e-4·max(1,|λ|) scale as the shift margin
        lam_new = similar(lam)
        if mode === :flip
            @inbounds for i in eachindex(lam)
                a = abs(lam[i])
                lam_new[i] = a < eps_floor ? eps_floor : a
            end
        else  # :clip
            @inbounds for i in eachindex(lam)
                lam_new[i] = lam[i] < eps_floor ? eps_floor : lam[i]
            end
        end
        Hd = F.vectors * LinearAlgebra.Diagonal(lam_new) * transpose(F.vectors)
        Hd = (Hd .+ transpose(Hd)) .* 0.5           # kill round-off asymmetry
        corr = 0.0
        @inbounds for i in eachindex(lam)
            corr = max(corr, lam_new[i] - lam[i])
        end
        if corr > 0.0
            lock(_PSD_DELTA_LOCK) do
                _PSD_DELTA_LAST[] = max(_PSD_DELTA_LAST[], corr)
            end
        end
        return SparseArrays.sparse(Hd)
    end
    delta = 0.0
    if n <= 800
        lam = LinearAlgebra.eigmin(LinearAlgebra.Symmetric(Matrix(H)))
        delta = lam < 0.0 ? (-lam + margin * max(1.0, abs(lam))) : (margin * max(1.0, lam))
    else
        dmax = 0.0
        @inbounds for c in 1:n
            dmax = max(dmax, abs(H[c, c]))
        end
        d = max(1e-6, 1e-8 * dmax)
        @inbounds for _ in 0:max_tries
            ok = try
                LinearAlgebra.cholesky(LinearAlgebra.Symmetric(H + d * LinearAlgebra.I))
                true
            catch
                false
            end
            if ok
                delta = d
                break
            end
            d *= 2.0
        end
        delta == 0.0 && (delta = d)
    end
    @inbounds for c in 1:n
        H[c, c] += delta
    end
    if delta > 0.0
        lock(_PSD_DELTA_LOCK) do
            _PSD_DELTA_LAST[] = max(_PSD_DELTA_LAST[], delta)
        end
    end
    return H
end

# ============================================================
# Sensitivity refresh hook
# ============================================================

const _ACTIVE_DEP_TOL = Ref(1e-10)

const _ACTIVE_DUAL_TOL = Ref(1e-3)

function refresh_lagrangian_hessian!(
    zone::AladinZoneDataAC,
    ws::AbstractALADINWorkspace,
)
    pflow = hasproperty(ws, :proximal_mode) && ws.proximal_mode === :pflow
    sigma_theta_eff = pflow ? ws.sigma_other : ws.sigma_theta
    sigma_vm_eff = pflow ? ws.sigma_other : ws.sigma_vm

    refresh_sensitivities_ac!(
        zone;
        rho = ws.rho[],
        sigma_theta = sigma_theta_eff,
        sigma_vm = sigma_vm_eff,
        sigma_pg = ws.sigma_pg,
        sigma_qg = ws.sigma_qg,
        sigma_other = ws.sigma_other,
        dual_tol = _ACTIVE_DUAL_TOL[],
        dep_tol = _ACTIVE_DEP_TOL[],
    )

    if pflow
        zone_pos = findfirst(z -> z === zone, ws.zones)::Int
        zone.hess = _add_flow_gn_curvature_ac(
            zone.hess, zone, ws, zone_pos, ws.rho[],
            ws.sigma_pflow, ws.sigma_qflow,
        )
    end
    zone.hess = _ensure_psd_ac!(zone.hess)   # :flip/:clip return a new matrix
    return nothing
end
