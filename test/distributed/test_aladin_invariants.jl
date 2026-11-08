# ============================================================
# ==========  ALADIN mathematical invariants  ================
# ============================================================
#
# Solver-independent checks of the ALADIN building blocks against
# Houska-Frasch-Diehl (2016), Alg. 2 / Alg. 3 — see the header below.
# Everything lives in a module so the helper names cannot clash with the
# other test files.

module AladinInvariants

# ============================================================
# ALADIN mathematical invariants (independent of Lirion's own QP/KKT code)
# Reference: Houska, Frasch, Diehl (2016), SIAM J. Optim. 26(2) — Alg. 2 / Alg. 3
# ============================================================
#
# For every (model × coupling × sectioning) cell, after K warm-up iterations:
#  T_H  residual-definition consistency: rows(y) vs ws.r_primal, ‖·‖₁ vs Φ's L1 term,
#       _dual_slope(ws, r) == ‖r‖²
#  T_C  LOCAL KKT (3.2):  ∇f + Aᵢᵀλ + ρΣ(y−x) − Σ dualⱼ∇hⱼ − bound duals = 0
#       with A := FD Jacobian of the consensus rows (independent of the hand
#       derivatives), f-gradient := zone.g. Validates the λ-term sign/orientation,
#       the prox assembly, the cost gradient (factor 2·c2·sbase²) and the
#       multiplier convention κ = −JuMP.dual used by Hessian / modified gradient.
#  T_B  COORDINATION QP (3.3) re-solved in dense linear algebra from
#       (Hᵢ, gᵢ, Cᵢ, A_FD, λ, μ) via the null-space method; Δy and λ_QP are
#       compared with Lirion's Clarabel output (zone.delta, state.nu).
#       Also the slack identity λ_QP = λ + μ s (HFD footnote 4).
#  T_D  AC Lagrangian Hessian vs central finite differences of ∇ₓL.
#  T_G  Danskin slope of V_ρ along d = λ_QP − λ vs FD of V_ρ (step (c) of Alg. 3).
#  T_N  SOCWR cone normal vs FD; cone distance sanity.

using Lirion, JuMP, LinearAlgebra, SparseArrays
const L = Lirion
const MOI = L.MOI

const RESULTS = Vector{NamedTuple}()
record!(cell, test, val, tol) = (push!(RESULTS, (; cell, test, val, tol, ok = val <= tol)); nothing)

# ------------------------------------------------------------
# workspace construction + manual stepping
# ------------------------------------------------------------
function build(file, model, coupling, sectioning; kw...)
    alg = L.aladin(; coupling, sectioning, max_iter = 50, kw...)
    inputs = L._build_inputs(String(file), alg; model, optimizer = nothing, outputs = false)
    cfg = inputs.model_type
    opt = L.resolve_optimizer(cfg, inputs.optimizer)
    tag = L.ALADIN{L._variant_type(Val(cfg.coupling)), cfg.base_model}()
    ws = L.build_workspace(inputs.system, cfg, opt, tag)
    return ws, L.build_params(cfg)
end

function half_iter!(ws, params, it)
    ws.iter = it
    L.local_solve!(ws)
    L.solve_coordinator_qp!(ws)
    L.compute_residuals!(ws)
    L.compute_tolerances!(ws, params.eps_abs, params.eps_rel)
end
function finish_iter!(ws, params, a3)
    ws.alpha1 = 1.0; ws.alpha3 = a3
    L._line_search_step!(ws)
    L.primal_update!(ws); L.dual_update!(ws); L.adapt_rho_mu!(ws, params)
end

# ------------------------------------------------------------
# consensus rows r(y) = Σ Aᵢ yᵢ − b  (Lirion's own row functions, arbitrary point)
# ------------------------------------------------------------
acc_y(z, c) = z.y[c]
rows_nt(ws::L._ACPhase2Workspace) = L._dual_residual(ws)
rows_nt(ws::L._DCPhase2Workspace) = L._dc_residual_rows(ws, acc_y)
rows_nt(ws::L._SOCWRPhase2Workspace) = L._socwr_consensus_rows(ws, acc_y)
flat(nt) = reduce(vcat, [Vector{Float64}(nt[k]) for k in keys(nt)]; init = Float64[])

offsets(ws) = cumsum([0; [z.n_vars for z in ws.zones]])
function getY(ws); reduce(vcat, [copy(z.y) for z in ws.zones]); end
function setY!(ws, Y); o = offsets(ws); for (i, z) in enumerate(ws.zones); copyto!(z.y, @view Y[o[i]+1:o[i+1]]); end; end
function rows_at(ws, Y)
    Y0 = getY(ws); setY!(ws, Y)
    nt = rows_nt(ws); @assert keys(nt) == keys(ws.vars) "row groups not in ws.vars order"
    r = flat(nt); setY!(ws, Y0); return r
end
function fd_jac(ws, Y; h = 1e-6)
    m = length(rows_at(ws, Y)); n = length(Y); A = zeros(m, n)
    for j in 1:n
        Yp = copy(Y); Ym = copy(Y); Yp[j] += h; Ym[j] -= h
        A[:, j] = (rows_at(ws, Yp) .- rows_at(ws, Ym)) ./ (2h)
    end
    return A
end
lam_vec(ws) = flat(map(s -> s.lambda, ws.vars))
nu_vec(ws) = flat(map(s -> s.nu, ws.vars))

# ------------------------------------------------------------
# prox gradient ρΣ(y−x) from quad_static (+ AC :pflow hook)
# ------------------------------------------------------------
function prox_grad(ws, zi)
    z = ws.zones[zi]; g = zeros(z.n_vars); v2c = z.var_to_col
    for (p, a) in z.quad_static.terms
        i = v2c[p.a]; j = v2c[p.b]
        if i == j; g[i] += 2a * (z.y[i] - z.x_prev[i])
        else; g[i] += a * (z.y[j] - z.x_prev[j]); g[j] += a * (z.y[i] - z.x_prev[i]); end
    end
    L._add_extra_prox_grad!(g, ws, zi, z)
    return g
end

# ------------------------------------------------------------
# −Σ dualⱼ ∇funcⱼ over ALL local constraints (scalar, vector cones, bounds)
# JuMP Min convention: L = f − Σ dualⱼ·(funcⱼ − setⱼ)
# ------------------------------------------------------------
# which = :all    → every constraint/bound (true KKT of (3.2))
# which = :strong → only what the coordinator enforces as CᵢΔy = 0 (equalities,
#                   |dual| > dual_tol rows, detected-active bounds): with the HFD
#                   modified gradient ON, g_mod + Aᵀλ + prox + (this) must vanish.
function constraint_term(z; which = :all, dual_tol = L._ACTIVE_DUAL_TOL[])
    model = z.local_model; n = z.n_vars; out = zeros(n)
    nlp = MOI.Nonlinear.Model(); rows = JuMP.ConstraintRef[]
    for (F, S) in JuMP.list_of_constraint_types(model)
        F <: JuMP.VariableRef && continue
        if F <: AbstractVector
            for ci in JuMP.all_constraints(model, F, S)
                d = JuMP.dual(ci); f = JuMP.constraint_object(ci).func
                for k in eachindex(f)
                    fk = f[k]; fk isa JuMP.AffExpr || error("non-affine cone arg")
                    for (v, c) in fk.terms; out[z.var_to_col[v]] -= d[k] * c; end
                end
            end
        else
            for ci in JuMP.all_constraints(model, F, S)
                if which === :strong && !(S === MOI.EqualTo{Float64}) && abs(JuMP.dual(ci)) <= dual_tol
                    continue
                end
                co = JuMP.constraint_object(ci)
                MOI.Nonlinear.add_constraint(nlp, co.func, co.set); push!(rows, ci)
            end
        end
    end
    if !isempty(rows)
        ev = MOI.Nonlinear.Evaluator(nlp, MOI.Nonlinear.SparseReverseMode(), [JuMP.index(v) for v in z.all_vars])
        MOI.initialize(ev, [:Jac])
        w = Float64[-JuMP.dual(ci) for ci in rows]
        tmp = zeros(n); MOI.eval_constraint_jacobian_transpose_product(ev, tmp, z.y, w); out .+= tmp
    end
    for (c, v) in enumerate(z.all_vars)
        which === :strong && !(c in z.active_bound_cols) && continue
        JuMP.has_lower_bound(v) && (out[c] -= JuMP.dual(JuMP.LowerBoundRef(v)))
        JuMP.has_upper_bound(v) && (out[c] -= JuMP.dual(JuMP.UpperBoundRef(v)))
        JuMP.is_fixed(v) && (out[c] -= JuMP.dual(JuMP.FixRef(v)))
    end
    return out
end

# independent ∇f from the JuMP cost expression (not from zone.pg_c2_sbase2)
function grad_obj_base(part, z)
    e = part.obj_base; g = zeros(z.n_vars)
    if e isa JuMP.GenericQuadExpr
        for (v, c) in e.aff.terms; g[z.var_to_col[v]] += c; end
        for (p, c) in e.terms
            i = z.var_to_col[p.a]; j = z.var_to_col[p.b]
            if i == j; g[i] += 2c * z.y[i] else g[i] += c * z.y[j]; g[j] += c * z.y[i] end
        end
    elseif e isa JuMP.GenericAffExpr
        for (v, c) in e.terms; g[z.var_to_col[v]] += c; end
    elseif e isa Real
        nothing                                   # zone without costed generators
    else
        error("obj_base of type $(typeof(e)) not handled")
    end
    return g
end

# ------------------------------------------------------------
# Hᵢ and Cᵢ exactly as the coordinator consumes them
# ------------------------------------------------------------
Hmat(z::L.AladinZoneDataAC) = Matrix(z.hess)
Hmat(z) = Matrix(Diagonal(z.hess_diag .+ L._QP_HESS_REG[]))
function Cmat(z)
    C = Matrix(L.build_active_jacobian(z))
    B = zeros(length(z.active_bound_cols), z.n_vars)
    for (r, c) in enumerate(z.active_bound_cols); B[r, c] = 1.0; end
    return vcat(C, B)
end

function audit_cell(cell, ws, params; warm = 3, do_hess = false, do_danskin = true)
    a3 = ws.alpha3
    for it in 1:warm; half_iter!(ws, params, it); finish_iter!(ws, params, a3); end
    half_iter!(ws, params, warm + 1)

    Y = getY(ws); o = offsets(ws); lam = lam_vec(ws); nu = nu_vec(ws); mu = ws.mu[]
    r0 = rows_at(ws, Y)

    # ---- T_H residual consistency
    record!(cell, "H1 ‖rows‖₂ == r_primal", abs(norm(r0) - ws.r_primal) / max(1, ws.r_primal), 1e-10)
    l1 = ws isa L._ACPhase2Workspace ? L._ac_resid_l1(ws, acc_y) :
         ws isa L._DCPhase2Workspace ? L._dc_resid_l1(ws, acc_y) : L._socwr_resid_l1(ws, acc_y)
    record!(cell, "H2 ‖rows‖₁ == Φ L1 term", abs(norm(r0, 1) - l1) / max(1, l1), 1e-10)
    rnt = ws isa L._ACPhase2Workspace ? L._dual_residual(ws) : rows_nt(ws)
    # _dual_slope reads live JuMP values in DC/SOCWR == zone.y right after local_solve
    record!(cell, "H3 slope(r) == ‖r‖²", abs(L._dual_slope(ws, L._dual_residual(ws)) - dot(r0, r0)) / max(1, dot(r0, r0)), 1e-8)

    # ---- coupling Jacobian by FD (independent of hand derivatives)
    A = fd_jac(ws, Y)

    # ---- T_C local KKT
    Atl = A' * lam; worst = 0.0; scale = 0.0; worstM = 0.0; worstG = 0.0
    for (i, z) in enumerate(ws.zones)
        gi = grad_obj_base(ws.partitions.parts[i], z)      # independent ∇f
        R = gi .+ Atl[o[i]+1:o[i+1]] .+ prox_grad(ws, i) .+ constraint_term(z)
        # lb == ub columns: Ipopt (fixed_variable_treatment=make_parameter) removes them and
        # reports NO bound multiplier, so their KKT row cannot be closed from JuMP duals.
        # They are always in active_bound_cols (Δy = 0), hence irrelevant to the QP.
        pinned = [k for k in 1:z.n_vars if (!isnan(z.lb[k]) && z.lb[k] == z.ub[k]) || JuMP.is_fixed(z.all_vars[k])]
        @assert all(k -> k in z.active_bound_cols || JuMP.is_fixed(z.all_vars[k]), pinned) "pinned col not in active set"
        R[pinned] .= 0.0
        worst = max(worst, norm(R, Inf)); scale = max(scale, norm(gi, Inf), norm(Atl, Inf))
        if norm(R, Inf) > 1e-4 && get(ENV, "AUDIT_DIAG", "0") == "1"
            idx = sortperm(abs.(R); rev = true)[1:min(6, end)]
            println("  [KKT diag] $cell zone $i: ", [(JuMP.name(z.all_vars[k]), round(R[k]; sigdigits = 4),
                     "g=", round(gi[k]; sigdigits = 3), "Atl=", round(Atl[o[i]+k]; sigdigits = 3),
                     "prox=", round(prox_grad(ws, i)[k]; sigdigits = 3), "con=", round(constraint_term(z)[k]; sigdigits = 3),
                     "fixed=", JuMP.is_fixed(z.all_vars[k])) for k in idx])
            println("  term=", JuMP.termination_status(z.local_model), " raw=", JuMP.raw_status(z.local_model))
        end
        if L._MODIFIED_GRADIENT[] && z isa L.AladinZoneDataAC
            Rm = z.g .+ Atl[o[i]+1:o[i+1]] .+ prox_grad(ws, i) .+ constraint_term(z; which = :strong)
            Rm[pinned] .= 0.0
            worstM = max(worstM, norm(Rm, Inf))
        else
            worstG = max(worstG, norm(z.g .- gi, Inf))
        end
    end
    record!(cell, "C  local KKT ‖·‖∞ / scale", worst / max(1, scale), ws isa L._SOCWRPhase2Workspace ? 1e-3 : 5e-5)
    if L._MODIFIED_GRADIENT[] && ws isa L._ACPhase2Workspace
        record!(cell, "M  g_mod + Aᵀλ + prox + C̃ᵀκ = 0", worstM / max(1, scale), 5e-5)
    else
        record!(cell, "C2 zone.g == ∇(obj_base)", worstG / max(1, scale), 1e-10)
    end

    # ---- T_B coordination QP by null-space method
    n = length(Y); H = zeros(n, n); g = zeros(n); Zs = Matrix{Float64}[]
    for (i, z) in enumerate(ws.zones)
        H[o[i]+1:o[i+1], o[i]+1:o[i+1]] = Hmat(z); g[o[i]+1:o[i+1]] = z.g
        push!(Zs, nullspace(Cmat(z); atol = 1e-9))
    end
    Z = zeros(n, sum(size(Zi, 2) for Zi in Zs)); c0 = 0
    for (i, Zi) in enumerate(Zs); Z[o[i]+1:o[i+1], c0+1:c0+size(Zi, 2)] = Zi; c0 += size(Zi, 2); end
    AZ = A * Z
    M = Symmetric(Z' * H * Z + mu * (AZ' * AZ)); rhs = -(Z' * g + AZ' * (lam + mu * r0))
    w = size(Z, 2) == 0 ? Float64[] : (M \ rhs)
    dy_ref = Z * w; s_ref = r0 + A * dy_ref; nu_ref = lam + mu * s_ref
    dy_lir = reduce(vcat, [copy(z.delta) for z in ws.zones])
    qobj(dy) = (s_ = r0 + A * dy; 0.5 * dot(dy, H * dy) + dot(g, dy) + dot(lam, s_) + 0.5 * mu * dot(s_, s_))
    # Δy is compared entry-wise only where the reduced QP is well conditioned. With
    # Hᵢ = ∇²fᵢ + 1e-6·I (DC / no-prox) the reduced Hessian has ~1e-6 eigenvalues along
    # directions that change neither the cost nor A·Δy; there the minimiser is
    # (local-solver noise)/1e-6 and is not a meaningful reference. The solver-independent
    # optimality statement is the REDUCED GRADIENT of Lirion's Δy (B0).
    w_lir = Z' * dy_lir
    record!(cell, "B0 reduced-gradient of Lirion Δy (rel)", size(Z, 2) == 0 ? 0.0 : norm(M * w_lir - rhs, Inf) / max(1.0, norm(rhs, Inf)), 1e-5)
    wellcond = size(Z, 2) == 0 || cond(Matrix(M)) < 1e9
    wellcond && record!(cell, "B1 Δy  QP(Lirion) vs dense", norm(dy_lir - dy_ref, Inf) / max(1e-3, norm(dy_ref, Inf)), 2e-3)
    record!(cell, "B1' QP objective(Lirion) − optimum (rel)", (qobj(dy_lir) - qobj(dy_ref)) / max(1.0, abs(qobj(dy_ref))), 1e-6)
    if norm(dy_lir - dy_ref, Inf) / max(1e-3, norm(dy_ref, Inf)) > 2e-3 && get(ENV, "AUDIT_DIAG", "0") == "1"
        for (i, z) in enumerate(ws.zones)
            d = (dy_lir - dy_ref)[o[i]+1:o[i+1]]
            th = hasproperty(z, :theta_cols) ? z.theta_cols : Int[]
            isempty(th) && (th = [k for k in 1:z.n_vars if startswith(JuMP.name(z.all_vars[k]), "theta")])
            rest = setdiff(1:z.n_vars, th)
            println("  [B1 diag] $cell zone $i: Δθ-diff min/max = ", extrema(d[th]), "  non-θ diff ∞ = ", isempty(rest) ? 0.0 : norm(d[rest], Inf),
                    "  ‖dy_ref‖∞=", norm(dy_ref[o[i]+1:o[i+1]], Inf), "  ‖dy_lir‖∞=", norm(dy_lir[o[i]+1:o[i+1]], Inf))
        end
    end
    record!(cell, "B1'' A·Δy(Lirion) vs dense", norm(A * (dy_lir - dy_ref), Inf) / max(1e-3, norm(A * dy_ref, Inf)), 2e-3)
    record!(cell, "B2 λ_QP QP(Lirion) vs dense", norm(nu - nu_ref, Inf) / max(1, norm(nu_ref, Inf)), 2e-3)
    s_lir = r0 + A * dy_lir
    record!(cell, "B3 λ_QP == λ + μ s (fn.4)", norm(nu - (lam + mu * s_lir), Inf) / max(1, norm(nu, Inf)), 2e-3)
    record!(cell, "B4 CΔy == 0", maximum([norm(Cmat(z) * z.delta, Inf) for z in ws.zones]), 1e-5)

    # ---- T_D AC Hessian vs FD of ∇L
    if do_hess && ws isa L._ACPhase2Workspace
        worstH = 0.0
        for z in ws.zones
            ev, rows, strong = L._build_ac_evaluator(z, L._ACTIVE_DUAL_TOL[])
            Hl = Matrix(L._extract_lagrangian_hessian_ac(ev, rows, z.y, z.n_vars))
            kap = Float64[-JuMP.dual(ci) for ci in rows]
            function gradL(y)
                gf = zeros(z.n_vars); isempty(z.pg_cols) || all(==(0.0), vcat(z.pg_c2_sbase2, z.pg_c1_sbase)) || MOI.eval_objective_gradient(ev, gf, y)
                jt = zeros(z.n_vars); isempty(rows) || MOI.eval_constraint_jacobian_transpose_product(ev, jt, y, kap)
                return gf .+ jt
            end
            Hfd = zeros(z.n_vars, z.n_vars); h = 1e-6
            for j in 1:z.n_vars
                yp = copy(z.y); ym = copy(z.y); yp[j] += h; ym[j] -= h
                Hfd[:, j] = (gradL(yp) .- gradL(ym)) ./ (2h)
            end
            worstH = max(worstH, norm(Hl - Hfd, Inf) / max(1, norm(Hfd, Inf)))
            worstH = max(worstH, norm(Hl - Hl', Inf))
        end
        record!(cell, "D  ∇²L vs FD (rel)", worstH, 1e-5)
    end

    # ---- T_G Danskin slope of V_ρ(x, λ + α d), d = λ_QP − λ  (Alg. 3 step (c))
    # DC zones default to OSQP (first-order): the FD of V_ρ is noise-limited there;
    # the DC Danskin slope has its own test in test_aladin_dc.jl.
    if do_danskin && !(ws isa L._DCPhase2Workspace)
        states = collect(values(ws.vars)); saved = [copy(s.lambda) for s in states]
        savedy = [copy(z.y) for z in ws.zones]
        dirs = map(s -> s.nu .- s.lambda, ws.vars); dv = collect(values(dirs))
        dn = sqrt(sum(sum(abs2, d) for d in dv; init = 0.0))
        if dn > 0
            slope = L._dual_slope(ws, dirs)
            eps = 0.05
            function V(a)
                for (s, s0, d) in zip(states, saved, dv); @. s.lambda = s0 + a * d; end
                return L._vrho_resolve!(ws)
            end
            vp = V(eps); vm = V(-eps); fd = (vp - vm) / (2eps)
            for (s, s0) in zip(states, saved); copyto!(s.lambda, s0); end
            for (z, y0) in zip(ws.zones, savedy); copyto!(z.y, y0); end
            noise = 1e-7 * max(1.0, abs(vp)) / eps      # solver objective accuracy / step
            get(ENV, "AUDIT_DIAG", "0") == "1" && println("  [Danskin] $cell slope=$slope fd=$fd noise=$noise")
            record!(cell, "G  Danskin slope vs FD(V_ρ)", max(0.0, abs(slope - fd) - 3noise) / max(noise, abs(slope), abs(fd)), 5e-2)
        end
    end
    return nothing
end

# ------------------------------------------------------------
# SOCWR cone normal / distance
# ------------------------------------------------------------
function audit_cones()
    for (set, q) in ((MOI.SecondOrderCone(4), v -> sum(abs2, v[2:end]) - v[1]^2),
                     (MOI.RotatedSecondOrderCone(4), v -> sum(abs2, v[3:end]) - 2v[1] * v[2]))
        v = [1.3, 0.7, -0.4, 0.9]; nrm = L._cone_active_normal(set, v); fd = similar(v)
        for k in eachindex(v); e = zeros(4); e[k] = 1e-6; fd[k] = (q(v + e) - q(v - e)) / 2e-6; end
        record!("cones", "N  normal vs FD $(typeof(set).name.name)", norm(nrm - fd, Inf), 1e-6)
    end
    record!("cones", "N  soc_dist interior == 0", L._soc_dist([2.0, 1.0, 1.0]), 0.0)
    record!("cones", "N  soc_dist formula", abs(L._soc_dist([0.0, 3.0, 4.0]) - 5 / sqrt(2)), 1e-12)
    # RSOC point (t,u,x): 2tu ≥ x²; (1,1,√2) on boundary
    record!("cones", "N  rsoc_dist boundary == 0", L._rsoc_dist([1.0, 1.0, sqrt(2.0)]), 1e-12)
end


function run_all(tdir::AbstractString)
    empty!(RESULTS)
    c14 = joinpath(tdir, "pglib_opf_case14_ieee.m")
    star = joinpath(tdir, "case4_3zone_star.m")
    cells = Any[]
    for (mname, model) in (("dc", L.dc()), ("ac", L.ac()), ("socwr", L.socwr()))
        for coupling in (:state, :power, :hybrid), sect in (:tie_line, :dummy_bus, :split_bus)
            coupling === :state && sect !== :tie_line && continue
            coupling === :power && sect === :split_bus && continue
            mname == "socwr" && sect === :dummy_bus && continue
            push!(cells, (mname, model, coupling, sect, c14, (;), false))
        end
    end
    push!(cells, ("ac", L.ac(), :hybrid, :split_bus, star, (;), false))
    push!(cells, ("dc", L.dc(), :hybrid, :split_bus, star, (;), false))
    push!(cells, ("ac", L.ac(), :hybrid, :dummy_bus, c14, (; proximal_mode = :pflow), false))
    push!(cells, ("dc", L.dc(), :power, :tie_line, c14, (; proximal_mode = :pflow), false))
    # canonical AC configuration (no prox in B + modified gradient)
    for (coupling, sect) in ((:state, :tie_line), (:power, :tie_line), (:hybrid, :dummy_bus))
        push!(cells, ("ac", L.ac(), coupling, sect, c14, (;), true))
    end
    saved = (L._AC_QP_NO_PROX[], L._MODIFIED_GRADIENT[], L._BOUND_DETECT_MODE[])
    try
        audit_cones()
        for (mname, model, coupling, sect, file, kw, canon) in cells
            cell = string(mname, "/", coupling, "/", sect, "/", splitext(basename(file))[1],
                          isempty(kw) ? "" : "/pflow", canon ? "/canon" : "")
            L._AC_QP_NO_PROX[] = canon; L._MODIFIED_GRADIENT[] = canon
            L._BOUND_DETECT_MODE[] = canon ? :primal_dual : :primal
            redirect_stdout(devnull) do
                ws, params = build(file, model, coupling, sect; rho = 100.0, mu = 1000.0,
                                   globalization = :full_step, kw...)
                audit_cell(cell, ws, params; do_hess = mname == "ac")
            end
        end
    finally
        L._AC_QP_NO_PROX[], L._MODIFIED_GRADIENT[], L._BOUND_DETECT_MODE[] = saved
    end
    return RESULTS
end

end # module

function run_aladin_invariant_tests()
    @testset "ALADIN mathematical invariants (HFD Alg. 2 / Alg. 3)" begin
        for r in AladinInvariants.run_all(@__DIR__)
            r.ok || @warn "ALADIN invariant violated" r.cell r.test r.val r.tol
            @test r.ok
        end
    end
    return nothing
end
