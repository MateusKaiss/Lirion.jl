# Basic ALADIN (Algorithm 1) — versão em Julia
# Tradução fiel do código Python fornecido pelo Prof. Laercio
# Autor(es) originais (Python): Prof. Laercio (adaptado)

using LinearAlgebra

# ------------------------------- Tipos básicos -------------------------------

struct LocalFunction
    """
    Representa f_i e (opcionalmente) seu gradiente e Hessiana.
    - dim: dimensão de x_i
    - f   : x_i -> escalar
    - grad: x_i -> vetor (ou nothing)
    - hess: x_i -> matriz (ou nothing)
    """
    dim::Int
    f::Function
    grad::Union{Function,Nothing}
    hess::Union{Function,Nothing}
end

LocalFunction(dim, f; grad=nothing, hess=nothing) = LocalFunction(dim, f, grad, hess)

Base.@kwdef struct AladinOptions
    max_iters::Int    = 200
    tol::Float64      = 1e-6
    sub_max_iters::Int = 200
    sub_tol::Float64   = 1e-8
    bt_beta::Float64   = 0.5
    bt_c::Float64      = 1e-4
end

mutable struct BasicALADIN
    """
    Implementação do Algorithm 1 (Basic ALADIN).
    - f_list: vetor de LocalFunction (f_i, ∇f_i, ∇²f_i)
    - A_list: vetor de matrizes A_i com shape (m, n_i)
    - b: vetor do acoplamento (m,)
    - Sigma_list: vetor de matrizes Σ_i ≻ 0 (se nada, usa identidade)
    """
    f_list::Vector{LocalFunction}
    A_list::Vector{Matrix{Float64}}
    b::Vector{Float64}
    m::Int
    N::Int
    Sigma_list::Vector{Matrix{Float64}}
    opts::AladinOptions
end

function BasicALADIN(
    f_list::Vector{LocalFunction},
    A_list::Vector{Matrix{Float64}},
    b::Vector{Float64};
    Sigma_list::Union{Vector{Matrix{Float64}},Nothing}=nothing,
    opts::AladinOptions=AladinOptions()
)
    m = length(b)
    N = length(f_list)
    Σ_list = Sigma_list === nothing ? [I*1.0 for fi in f_list] : Sigma_list

    @assert length(A_list) == N "A_list deve ter N elementos"
    @assert length(Σ_list) == N "Sigma_list deve ter N elementos"

    # Checagens rápidas de conformidade e SPD de Σ_i
    for (i, (fi, Ai, Si)) in enumerate(zip(f_list, A_list, Σ_list))
        @assert size(Ai) == (m, fi.dim) "A_$(i) com shape inválido: $(size(Ai)) ≠ ($(m), $(fi.dim))"
        # checagem SPD simples
        try
            cholesky(Symmetric(Si); check=true)
        catch
            error("Σ_$i não é SPD")
        end
    end

    return BasicALADIN(f_list, A_list, b, m, N, Σ_list, opts)
end

# -------------------- Utilitários numéricos internos --------------------

_norm_S(y_minus_x::AbstractVector, S::AbstractMatrix) = sqrt(dot(y_minus_x, S * y_minus_x))

function _regularize_SPD(H::AbstractMatrix)
    # Garante SPD por regularização espectral leve
    Hsym = Symmetric((H + H')/2)  # força simetria numérica
    try
        cholesky(Hsym; check=true)
        return Matrix(Hsym)
    catch
        E = eigen(Hsym)
        vals = map(v -> max(v, 1e-8), E.values)
        return E.vectors * Diagonal(vals) * E.vectors'
    end
end

# ----------------- Solver de subproblema local (gradiente + BT) -----------------

function _solve_local_subproblem!(
    M::BasicALADIN,
    i::Int,
    x_i::Vector{Float64},
    lam::Vector{Float64}
)
    fi = M.f_list[i]
    Ai = M.A_list[i]
    Si = M.Sigma_list[i]  # matriz de ponderação

    F(y) = fi.f(y) + dot(lam, Ai * y) + 0.5 * dot(y .- x_i, Si * (y .- x_i))

    function gradF(y)
        if fi.grad === nothing
            error("grad f_i necessário para a variante básica; use deriv-free se não tiver grad.")
        end
        return fi.grad(y) + Ai' * lam + Si * (y .- x_i)
    end

    # inicialização
    y = copy(x_i)

    # descida por gradiente com backtracking
    for _ in 1:M.opts.sub_max_iters
        g = gradF(y)
        ng = norm(g)
        if ng <= M.opts.sub_tol
            break
        end
        d = -g
        t = 1.0
        f_y = F(y)
        while true
            y_new = y .+ t .* d
            if F(y_new) <= f_y + M.opts.bt_c * t * dot(g, d)
                y = y_new
                break
            end
            t *= M.opts.bt_beta
            if t < 1e-16
                y = y_new
                break
            end
        end
    end

    gi = fi.grad === nothing ? nothing : fi.grad(y)
    Hi = fi.hess === nothing ? nothing : fi.hess(y)
    if Hi === nothing
        Hi = I * 1.0
        Hi = Matrix(Hi, fi.dim, fi.dim)
    end
    Hi = _regularize_SPD(Hi)

    return y, gi, Hi
end

# ------------------------------- Loop principal -------------------------------

function solve(
    M::BasicALADIN,
    x0_list::Vector{Vector{Float64}};
    lam0::Union{Vector{Float64},Nothing}=nothing,
    callback::Union{Function,Nothing}=nothing
)
    @assert length(x0_list) == M.N
    x_list = [copy(xi) for xi in x0_list]
    lam = lam0 === nothing ? zeros(M.m) : copy(lam0)

    info = Dict{String,Any}("nit"=>0, "converged"=>false, "hist_res"=>Float64[])

    for k in 1:M.opts.max_iters
        # 1) subproblemas locais
        y_list = Vector{Vector{Float64}}()
        g_list = Vector{Union{Vector{Float64},Nothing}}()
        H_list = Vector{Matrix{Float64}}()
        for i in 1:M.N
            y_i, g_i, H_i = _solve_local_subproblem!(M, i, x_list[i], lam)
            push!(y_list, y_i)
            push!(g_list, g_i)
            push!(H_list, H_i)
        end

        # 2) critério de parada: max_i ||y_i - x_i||_{Σ_i}
        res = maximum(_norm_S(y_list[i] .- x_list[i], M.Sigma_list[i]) for i in 1:M.N)
        push!(info["hist_res"], res)
        if callback !== nothing
            callback(k, x_list, y_list, lam)
        end
        if res <= M.opts.tol
            info["nit"] = k
            info["converged"] = true
            return x_list, lam, info
        end

        # 3) QP acoplado via forma fechada
        #    Mmat = sum A_i H_i^{-1} A_i^T
        #    R    = sum A_i (y_i - H_i^{-1} g_i) - b
        #    lam⁺ = Mmat^{-1} R
        #    Δy_i = - H_i^{-1}(g_i + A_i^T lam⁺)
        Mmat = zeros(M.m, M.m)
        R = -copy(M.b)

        for i in 1:M.N
            Ai = M.A_list[i]
            Hi = H_list[i]
            gi = g_list[i]
            invHi_Ait = Hi \ (Ai')               # (n_i x m)
            Mmat .+= Ai * invHi_Ait                # (m x m)
            invHi_gi = gi === nothing ? zeros(size(Ai,2)) : (Hi \ gi)
            R .+= Ai * (y_list[i] .- invHi_gi)
        end

        lam_plus = Mmat \ R

        # 4) passo primal
        new_x_list = Vector{Vector{Float64}}(undef, M.N)
        for i in 1:M.N
            Ai = M.A_list[i]
            Hi = H_list[i]
            gi = g_list[i]
            rhs = (gi === nothing ? zeros(size(Ai,2)) : gi) .+ (Ai' * lam_plus)
            delta_y_i = -(Hi \ rhs)
            new_x_list[i] = y_list[i] .+ delta_y_i
        end

        x_list = new_x_list
        lam = lam_plus
    end

    info["nit"] = M.opts.max_iters
    return x_list, lam, info
end

# ------------------------- Exemplo mínimo (quadrático) -------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    # min 0.5*q1*x1^2 + 0.5*q2*(x2-1)^2  s.t.  x1 - x2 = 0
    q1, q2 = 0.1, 10.0
    m = 1
    A1 = [1.0;;]         # (1 x 1)
    A2 = [-1.0;;]        # (1 x 1)
    b  = [0.0]

    f1 = LocalFunction(1, x -> 0.5*q1*(x[1]^2);
        grad = x -> [q1*x[1]],
        hess = x -> [q1;;]
    )

    f2 = LocalFunction(1, x -> 0.5*q2*((x[1]-1.0)^2);
        grad = x -> [q2*(x[1]-1.0)],
        hess = x -> [q2;;]
    )

    Σ1 = [1.0;;]
    Σ2 = [1.0;;]

    model = BasicALADIN([f1, f2], [A1, A2], b; Sigma_list=[Σ1, Σ2], opts=AladinOptions(max_iters=100, tol=1e-10))

    x0 = [[0.0], [0.0]]
    lam0 = [0.0]
    x_sol, lam_sol, info = solve(model, x0; lam0=lam0)

    x1, x2 = x_sol[1][1], x_sol[2][1]
    println("x1, x2, λ = ", (x1, x2, lam_sol[1]))
    println("x* = ", q2/(q1+q2), "  λ* = ", -q1*q2/(q1+q2))
    println("converged = ", info["converged"], "  iters = ", info["nit"]) 
end
