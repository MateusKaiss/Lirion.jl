# =========================
# ALADIN - Augmented Lagrangian based Alternating Direction Inexact Newton
# =========================

"""
    AladinOptions

Opções de configuração para o algoritmo ALADIN.

# Campos
- `max_iter::Int`: Número máximo de iterações (padrão: 100)
- `tol::Float64`: Tolerância para convergência (padrão: 1e-4)
- `rho::Float64`: Parâmetro de penalidade para o Lagrangiano aumentado (padrão: 1.0)
- `mu::Float64`: Parâmetro de regularização para o QP (padrão: 1e-3)
- `verbose::Bool`: Exibir informações de progresso (padrão: true)
- `hessian_approximation::Symbol`: Tipo de aproximação Hessiana (:exact, :gauss_newton, :bfgs) (padrão: :gauss_newton)
"""
Base.@kwdef mutable struct AladinOptions
    max_iter::Int = 100
    tol::Float64 = 1e-4
    rho::Float64 = 1.0
    mu::Float64 = 1e-3
    verbose::Bool = true
    hessian_approximation::Symbol = :gauss_newton
end

"""
    AladinRegion

Representa uma região (subproblema) no problema de otimização distribuída.

# Campos
- `id::Int`: Identificador da região (1, 2, 3...)
- `data::Dict`: com os dados daquela região (barramentos, geradores, ramos, etc.).
- `coupling_vars::Vector{Symbol}`: variáveis dessa região que participam das restrições de consenso (por exemplo, ângulos/voltagens em barras de fronteira).
- `A_consensus::Matrix{Float64}`: matriz local de consenso (Aᵢ); pode ser usada para compor a matriz global ∑Aᵢxᵢ = b.

O construtor simplificado (AladinRegion(id, data)) só preenche id e data e deixa coupling_vars vazio 
e A_consensus = nothing. A divisão da rede em regiões e a montagem de A_consensus não está neste 
arquivo – isso é feito “fora”, quando você monta o problema.

"""
mutable struct AladinRegion
    id::Int
    data::Dict
    coupling_vars::Vector{Symbol}
    A_consensus::Union{Nothing, Matrix{Float64}}
    
    function AladinRegion(id::Int, data::Dict)
        new(id, data, Symbol[], nothing)
    end
end

"""
    AladinProblem

Estrutura principal do problema ALADIN.

# Campos
- `regions::Vector{AladinRegion}`: Lista de regiões (subproblemas)  ( As regiões da rede)
- `consensus_matrix::Matrix{Float64}`: Matriz de consenso global ∑Aᵢxᵢ = b
- `consensus_rhs::Vector{Float64}`: Lado direito da restrição de consenso b  (vetor b da restrição de consenso)
- `model_type::Type{<:OPFModel}`: Tipo de modelo OPF (DC, AC, SOCP)
- `optimizer`: Otimizador a ser usado nos subproblemas (o solver de otimização usado pelo JuMP)
- `options::AladinOptions`: Opções do algoritmo

Variáveis de iteração:
lambda: 
x_current[i]: 
y_current[i]: solução local da região i (saída dos NLPs desacoplados).

O construtor cuida de inicializar lambda com zeros e cria vetores vazios para x_current/y_current.

"""
mutable struct AladinProblem
    regions::Vector{AladinRegion}
    consensus_matrix::Matrix{Float64}
    consensus_rhs::Vector{Float64}
    model_type::Type{<:OPFModel}
    optimizer
    options::AladinOptions
    
    # Variáveis de iteração
    lambda::Vector{Float64}  # vetor de multiplicadores de Lagrange associados às restrições de consenso.
    x_current::Vector{Vector{Float64}}  # soluçções primal “consistente” (após o QP) da região i.
    y_current::Vector{Vector{Float64}}  # local da região i (saída dos NLPs desacoplados).

    #O construtor cuida de inicializar lambda com zeros e cria vetores vazios para x_current/y_current.
    
    function AladinProblem(
        regions::Vector{AladinRegion},
        consensus_matrix::Matrix{Float64},
        consensus_rhs::Vector{Float64},
        model_type::Type{<:OPFModel},
        optimizer;
        options::AladinOptions = AladinOptions()
    )
        n_regions = length(regions)
        lambda = zeros(size(consensus_matrix, 1))
        x_current = [Float64[] for _ in 1:n_regions]
        y_current = [Float64[] for _ in 1:n_regions]
        
        new(regions, consensus_matrix, consensus_rhs, model_type, optimizer, 
            options, lambda, x_current, y_current)
    end
end

"""
    solve_aladin(problem::AladinProblem)

Resolve o problema de otimização distribuída usando o algoritmo ALADIN.

# Algoritmo
1. Resolver NLPs desacoplados com Lagrangiano aumentado
2. Computar gradientes e aproximações Hessianas
3. Resolver QP acoplado para atualização de consenso
4. Atualizar multiplicadores de Lagrange
5. Verificar convergência

# Retorna
- `results::Dict`: Dicionário com resultados da otimização

Esta é uma implementação do ALADIN para OPF distribído

"""

#função principal
function solve_aladin(problem::AladinProblem)
    opts = problem.options
    
    if opts.verbose
        println("="^70)
        println("ALADIN - Augmented Lagrangian based Alternating Direction Inexact Newton")
        println("="^70)
        println("Número de regiões: ", length(problem.regions))
        println("Tipo de modelo: ", problem.model_type)
        println("Tolerância: ", opts.tol)
        println("Máximo de iterações: ", opts.max_iter)
        println("="^70)
    end
    
    """
    """

    # Inicialização
    initialize_aladin!(problem)   # Resolve cada região isoladamente (sem acoplamento) para obter um ponto inicial x_current[i] e y_current[i].
    
    # Loop principal do ALADIN
    converged = false
    iteration = 0
    
    for k in 1:opts.max_iter    #loop principal
        iteration = k
        
        if opts.verbose
            println("\n--- Iteração $k ---")
        end
        
        # Passo 1: Resolver NLPs desacoplados
        # Cada região resolve seu subproblema com Lagrangiano aumentado.
        solve_decoupled_nlps!(problem)
        
        # Passo 2: Calcular gradientes e Hessianas
        gradients, hessians = compute_derivatives!(problem)
        
        # Passo 3: Resolver QP acoplado
        # Esse QP coordena todas as regiões para satisfazer o consenso.
        solve_coupled_qp!(problem, gradients, hessians)
        
        # Passo 4: Atualizar multiplicadores
        update_multipliers!(problem)
        
        # Passo 5: Verificar convergência
        # Convergência quando max_i ||yᵢ - xᵢ|| < tol.
        residual = check_convergence(problem)
        
        if opts.verbose
            println("Resíduo de consenso: ", residual)
        end
        
        #Se atingir tolerância → imprime mensagem de convergência
        if residual < opts.tol
            converged = true
            if opts.verbose
                println("\n✓ Convergência alcançada!")
            end
            break
        end
    end
    
    # Se NÃO atingir tolerância → imprime mensagem que chegou em max_iter.
    if !converged && opts.verbose
        println("\n⚠ Número máximo de iterações atingido sem convergência.")
    end
    
    # Compilar resultados
    results = compile_results(problem, converged, iteration)
    
    return results
end

"""
    initialize_aladin!(problem::AladinProblem)

Inicializa as variáveis do ALADIN resolvendo os subproblemas desacoplados.
"""
function initialize_aladin!(problem::AladinProblem)
    if problem.options.verbose
        println("\nInicializando ALADIN...")
    end
    
    # Inicializar multiplicadores com zeros  (lambda)
    problem.lambda = zeros(size(problem.consensus_matrix, 1))
    
    # Resolver cada região independentemente para obter ponto inicial
    # Para cada região:
    #       Chama build_and_optimize_model (definida em outro lugar) para montar e resolver o OPF local.
    #       Extrai a solução primal com extract_primal_solution.
    #       Usa essa solução como ponto inicial tanto para x_current[i] quanto para y_current[i].


    for (i, region) in enumerate(problem.regions)
        model = build_and_optimize_model(
            region.data,
            problem.model_type;
            optimizer = problem.optimizer,
            build_only = false
        )
        
        # Extrair solução primal
        problem.x_current[i] = extract_primal_solution(model, region)
        problem.y_current[i] = copy(problem.x_current[i])
    end
    
    if problem.options.verbose
        println("Inicialização concluída.")
    end
end

"""
    solve_decoupled_nlps!(problem::AladinProblem)

Resolve os NLPs desacoplados com termo de Lagrangiano aumentado.
Minimiza: fᵢ(yᵢ) + λᵀAᵢyᵢ + (1/2)||yᵢ - xᵢ||²   (???_Σᵢ)
"""
function solve_decoupled_nlps!(problem::AladinProblem)
    n_regions = length(problem.regions)
    
    # Resolver em paralelo (se disponível)
    for i in 1:n_regions
        region = problem.regions[i]
        
        # Construir modelo com Lagrangiano aumentado
        # Esse mecanismo é um esqueleto, que será implementado a adição do tempo de Langragiano 
        #     aumento na função objetivo
        # Se o solver local converge, a nova solução local y_current[i] é atualizada.
        model = build_augmented_lagrangian_model(
            region,
            problem.model_type,
            problem.optimizer,
            problem.lambda,
            problem.x_current[i],
            problem.options.rho
        )
        
        # Otimizar
        JuMP.optimize!(model)
        
        # Extrair solução
        if JuMP.termination_status(model) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
            problem.y_current[i] = extract_primal_solution(model, region)
        else
            @warn "Região $i não convergiu. Status: $(JuMP.termination_status(model))"
        end
    end
end

"""
    compute_derivatives!(problem::AladinProblem)

Computa gradientes e aproximações Hessianas para cada subproblema.
"""
function compute_derivatives!(problem::AladinProblem)
    n_regions = length(problem.regions)
    gradients = Vector{Vector{Float64}}(undef, n_regions)
    hessians = Vector{Matrix{Float64}}(undef, n_regions)
    
    for i in 1:n_regions
        region = problem.regions[i]
        y_i = problem.y_current[i]
        
        # Computar gradiente: ∇fᵢ(yᵢ) + Aᵢᵀλ
        # Gradiente do Lagrangiano aumentado (ainda retorna zeros).
        gradients[i] = compute_gradient(region, y_i, problem)
        
        # Computar aproximação Hessiana
        # Hessiana (exata, Gauss-Newton ou BFGS). No momento, todos esses são placeholders:
        #       compute_gauss_newton_hessian: retorna identidade.
        #       compute_exact_hessian: retorna identidade.
        #       compute_bfgs_hessian: também identidade.
        if problem.options.hessian_approximation == :exact
            hessians[i] = compute_exact_hessian(region, y_i)
        elseif problem.options.hessian_approximation == :gauss_newton
            hessians[i] = compute_gauss_newton_hessian(region, y_i)
        elseif problem.options.hessian_approximation == :bfgs
            hessians[i] = compute_bfgs_hessian(region, y_i, problem)
        else
            error("Aproximação Hessiana desconhecida: $(problem.options.hessian_approximation)")
        end
    end
    
    return gradients, hessians
end

"""
    solve_coupled_qp!(problem::AladinProblem, gradients, hessians)

Resolve o QP acoplado:
min_Δy ∑ᵢ {(1/2)ΔyᵢᵀHᵢΔyᵢ + gᵢᵀΔyᵢ}
s.t. ∑ᵢ Aᵢ(yᵢ + Δyᵢ) = b


"""
function solve_coupled_qp!(problem::AladinProblem, gradients, hessians)
    n_regions = length(problem.regions)
    
    # Criar modelo QP no JuMP
    qp_model = JuMP.Model(problem.optimizer)
    
    # Variáveis: incrementos Δy para cada região
    # Para cada região i, cria variáveis de incremento Δyᵢ
    n_vars = [length(problem.y_current[i]) for i in 1:n_regions]
    @variable(qp_model, Δy[i=1:n_regions, j=1:n_vars[i]])
    
    # Função objetivo quadrática
    #  
    @objective(qp_model, Min,
        sum(
            0.5 * sum(Δy[i,j] * sum(hessians[i][j,k] * Δy[i,k] for k in 1:n_vars[i]) 
                     for j in 1:n_vars[i]) +
            sum(gradients[i][j] * Δy[i,j] for j in 1:n_vars[i])
            for i in 1:n_regions
        )
    )
    
    # Restrição de consenso: ∑ᵢ Aᵢ(yᵢ + Δyᵢ) = b
    # consensus_matrix é a matriz A da restrição de consenso
    # get_global_index(i, j, n_vars) mapeia o índice local (região i, variável j) para o índice global em x_all = vcat(x₁, x₂, ..., x_R)
    n_consensus = size(problem.consensus_matrix, 1)
    @constraint(qp_model, consensus[c=1:n_consensus],
        sum(
            sum(problem.consensus_matrix[c, get_global_index(i, j, n_vars)] * 
                (problem.y_current[i][j] + Δy[i,j])
                for j in 1:n_vars[i])
            for i in 1:n_regions
        ) == problem.consensus_rhs[c]
    )
    
    # Resolver QP
    JuMP.optimize!(qp_model)
    
    # Atualizar x = y + Δy
    # se o QP convergir, atualiza x  (verificar se está atualizando corretamente)
    if JuMP.termination_status(qp_model) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        for i in 1:n_regions
            problem.x_current[i] = problem.y_current[i] + JuMP.value.(Δy[i,:])
        end
    else
        @warn "QP acoplado não convergiu. Status: $(JuMP.termination_status(qp_model))"
    end
end

"""
    update_multipliers!(problem::AladinProblem)

Atualiza os multiplicadores de Lagrange: λ⁺ = λ + ρ(∑Aᵢxᵢ - b)
"""
function update_multipliers!(problem::AladinProblem)
    # Calcular violação de consenso
    consensus_violation = compute_consensus_violation(problem)
    
    # Atualizar multiplicadores
    problem.lambda += problem.options.rho * consensus_violation
end

"""
    check_convergence(problem::AladinProblem)

Verifica convergência medindo ||yᵢ - xᵢ|| para todas as regiões.

Mede o desacoplamento entre as soluções locais (yᵢ) e a solução coordenada (xᵢ)
Se isso for pequeno em todas as regiões, considera que o algoritmo convergiu

"""

function check_convergence(problem::AladinProblem)
    max_residual = 0.0
    
    for i in 1:length(problem.regions)
        residual = norm(problem.y_current[i] - problem.x_current[i])
        max_residual = max(max_residual, residual)
    end
    
    return max_residual
end

# =========================
# Funções auxiliares
# =========================

"""
    extract_primal_solution(model::JuMP.Model, region::AladinRegion)

Extrai a solução primal do modelo JuMP como vetor.

Extrai as variáveis relevantes do modelo JuMP em um único vetor de Float64
Aqui assume que as variáveis do modelo se chamam :theta, :v_magnitude, :v_magnitude_sq, :p_generation, :q_generation
A ordem e o tamanho dessas variáveis precisam ser compatíveis com a forma como você monta a consensus_matrix

"""
function extract_primal_solution(model::JuMP.Model, region::AladinRegion)
    solution = Float64[]
    
    # Extrair variáveis relevantes (theta, v_magnitude, p_generation, etc.)
    if haskey(model, :theta)
        append!(solution, vec(JuMP.value.(model[:theta])))
    end
    
    if haskey(model, :v_magnitude)
        append!(solution, vec(JuMP.value.(model[:v_magnitude])))
    end
    
    if haskey(model, :v_magnitude_sq)
        append!(solution, vec(JuMP.value.(model[:v_magnitude_sq])))
    end
    
    if haskey(model, :p_generation)
        append!(solution, vec(JuMP.value.(model[:p_generation])))
    end
    
    if haskey(model, :q_generation)
        append!(solution, vec(JuMP.value.(model[:q_generation])))
    end
    
    return solution
end

"""
    compute_consensus_violation(problem::AladinProblem)

Calcula a violação da restrição de consenso: ∑Aᵢxᵢ - b
"""
function compute_consensus_violation(problem::AladinProblem)
    n_regions = length(problem.regions)
    n_consensus = length(problem.consensus_rhs)
    violation = zeros(n_consensus)
    
    # Concatenar todas as variáveis primais
    # Concatena x_current[1], x_current[2], ..., x_current[R] em um único vetor x_all
    x_all = reduce(vcat, problem.x_current)
    
    # Calcular Ax - b
    violation = problem.consensus_matrix * x_all - problem.consensus_rhs
    
    return violation
end

"""
    get_global_index(region_id::Int, local_index::Int, n_vars::Vector{Int})

    Converte o índice local (dentro da região) em índice global na concatenação x_all
"""
function get_global_index(region_id::Int, local_index::Int, n_vars::Vector{Int})
    offset = sum(n_vars[1:region_id-1])
    return offset + local_index
end

"""
    compile_results(problem::AladinProblem, converged::Bool, iterations::Int)

    Compila os resultados finais da otimização.

    Junta tudo num dicionário:
            Se convergiu ou não
            Quantas iterações
            Soluções por região
            Multiplicadores finais
            Norma da violação de consenso
"""
function compile_results(problem::AladinProblem, converged::Bool, iterations::Int)
    results = Dict{String, Any}()
    
    results["converged"] = converged
    results["iterations"] = iterations
    results["regions"] = length(problem.regions)
    
    # Soluções por região
    results["solutions"] = Dict{Int, Vector{Float64}}()
    for i in 1:length(problem.regions)
        results["solutions"][i] = problem.x_current[i]
    end
    
    # Multiplicadores finais
    results["lambda"] = problem.lambda
    
    # Violação de consenso final
    results["consensus_violation"] = norm(compute_consensus_violation(problem))
    
    return results
end

# =========================
# Funções específicas do modelo (a serem implementadas)
# =========================

"""
    build_augmented_lagrangian_model(region, model_type, optimizer, lambda, x_prev, rho)

Constrói modelo com Lagrangiano aumentado para o subproblema.
"""
function build_augmented_lagrangian_model(
    region::AladinRegion,
    model_type::Type{<:OPFModel},
    optimizer,
    lambda::Vector{Float64},
    x_prev::Vector{Float64},
    rho::Float64
)
    # Construir modelo base
    model = build_and_optimize_model(
        region.data,
        model_type;
        optimizer = optimizer,
        build_only = true
    )
    
    # Adicionar termos do Lagrangiano aumentado à função objetivo
    # TODO: Implementar adição dos termos λᵀAᵢyᵢ + (ρ/2)||yᵢ - xᵢ||²
    
    return model
end

"""
    compute_gradient(region, y, problem)

Computa o gradiente do Lagrangiano aumentado.
"""
function compute_gradient(
    region::AladinRegion,
    y::Vector{Float64},
    problem::AladinProblem
)
    # TODO: Implementar cálculo do gradiente
    # Por enquanto, retorna vetor zero
    return zeros(length(y))
end

"""
    compute_gauss_newton_hessian(region, y)

Computa aproximação Hessiana usando Gauss-Newton.
"""
function compute_gauss_newton_hessian(region::AladinRegion, y::Vector{Float64})
    n = length(y)
    # Aproximação diagonal simples
    return Matrix{Float64}(I, n, n)
end

"""
    compute_exact_hessian(region, y)

Computa Hessiana exata (se disponível).
"""
function compute_exact_hessian(region::AladinRegion, y::Vector{Float64})
    # TODO: Implementar cálculo da Hessiana exata
    n = length(y)
    return Matrix{Float64}(I, n, n)
end

"""
    compute_bfgs_hessian(region, y, problem)

Computa aproximação Hessiana usando BFGS.
"""
function compute_bfgs_hessian(
    region::AladinRegion,
    y::Vector{Float64},
    problem::AladinProblem
)
    # TODO: Implementar BFGS
    n = length(y)
    return Matrix{Float64}(I, n, n)
end
