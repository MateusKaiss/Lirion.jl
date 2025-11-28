# =========================
# ALADIN - Exemplo de Uso
# =========================

"""
Exemplo de uso do algoritmo ALADIN para otimização distribuída
de fluxo de potência ótimo (OPF).

Este script demonstra:
1. Como carregar dados de uma rede elétrica
2. Como decompor a rede em regiões
3. Como configurar e resolver usando ALADIN
4. Como analisar os resultados

Tem um exemplo centralizado e outro distribuído 

Características da rede:
6 barras
7 branches
2 geradores
barras com carga:  todas (exceto a 1)

Base de potência: 100 MVA

| Barra | Tipo                              | Tensão inicial | Demanda P (pu) | Demanda Q (pu) |
| ----- | --------------------------------- | -------------- | -------------- | -------------- |
| 1     | **Barra de referência** + gerador | 1.0            | 0.0            | 0.0            |
| 2     | Carga                             | 1.0            | 0.5            | 0.1            |
| 3     | Carga                             | 1.0            | 0.3            | 0.1            |
| 4     | Carga + gerador                   | 1.0            | 0.4            | 0.1            |
| 5     | Carga                             | 1.0            | 0.2            | 0.05           |
| 6     | Carga                             | 1.0            | 0.3            | 0.05           |

todos os limites de tensão
0.95 ≤ V ≤ 1.05 pu


os geradores, coloquei nas barras 1 e 4

| Gerador | Barra | Pmin–Pmax (pu) | Qmin–Qmax (pu) |
| ------- | ----- | -------------- | -------------- |
| G1      | 1     | 0.0 – 2.0      | -0.5 – 0.5     |
| G2      | 4     | 0.0 – 1.5      | -0.3 – 0.3     |


Custo de geração (quadrático):

G1: 10 · Pg

G2: 15 · Pg

como o G2 é mais caro → o OPF prefere acionar mais G1.


Linhas de transmissão
| Linha | De → Para | Capacidade (pu) | R    | X   |
| ----- | --------- | --------------- | ---- | --- |
| 1     | 1 → 2     | 1.0             | 0.01 | 0.1 |
| 2     | 1 → 3     | 1.0             | 0.01 | 0.1 |
| 3     | 2 → 3     | 0.8             | 0.01 | 0.1 |
| 4     | 2 → 4     | 1.0             | 0.01 | 0.1 |
| 5     | 3 → 5     | 0.8             | 0.01 | 0.1 |
| 6     | 4 → 5     | 1.0             | 0.01 | 0.1 |
| 7     | 5 → 6     | 0.8             | 0.01 | 0.1 |


No exemplo fiz uma distribuição simples e lógica (vou ver com o Daniel), mas acredito que 
seja o suficiente para este exemplo simples

   Região 1                         Região 2
───────────────               ─────────────────
 Barras: 1, 2, 3               Barras: 4, 5, 6
 Gerador em 1                  Gerador em 4
 Mais carga nas barras 2-3     Carga nas barra 5-6

ou seja, separei somente em duas regiões. E na barra 5 considerei como sendo 
a variável de consenso, onde as variáveis são compartilhadas (na barra 5 ocorre o consenso)

Com isso, cada região deve otimizar localmente:
| Região | Variáveis internas | Variáveis compartilhadas      |
| ------ | ------------------ | ----------------------------- |
| 1      | θ₁, θ₂, θ₃, Pg₁    | θ₅ (na fronteira 3–5/2–5)     |
| 2      | θ₄, θ₅, θ₆, Pg₂    | θ₅ (mas calculado localmente) |


Pelo que verifiquei no artigo, o aladin exige:
θ₅(Região 1) = θ₅(Região 2)

ou seja se não forem iguais vai haver violação de consenso, e o 
aladin vai interagir té convergir (verificar se é isso mesmo ???)

a principio, o propósito desse exemplo de teste seria para:
* considerar geração em ambas as regiões
* haver congestionamento possível entre elas
* uma única barra de fronte para ajudar a visualizar o consenso
* seria uma estrutura típica de VPP com dois controladores locais ??? (verificar se posso afirmar isso)
"""

using JuMP
using Ipopt

# Incluir arquivos do projeto
include("helpers.jl")
include("build_model.jl")
include("bus.jl")
include("branch.jl")
include("generator.jl")
include("battery.jl")
include("prosumer.jl")
include("active_power_balance.jl")
include("reactive_power_balance.jl")
include("objective_function.jl")
include("aladin.jl")
include("aladin_decomposition.jl")

# Definir tipos de modelo OPF
abstract type OPFModel end
struct DC <: OPFModel end
struct AC <: OPFModel end
struct SOCP <: OPFModel end

# Tipo para construção de subproblema
struct SubproblemBuild end

"""
    create_simple_test_case()

Cria um caso de teste simples para demonstração.
"""
function create_simple_test_case()
    # Criar dados de exemplo para uma rede pequena
    # Este é um exemplo simplificado - na prática, você carregaria dados reais
    
    matpower_data = Dict{String, Any}()
    
    # Parâmetros de barramento
    bus_params = (
        bus_ids = [1, 2, 3, 4, 5, 6],
        vmin = Dict(i => 0.95 for i in 1:6),
        vmax = Dict(i => 1.05 for i in 1:6),
        vstart = Dict(i => 1.0 for i in 1:6),
        initial_voltage_angle = Dict(i => 0.0 for i in 1:6),
        pd = Dict(1 => 0.0, 2 => 0.5, 3 => 0.3, 4 => 0.4, 5 => 0.2, 6 => 0.3),
        qd = Dict(1 => 0.0, 2 => 0.1, 3 => 0.1, 4 => 0.1, 5 => 0.05, 6 => 0.05),
        reference_bus = 1
    )
    
    # Parâmetros de geradores
    gen_params = (
        number_of_generators = 2,
        bus = [1, 4],
        pmin = [0.0, 0.0],
        pmax = [2.0, 1.5],
        qmin = [-0.5, -0.3],
        qmax = [0.5, 0.3],
        generator_status = [1, 1],
        plant_type = ["C", "C"]  # C = Convencional
    )
    
    # Parâmetros de linhas
    branch_params = (
        number_of_branches = 7,
        bus_from = [1, 1, 2, 2, 3, 4, 5],
        bus_to = [2, 3, 3, 4, 5, 5, 6],
        capacity = [1.0, 1.0, 0.8, 1.0, 0.8, 1.0, 0.8],
        resistance = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01],
        reactance = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
        b = [10.0, 10.0, 8.0, 10.0, 8.0, 10.0, 8.0],
        existing = ones(Int, 7),
        min_ang = fill(NaN, 7),
        max_ang = fill(NaN, 7),
        tap_angle = zeros(7),
        inv_tap = ones(7),
        g = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
        b_shunt = zeros(7),
        active = 1:7,
        bus_idx_map = Dict(i => i for i in 1:6)
    )
    
    # Parâmetros de custo
    objective_params = (
        generator_costs = Dict(
            1 => Dict(1 => [0.0, 10.0, 0.0], 2 => [0.0, 15.0, 0.0])
        ),
        generator_status = [1, 1],
        number_of_generators = 2,
        plant_type = ["C", "C"],
        sbase = 100.0
    )
    
    matpower_data["bus_params"] = bus_params
    matpower_data["generator_params"] = gen_params
    matpower_data["branch_params"] = branch_params
    matpower_data["objective_params"] = objective_params
    matpower_data["sbase"] = 100.0
    
    return matpower_data
end

"""
    run_centralized_opf(matpower_data, model_type)

Executa OPF centralizado para comparação.
"""
function run_centralized_opf(matpower_data, model_type)
    println("\n" * "="^70)
    println("EXECUTANDO OPF CENTRALIZADO")
    println("="^70)
    
    optimizer = Ipopt.Optimizer
    model = build_model(matpower_data, model_type; optimizer = optimizer)
    
    status = JuMP.termination_status(model)
    
    if status in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        obj_value = JuMP.objective_value(model)
        println("✓ OPF Centralizado convergiu")
        println("  Valor objetivo: $obj_value")
        return model, obj_value
    else
        println("✗ OPF Centralizado falhou: $status")
        return model, nothing
    end
end

"""
    run_distributed_aladin(matpower_data, model_type, n_regions)

Executa OPF distribuído usando ALADIN.
"""
function run_distributed_aladin(matpower_data, model_type, n_regions)
    println("\n" * "="^70)
    println("EXECUTANDO OPF DISTRIBUÍDO COM ALADIN")
    println("="^70)
    
    # Passo 1: Decompor a rede
    println("\nDecompondo a rede em $n_regions regiões...")
    regions, consensus_matrix, consensus_rhs = decompose_network(
        matpower_data,
        n_regions;
        method = :geographical
    )
    
    visualize_decomposition(regions)
    
    # Passo 2: Configurar opções do ALADIN
    aladin_opts = AladinOptions(
        max_iter = 50,
        tol = 1e-4,
        rho = 1.0,
        mu = 1e-3,
        verbose = true,
        hessian_approximation = :gauss_newton
    )
    
    # Passo 3: Criar problema ALADIN
    optimizer = Ipopt.Optimizer
    problem = AladinProblem(
        regions,
        consensus_matrix,
        consensus_rhs,
        model_type,
        optimizer;
        options = aladin_opts
    )
    
    # Passo 4: Resolver
    results = solve_aladin(problem)
    
    return results
end

"""
    compare_results(centralized_obj, distributed_results)

Compara resultados centralizados e distribuídos.
"""
function compare_results(centralized_obj, distributed_results)
    println("\n" * "="^70)
    println("COMPARAÇÃO DE RESULTADOS")
    println("="^70)
    
    if centralized_obj !== nothing
        println("Valor objetivo centralizado: $centralized_obj")
    end
    
    if distributed_results["converged"]
        println("ALADIN convergiu em $(distributed_results["iterations"]) iterações")
        println("Violação de consenso final: $(distributed_results["consensus_violation"])")
        
        if centralized_obj !== nothing
            # Calcular diferença relativa
            # TODO: Calcular valor objetivo do ALADIN
            println("\nNota: Comparação detalhada de valores objetivos requer implementação adicional")
        end
    else
        println("ALADIN não convergiu dentro do limite de iterações")
    end
    
    println("="^70)
end

"""
    main()

Função principal que executa o exemplo completo.
"""
function main()
    println("\n")
    println("╔" * "="^68 * "╗")
    println("║" * " "^15 * "ALADIN - Exemplo de Uso" * " "^30 * "║")
    println("║" * " "^10 * "Otimização Distribuída de Fluxo de Potência" * " "^15 * "║")
    println("╚" * "="^68 * "╝")
    
    # Configurações
    model_type = DC  # Pode ser DC, AC ou SOCP
    n_regions = 2    # Número de regiões para decomposição
    
    # Criar caso de teste
    println("\nCriando caso de teste...")
    matpower_data = create_simple_test_case()
    
    # Executar OPF centralizado
    centralized_model, centralized_obj = run_centralized_opf(
        matpower_data,
        model_type
    )
    
    # Executar ALADIN distribuído
    distributed_results = run_distributed_aladin(
        matpower_data,
        model_type,
        n_regions
    )
    
    # Comparar resultados
    compare_results(centralized_obj, distributed_results)
    
    println("\n✓ Exemplo concluído!\n")
    
    return centralized_model, distributed_results
end

# Executar exemplo se este arquivo for executado diretamente
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
