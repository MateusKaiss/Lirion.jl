# =========================
# ALADIN - Funções de Decomposição de Rede
# =========================

"""
    decompose_network(matpower_data::Dict, n_regions::Int; method::Symbol=:metis)

Decompõe a rede elétrica em regiões para otimização distribuída.

# Argumentos
- `matpower_data::Dict`: Dados da rede no formato MATPOWER
- `n_regions::Int`: Número de regiões desejadas
- `method::Symbol`: Método de decomposição (:metis, :geographical, :manual)

# Retorna
- `regions::Vector{AladinRegion}`: Vetor de regiões decompostas
- `consensus_matrix::Matrix{Float64}`: Matriz de consenso global
- `consensus_rhs::Vector{Float64}`: Lado direito das restrições de consenso
"""
function decompose_network(
    matpower_data::Dict,
    n_regions::Int;
    method::Symbol = :geographical
)
    if method == :geographical
        return decompose_geographical(matpower_data, n_regions)
    elseif method == :manual
        return decompose_manual(matpower_data, n_regions)
    else
        error("Método de decomposição '$method' não implementado")
    end
end

"""
    decompose_geographical(matpower_data::Dict, n_regions::Int)

Decompõe a rede baseado em proximidade geográfica ou topológica.
"""
function decompose_geographical(matpower_data::Dict, n_regions::Int)
    bus_params = matpower_data["bus_params"]
    branch_params = matpower_data["branch_params"]
    gen_params = matpower_data["generator_params"]
    
    n_buses = length(bus_params.bus_ids)
    
    # Atribuir buses às regiões de forma balanceada
    buses_per_region = div(n_buses, n_regions)
    region_assignments = Dict{Int, Int}()
    
    for (idx, bus_id) in enumerate(bus_params.bus_ids)
        region_id = min(div(idx - 1, buses_per_region) + 1, n_regions)
        region_assignments[bus_id] = region_id
    end
    
    # Criar regiões
    regions = AladinRegion[]
    
    for r in 1:n_regions
        # Identificar buses desta região
        region_buses = [bus_id for (bus_id, reg) in region_assignments if reg == r]
        
        # Criar dados da região
        region_data = extract_region_data(
            matpower_data,
            region_buses,
            region_assignments
        )
        
        # Criar objeto de região
        region = AladinRegion(r, region_data)
        push!(regions, region)
    end
    
    # Construir matriz de consenso para variáveis de fronteira
    consensus_matrix, consensus_rhs = build_consensus_constraints(
        regions,
        matpower_data,
        region_assignments
    )
    
    return regions, consensus_matrix, consensus_rhs
end

"""
    extract_region_data(matpower_data::Dict, region_buses::Vector{Int}, 
                       region_assignments::Dict{Int,Int})

Extrai os dados de uma região específica.
"""
function extract_region_data(
    matpower_data::Dict,
    region_buses::Vector{Int},
    region_assignments::Dict{Int, Int}
)
    region_data = deepcopy(matpower_data)
    
    # Filtrar buses
    bus_params = region_data["bus_params"]
    bus_mask = [bus_id in region_buses for bus_id in bus_params.bus_ids]
    
    # Atualizar parâmetros de bus
    region_data["bus_params"] = filter_bus_params(bus_params, region_buses)
    
    # Filtrar geradores que pertencem aos buses da região
    gen_params = region_data["generator_params"]
    region_data["generator_params"] = filter_generator_params(
        gen_params,
        region_buses
    )
    
    # Filtrar branches (incluir branches de fronteira)
    branch_params = region_data["branch_params"]
    region_data["branch_params"] = filter_branch_params(
        branch_params,
        region_buses,
        region_assignments
    )
    
    # Filtrar outros componentes se existirem
    if haskey(region_data, "battery_params")
        region_data["battery_params"] = filter_battery_params(
            region_data["battery_params"],
            region_buses
        )
    end
    
    if haskey(region_data, "prosumer_params")
        region_data["prosumer_params"] = filter_prosumer_params(
            region_data["prosumer_params"],
            region_buses
        )
    end
    
    return region_data
end

"""
    filter_bus_params(params, region_buses::Vector{Int})

Filtra parâmetros de barramento para incluir apenas buses da região.
"""
function filter_bus_params(params, region_buses::Vector{Int})
    # Criar novo objeto com apenas os buses da região
    filtered = deepcopy(params)
    
    # Filtrar cada campo
    bus_mask = [bus_id in region_buses for bus_id in params.bus_ids]
    
    filtered.bus_ids = params.bus_ids[bus_mask]
    filtered.vmin = Dict(k => v for (k, v) in params.vmin if k in region_buses)
    filtered.vmax = Dict(k => v for (k, v) in params.vmax if k in region_buses)
    filtered.vstart = Dict(k => v for (k, v) in params.vstart if k in region_buses)
    
    if hasfield(typeof(params), :pd)
        filtered.pd = Dict(k => v for (k, v) in params.pd if k in region_buses)
    end
    
    if hasfield(typeof(params), :qd)
        filtered.qd = Dict(k => v for (k, v) in params.qd if k in region_buses)
    end
    
    return filtered
end

"""
    filter_generator_params(params, region_buses::Vector{Int})

Filtra parâmetros de geradores para incluir apenas os conectados aos buses da região.
"""
function filter_generator_params(params, region_buses::Vector{Int})
    filtered = deepcopy(params)
    
    # Identificar geradores conectados aos buses da região
    gen_mask = [params.bus[g] in region_buses for g in 1:params.number_of_generators]
    
    # Filtrar arrays e dicionários
    filtered.number_of_generators = sum(gen_mask)
    
    if hasfield(typeof(params), :bus)
        filtered.bus = params.bus[gen_mask]
    end
    
    if hasfield(typeof(params), :pmin)
        filtered.pmin = params.pmin[gen_mask]
    end
    
    if hasfield(typeof(params), :pmax)
        filtered.pmax = params.pmax[gen_mask]
    end
    
    if hasfield(typeof(params), :qmin)
        filtered.qmin = params.qmin[gen_mask]
    end
    
    if hasfield(typeof(params), :qmax)
        filtered.qmax = params.qmax[gen_mask]
    end
    
    if hasfield(typeof(params), :generator_status)
        filtered.generator_status = params.generator_status[gen_mask]
    end
    
    return filtered
end

"""
    filter_branch_params(params, region_buses::Vector{Int}, 
                        region_assignments::Dict{Int,Int})

Filtra parâmetros de linhas. Inclui linhas internas e de fronteira.
"""
function filter_branch_params(
    params,
    region_buses::Vector{Int},
    region_assignments::Dict{Int, Int}
)
    filtered = deepcopy(params)
    
    # Identificar branches que conectam buses da região
    # Incluir branches de fronteira (conectam a outras regiões)
    branch_mask = Bool[]
    
    for i in 1:params.number_of_branches
        from_bus = params.bus_from[i]
        to_bus = params.bus_to[i]
        
        # Incluir se pelo menos um dos buses está na região
        include_branch = (from_bus in region_buses) || (to_bus in region_buses)
        push!(branch_mask, include_branch)
    end
    
    # Filtrar
    filtered.number_of_branches = sum(branch_mask)
    
    if hasfield(typeof(params), :bus_from)
        filtered.bus_from = params.bus_from[branch_mask]
    end
    
    if hasfield(typeof(params), :bus_to)
        filtered.bus_to = params.bus_to[branch_mask]
    end
    
    if hasfield(typeof(params), :capacity)
        filtered.capacity = params.capacity[branch_mask]
    end
    
    if hasfield(typeof(params), :resistance)
        filtered.resistance = params.resistance[branch_mask]
    end
    
    if hasfield(typeof(params), :reactance)
        filtered.reactance = params.reactance[branch_mask]
    end
    
    if hasfield(typeof(params), :b)
        filtered.b = params.b[branch_mask]
    end
    
    if hasfield(typeof(params), :existing)
        filtered.existing = params.existing[branch_mask]
    end
    
    return filtered
end

"""
    filter_battery_params(params, region_buses::Vector{Int})

Filtra parâmetros de baterias.
"""
function filter_battery_params(params, region_buses::Vector{Int})
    # Implementação similar aos geradores
    return deepcopy(params)
end

"""
    filter_prosumer_params(params, region_buses::Vector{Int})

Filtra parâmetros de prosumers.
"""
function filter_prosumer_params(params, region_buses::Vector{Int})
    # Implementação similar aos geradores
    return deepcopy(params)
end

"""
    build_consensus_constraints(regions::Vector{AladinRegion}, 
                                matpower_data::Dict,
                                region_assignments::Dict{Int,Int})

Constrói as restrições de consenso para variáveis de fronteira.

As restrições de consenso garantem que:
1. Fluxos de potência em branches de fronteira sejam consistentes
2. Ângulos de tensão em buses de fronteira sejam iguais
3. Magnitudes de tensão em buses de fronteira sejam iguais
"""


# =================================================================================
# IMPLEMENTAÇÃO COMPLETA DA FUNÇÃO build_consensus_constraints
# Foco: Consenso de Tensão (Ângulo e Magnitude) em Barras de Fronteira
# =================================================================================

using LinearAlgebra

# Estruturas de dados simuladas para fins de demonstração
# (Assumindo que estas estruturas existem e são preenchidas corretamente)
mutable struct AladinRegion
    id::Int
    data::Dict
    coupling_vars::Vector{Symbol}
    A_consensus::Union{Nothing, Matrix{Float64}}
    # ... outros campos
end

# Função auxiliar para obter o índice local de uma barra no vetor de variáveis
# A ordem das variáveis é: [theta_1, ..., theta_N, V_1, ..., V_N, ...]
function get_local_bus_var_index(region::AladinRegion, bus_id::Int, var_type::Symbol)
    bus_ids = region.data["bus_params"].bus_ids
    local_idx = findfirst(==(bus_id), bus_ids)
    
    if isnothing(local_idx)
        return nothing
    end
    
    n_buses = length(bus_ids)
    T = get(region.data, "time_periods", 1)
    
    # Assumindo que as variáveis são agrupadas por tipo e depois por tempo
    # Ex: [theta_t1_b1, ..., theta_t1_bN, theta_t2_b1, ...]
    
    # Para simplificar, vamos assumir T=1 e a ordem: [theta_b1, ..., theta_bN, V_b1, ..., V_bN]
    
    if var_type == :theta
        # theta está nas primeiras N posições
        return local_idx
    elseif var_type == :v_magnitude || var_type == :v_magnitude_sq
        # V está nas posições N+1 a 2N
        return n_buses * T + local_idx
    else
        error("Tipo de variável desconhecido para consenso de barra: $var_type")
    end
end

# Função auxiliar para obter o tamanho total do vetor de variáveis de uma região
function get_region_var_size(region::AladinRegion)
    n_buses = length(region.data["bus_params"].bus_ids)
    n_gens = region.data["generator_params"].number_of_generators
    T = get(region.data, "time_periods", 1)
    
    # Assumindo AC/SOCP: theta (N*T) + V (N*T) + Pgen (G*T) + Qgen (G*T)
    return 2 * n_buses * T + 2 * n_gens * T
end

"""
    build_consensus_constraints(regions::Vector{AladinRegion}, 
                                matpower_data::Dict,
                                region_assignments::Dict{Int,Int})

Constrói as restrições de consenso para variáveis de fronteira (tensão).

As restrições de consenso garantem que:
1. Ângulos de tensão em barras de fronteira sejam iguais: theta_r1 - theta_r2 = 0
2. Magnitudes de tensão em barras de fronteira sejam iguais: V_r1 - V_r2 = 0

# Argumentos
- `regions::Vector{AladinRegion}`: Vetor de regiões decompostas
- `matpower_data::Dict`: Dados completos da rede
- `region_assignments::Dict{Int,Int}`: Mapeamento de barras para regiões

# Retorna
- `consensus_matrix::Matrix{Float64}`: Matriz de consenso global A
- `consensus_rhs::Vector{Float64}`: Lado direito das restrições de consenso b
"""
function build_consensus_constraints(
    regions::Vector{AladinRegion},
    matpower_data::Dict,
    region_assignments::Dict{Int, Int}
)
    branch_params = matpower_data["branch_params"]
    
    # 1. Identificar Barras de Fronteira Duplicadas
    # Uma barra é de fronteira se estiver conectada a um ramo que liga duas regiões.
    # No ALADIN, as barras de fronteira são duplicadas e pertencem a ambas as regiões.
    
    # Mapeamento: (bus_id) -> Set of region_ids
    bus_to_regions = Dict{Int, Set{Int}}()
    
    for r in regions
        for bus_id in r.data["bus_params"].bus_ids
            if !haskey(bus_to_regions, bus_id)
                bus_to_regions[bus_id] = Set{Int}()
            end
            push!(bus_to_regions[bus_id], r.id)
        end
    end
    
    # Identificar pares de regiões que compartilham uma barra de fronteira
    # Mapeamento: (bus_id) -> Vector of (region_id, local_index)
    boundary_bus_pairs = Dict{Int, Vector{Tuple{Int, Int}}}()
    
    for bus_id in keys(bus_to_regions)
        regions_sharing = collect(bus_to_regions[bus_id])
        
        if length(regions_sharing) > 1
            # Esta é uma barra de fronteira (duplicada)
            
            # Para cada par de regiões que compartilha a barra, criamos restrições
            # (r1, r2) onde r1 < r2
            for i in 1:length(regions_sharing)
                for j in i+1:length(regions_sharing)
                    r1_id = regions_sharing[i]
                    r2_id = regions_sharing[j]
                    
                    # Encontrar os objetos AladinRegion
                    r1 = regions[r1_id]
                    r2 = regions[r2_id]
                    
                    # Encontrar o índice local da barra nas variáveis de cada região
                    local_idx_r1 = get_local_bus_var_index(r1, bus_id, :theta) # Usando theta como proxy
                    local_idx_r2 = get_local_bus_var_index(r2, bus_id, :theta)
                    
                    if !isnothing(local_idx_r1) && !isnothing(local_idx_r2)
                        if !haskey(boundary_bus_pairs, bus_id)
                            boundary_bus_pairs[bus_id] = []
                        end
                        # Armazenar o par de regiões e o índice local da variável theta
                        push!(boundary_bus_pairs[bus_id], (r1_id, local_idx_r1))
                        push!(boundary_bus_pairs[bus_id], (r2_id, local_idx_r2))
                    end
                end
            end
        end
    end
    
    # 2. Calcular Dimensões
    
    # O número de restrições de consenso é 2 (theta e V) * número de pares de regiões
    # que compartilham uma barra de fronteira.
    n_consensus = 0
    for bus_id in keys(boundary_bus_pairs)
        # Cada barra de fronteira compartilhada por N regiões gera N-1 restrições de consenso
        # para cada variável (theta e V).
        # Simplificando: 2 restrições (theta, V) para cada par de regiões (r1, r2)
        n_consensus += 2 * (length(bus_to_regions[bus_id]) - 1)
    end
    
    # Calcular o número total de variáveis em todas as regiões (tamanho do vetor x_all)
    n_total_vars = sum(get_region_var_size(r) for r in regions)
    
    # 3. Inicializar Matrizes
    consensus_matrix = zeros(n_consensus, n_total_vars)
    consensus_rhs = zeros(n_consensus)
    
    # 4. Preencher a Matriz de Consenso
    
    # Mapeamento de índice global
    var_offsets = [0]
    for r in regions
        push!(var_offsets, var_offsets[end] + get_region_var_size(r))
    end
    
    constraint_idx = 1
    
    for bus_id in keys(bus_to_regions)
        regions_sharing = collect(bus_to_regions[bus_id])
        
        if length(regions_sharing) > 1
            # Criar restrições de consenso para theta e V
            
            # Iterar sobre os pares de regiões (r1, r2) para criar as restrições
            for i in 1:length(regions_sharing)
                r1_id = regions_sharing[i]
                r1 = regions[r1_id]
                
                # Criar restrições apenas com o primeiro par (r1, r_base)
                if i == 1
                    continue # r1 será a região base
                end
                
                r_base_id = regions_sharing[1]
                r_base = regions[r_base_id]
                
                # --- Consenso de Ângulo (theta_r1 - theta_r_base = 0) ---
                
                # Posição da variável theta na região r1
                local_idx_theta_r1 = get_local_bus_var_index(r1, bus_id, :theta)
                global_idx_theta_r1 = var_offsets[r1_id] + local_idx_theta_r1
                
                # Posição da variável theta na região r_base
                local_idx_theta_r_base = get_local_bus_var_index(r_base, bus_id, :theta)
                global_idx_theta_r_base = var_offsets[r_base_id] + local_idx_theta_r_base
                
                # Restrição: 1 * theta_r1 - 1 * theta_r_base = 0
                consensus_matrix[constraint_idx, global_idx_theta_r1] = 1.0
                consensus_matrix[constraint_idx, global_idx_theta_r_base] = -1.0
                consensus_rhs[constraint_idx] = 0.0
                constraint_idx += 1
                
                # --- Consenso de Magnitude (V_r1 - V_r_base = 0) ---
                
                # Posição da variável V na região r1
                local_idx_v_r1 = get_local_bus_var_index(r1, bus_id, :v_magnitude)
                global_idx_v_r1 = var_offsets[r1_id] + local_idx_v_r1
                
                # Posição da variável V na região r_base
                local_idx_v_r_base = get_local_bus_var_index(r_base, bus_id, :v_magnitude)
                global_idx_v_r_base = var_offsets[r_base_id] + local_idx_v_r_base
                
                # Restrição: 1 * V_r1 - 1 * V_r_base = 0
                consensus_matrix[constraint_idx, global_idx_v_r1] = 1.0
                consensus_matrix[constraint_idx, global_idx_v_r_base] = -1.0
                consensus_rhs[constraint_idx] = 0.0
                constraint_idx += 1
            end
        end
    end
    
    # 5. Armazenar a matriz de consenso local em cada região (A_consensus)
    for r_id in 1:length(regions)
        r = regions[r_id]
        offset_start = var_offsets[r_id] + 1
        offset_end = var_offsets[r_id+1]
        
        # Extrair as colunas da matriz global que correspondem às variáveis desta região
        r.A_consensus = consensus_matrix[:, offset_start:offset_end]
    end
    
    return consensus_matrix, consensus_rhs
end


# 1. Esta implementação assume que as barras de fronteira são duplicadas em todas
#    as regiões que compartilham o ramo de acoplamento.
# 2. A ordem das variáveis no vetor local (y) é crucial e deve ser a mesma
#    definida em `extract_primal_solution` e `get_primal_variables_vector`.
# 3. A lógica de consenso de tensão (theta e V) é a mais comum para decomposição
#    por barras. O consenso de fluxo de potência (P_flow_r1 + P_flow_r2 = 0)
#    também é possível, mas requer um mapeamento mais complexo das variáveis de fluxo.
# 4. As funções auxiliares `get_local_bus_var_index` e `get_region_var_size`
#    são simplificadas e precisam ser ajustadas para a estrutura exata de dados
#    e o número de períodos de tempo (T) do seu projeto.


"""
    decompose_manual(matpower_data::Dict, n_regions::Int)

Permite decomposição manual especificando quais buses pertencem a cada região.
"""
function decompose_manual(matpower_data::Dict, n_regions::Int)
    # Placeholder para decomposição manual
    # Usuário pode especificar region_assignments manualmente
    error("Decomposição manual não implementada. Use decompose_geographical.")
end

"""
    identify_boundary_variables(regions::Vector{AladinRegion})

Identifica variáveis que aparecem em múltiplas regiões (variáveis de fronteira).
"""
function identify_boundary_variables(regions::Vector{AladinRegion})
    boundary_vars = Dict{Symbol, Vector{Int}}()
    
    # Para cada tipo de variável (theta, v_magnitude, p_flow, etc.)
    # identificar quais aparecem em múltiplas regiões
    
    # TODO: Implementar lógica de identificação
    
    return boundary_vars
end

"""
    create_consensus_matrix_for_power_flow(regions, boundary_branches)

Cria matriz de consenso específica para fluxo de potência.

Para cada branch de fronteira (i,j) conectando região r1 e r2:
- Variável p_flow em r1 deve ser negativa da variável p_flow em r2
- Similar para q_flow
"""
function create_consensus_matrix_for_power_flow(
    regions::Vector{AladinRegion},
    boundary_branches::Vector{Int}
)
    # TODO: Implementar construção detalhada da matriz
    
    n_consensus = length(boundary_branches) * 2
    n_vars = sum(length(r.data["bus_params"].bus_ids) for r in regions)
    
    A = zeros(n_consensus, n_vars)
    b = zeros(n_consensus)
    
    return A, b
end

"""
    visualize_decomposition(regions::Vector{AladinRegion})

Visualiza a decomposição da rede (opcional, requer pacotes de plotagem).
"""
function visualize_decomposition(regions::Vector{AladinRegion})
    println("\nDecomposição da rede:")
    println("="^50)
    
    for region in regions
        n_buses = length(region.data["bus_params"].bus_ids)
        n_gens = region.data["generator_params"].number_of_generators
        n_branches = region.data["branch_params"].number_of_branches
        
        println("Região $(region.id):")
        println("  - Buses: $n_buses")
        println("  - Geradores: $n_gens")
        println("  - Linhas: $n_branches")
    end
    
    println("="^50)
end
