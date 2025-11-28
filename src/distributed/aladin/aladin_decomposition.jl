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
function build_consensus_constraints(
    regions::Vector{AladinRegion},
    matpower_data::Dict,
    region_assignments::Dict{Int, Int}
)
    branch_params = matpower_data["branch_params"]
    
    # Identificar branches de fronteira (conectam regiões diferentes)
    boundary_branches = Int[]
    
    for i in 1:branch_params.number_of_branches
        from_bus = branch_params.bus_from[i]
        to_bus = branch_params.bus_to[i]
        
        from_region = region_assignments[from_bus]
        to_region = region_assignments[to_bus]
        
        if from_region != to_region
            push!(boundary_branches, i)
        end
    end
    
    # Para cada branch de fronteira, criar restrições de consenso
    # Exemplo: P_flow(from) + P_flow(to) = 0 (conservação de potência)
    
    n_consensus = length(boundary_branches) * 2  # P e Q para cada branch
    n_total_vars = sum(length(r.data["bus_params"].bus_ids) for r in regions)
    
    # Matriz de consenso (simplificada por enquanto)
    consensus_matrix = zeros(n_consensus, n_total_vars)
    consensus_rhs = zeros(n_consensus)
    
    # TODO: Preencher matriz de consenso baseado nas variáveis de fronteira
    
    return consensus_matrix, consensus_rhs
end

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
