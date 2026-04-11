module GurobiExt

using Lirion
using Gurobi
using JuMP: optimizer_with_attributes

const _GUROBI_SILENT_FACTORY = optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0)

Lirion._silent_factory(::Type{Gurobi.Optimizer})  = _GUROBI_SILENT_FACTORY
Lirion._gurobi_socp_factory()                     = _GUROBI_SILENT_FACTORY
Lirion._default_socp_optimizer()                  = Gurobi.Optimizer

end # module
