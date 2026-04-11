module SCSExt

using Lirion
using SCS
using JuMP: optimizer_with_attributes

const _SCS_SILENT_FACTORY = optimizer_with_attributes(SCS.Optimizer, "verbose" => 0)

Lirion._silent_factory(::Type{SCS.Optimizer}) = _SCS_SILENT_FACTORY
Lirion._scs_socp_factory()                   = _SCS_SILENT_FACTORY
Lirion._default_socp_optimizer()             = SCS.Optimizer

end # module
