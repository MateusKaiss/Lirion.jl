module GLPKExt

using Lirion
using GLPK
using JuMP: optimizer_with_attributes

const _GLPK_SILENT_FACTORY = optimizer_with_attributes(
    GLPK.Optimizer,
    "msg_lev" => 0,
)

Lirion._silent_factory(::Type{GLPK.Optimizer}) = _GLPK_SILENT_FACTORY

end # module
