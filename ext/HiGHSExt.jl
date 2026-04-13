module HiGHSExt

using Lirion
using HiGHS
using JuMP: optimizer_with_attributes

const _HIGHS_SILENT_FACTORY = optimizer_with_attributes(
    HiGHS.Optimizer,
    "log_to_console" => false,
    "output_flag"    => false,
)

Lirion._silent_factory(::Type{HiGHS.Optimizer}) = _HIGHS_SILENT_FACTORY

end # module
