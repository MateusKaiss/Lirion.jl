
@inline formulation_of(::DCStateWorkspace) = DC
@inline formulation_of(::DCPowerWorkspace) = DC
@inline formulation_of(::DCHybridWorkspace) = DC
@inline formulation_of(::ACStateWorkspace) = AC
@inline formulation_of(::ACPowerWorkspace) = AC
@inline formulation_of(::ACHybridWorkspace) = AC
@inline formulation_of(::SOCWRStateWorkspace) = SOCWR
@inline formulation_of(::SOCWRPowerWorkspace) = SOCWR
@inline formulation_of(::SOCWRHybridWorkspace) = SOCWR

@inline coupling_of(::DCStateWorkspace) = State
@inline coupling_of(::DCPowerWorkspace) = Power
@inline coupling_of(::DCHybridWorkspace) = Hybrid
@inline coupling_of(::ACStateWorkspace) = State
@inline coupling_of(::ACPowerWorkspace) = Power
@inline coupling_of(::ACHybridWorkspace) = Hybrid
@inline coupling_of(::SOCWRStateWorkspace) = State
@inline coupling_of(::SOCWRPowerWorkspace) = Power
@inline coupling_of(::SOCWRHybridWorkspace) = Hybrid

@inline formulation_of(::DCStateWorkspaceALADIN) = DC
@inline formulation_of(::DCPowerWorkspaceALADIN) = DC
@inline formulation_of(::DCHybridWorkspaceALADIN) = DC
@inline formulation_of(::ACStateWorkspaceALADIN) = AC
@inline formulation_of(::ACPowerWorkspaceALADIN) = AC
@inline formulation_of(::ACHybridWorkspaceALADIN) = AC
@inline formulation_of(::SOCWRStateWorkspaceALADIN) = SOCWR
@inline formulation_of(::SOCWRPowerWorkspaceALADIN) = SOCWR
@inline formulation_of(::SOCWRHybridWorkspaceALADIN) = SOCWR

@inline coupling_of(::DCStateWorkspaceALADIN) = State
@inline coupling_of(::DCPowerWorkspaceALADIN) = Power
@inline coupling_of(::DCHybridWorkspaceALADIN) = Hybrid
@inline coupling_of(::ACStateWorkspaceALADIN) = State
@inline coupling_of(::ACPowerWorkspaceALADIN) = Power
@inline coupling_of(::ACHybridWorkspaceALADIN) = Hybrid
@inline coupling_of(::SOCWRStateWorkspaceALADIN) = State
@inline coupling_of(::SOCWRPowerWorkspaceALADIN) = Power
@inline coupling_of(::SOCWRHybridWorkspaceALADIN) = Hybrid

# ----------------------------------------------------------------------
# Algorithm family discriminators
# ----------------------------------------------------------------------

@inline is_admm(::AbstractADMMWorkspace) = true
@inline is_admm(::AbstractALADINWorkspace) = false

@inline is_aladin(::AbstractADMMWorkspace) = false
@inline is_aladin(::AbstractALADINWorkspace) = true

@inline function tag_of(ws::AbstractADMMWorkspace)
    return ADMM{coupling_of(ws), formulation_of(ws)}()
end

@inline function tag_of(ws::AbstractALADINWorkspace)
    return ALADIN{coupling_of(ws), formulation_of(ws)}()
end

# ----------------------------------------------------------------------
# Workspace → consensus variables (sugar over traits.jl)
# ----------------------------------------------------------------------

@inline consensus_vars(ws::AbstractADMMWorkspace) =
    consensus_vars(formulation_of(ws), coupling_of(ws))

@inline consensus_vars(ws::AbstractALADINWorkspace) =
    consensus_vars(formulation_of(ws), coupling_of(ws))
