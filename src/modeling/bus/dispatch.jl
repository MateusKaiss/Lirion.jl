
# =======================================================================
# Dispatch by model type (DC / AC / SOCP)
# =======================================================================

function buses!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{DC},
)
    bus = sys.bus
    create_voltage_angle_variables!(model, bus)
    return nothing
end

function buses!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{AC},
)
    bus = sys.bus
    create_voltage_angle_variables!(model, bus)
    create_voltage_magnitude_variables!(model, bus)
    return nothing
end

function buses!(
    model::JuMP.Model,
    sys::SystemData,
    ::Type{SOCP},
)
    bus = sys.bus
    create_voltage_angle_variables!(model, bus)
    create_socp_bus_variables!(model, bus)
    return nothing
end
