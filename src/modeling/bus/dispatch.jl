function buses!(
    model::JuMP.Model,
    sys::SystemData,
    model_type::Type{<:OPFModel},
)
    bus = sys.bus

    add_bus_variables!(model, bus, model_type)

    return nothing
end

# =======================================================================
# VARIABLES
# =======================================================================

function add_bus_variables!(model::JuMP.Model, bus::BusData, ::Type{DC})
    return add_voltage_angle_variables!(model, bus)
end

function add_bus_variables!(model::JuMP.Model, bus::BusData, ::Type{AC})
    add_voltage_angle_variables!(model, bus)
    return add_voltage_magnitude_variables!(model, bus)
end

function add_bus_variables!(model::JuMP.Model, bus::BusData, ::Type{ACR})
    return add_rectangular_bus_variables!(model, bus)
end

function add_bus_variables!(model::JuMP.Model, bus::BusData, ::Type{IVR})
    return add_rectangular_bus_variables!(model, bus)
end

function add_bus_variables!(model::JuMP.Model, bus::BusData, ::Type{LPAC})
    add_voltage_angle_variables!(model, bus)
    return add_lpac_bus_variables!(model, bus)
end

function add_bus_variables!(model::JuMP.Model, bus::BusData, ::Type{SOCWR})
    add_voltage_angle_variables!(model, bus)
    return add_socwr_bus_variables!(model, bus)
end

function add_bus_variables!(model::JuMP.Model, bus::BusData, ::Type{SOCBF})
    return add_socwr_bus_variables!(model, bus)
end
