function storages!(
    model::JuMP.Model,
    sys::SystemData,
    model_type::Type{<:OPFModel},
)
    stor = sys.storage

    if stor.n == 0
        @info "No storage information in SystemData. Skipping storage modeling."
        return nothing
    end

    add_storage_variables!(model, stor)
    add_storage_constraints!(model, stor, model_type)

    return nothing
end

# =======================================================================
# VARIABLES
# =======================================================================

function add_storage_variables!(model::JuMP.Model, stor::StorageData)
    add_storage_energy_variables!(model, stor)
    add_charge_discharge_variables!(model, stor)
    add_storage_power_variables!(model, stor)
    return nothing
end

# =======================================================================
# CONSTRAINTS
# =======================================================================

function add_storage_constraints!(model::JuMP.Model, stor::StorageData, ::Type{DC})
    return add_dc_storage_constraints!(model, stor)
end

function add_storage_constraints!(model::JuMP.Model, stor::StorageData, ::Type{AC})
    return add_ac_storage_constraints!(model, stor)
end

function add_storage_constraints!(model::JuMP.Model, stor::StorageData, ::Type{ACR})
    return add_acr_storage_constraints!(model, stor)
end

function add_storage_constraints!(model::JuMP.Model, stor::StorageData, ::Type{IVR})
    return add_acr_storage_constraints!(model, stor)
end

function add_storage_constraints!(model::JuMP.Model, stor::StorageData, ::Type{LPAC})
    return add_dc_storage_constraints!(model, stor)
end

function add_storage_constraints!(model::JuMP.Model, stor::StorageData, ::Type{SOCWR})
    return add_socwr_storage_constraints!(model, stor)
end

function add_storage_constraints!(model::JuMP.Model, stor::StorageData, ::Type{SDPWRM})
    return add_socwr_storage_constraints!(model, stor)
end

function add_storage_constraints!(model::JuMP.Model, stor::StorageData, ::Type{QC})
    return add_socwr_storage_constraints!(model, stor)
end

function add_storage_constraints!(model::JuMP.Model, stor::StorageData, ::Type{SOCBF})
    return add_socwr_storage_constraints!(model, stor)
end

function add_storage_constraints!(model::JuMP.Model, stor::StorageData, ::Type{BFA})
    return add_dc_storage_constraints!(model, stor)
end

function add_storage_constraints!(model::JuMP.Model, stor::StorageData, ::Type{NFA})
    return add_dc_storage_constraints!(model, stor)
end
