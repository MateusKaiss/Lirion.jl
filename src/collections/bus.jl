struct BusData
    n::Int
    id::Vector{Int}
    bus_type::Vector{Int}
    pd::Vector{Float64}
    qd::Vector{Float64}
    gs::Vector{Float64}
    bs::Vector{Float64}
    area::Vector{Int}
    vm::Vector{Float64}
    va::Vector{Float64}
    base_kv::Vector{Float64}
    zone::Vector{Int}
    vmax::Vector{Float64}
    vmin::Vector{Float64}
    
    reference_bus::Union{Nothing,Int}
    vmin_map::Dict{Int,Float64}
    vmax_map::Dict{Int,Float64}
    vstart_map::Dict{Int,Union{Nothing,Float64}}
    vastart_map::Dict{Int,Union{Nothing,Float64}}
end

# ---------------------------------------------------------------------
# Getting Bus Parameters
# ---------------------------------------------------------------------
function get_bus_parameters(bus_data)
    bus_ids = bus_data["bus_i"]

    reference_bus_idx = findfirst(==(3), bus_data["bus_type"])
    reference_bus = nothing
    if reference_bus_idx !== nothing
        reference_bus = bus_ids[reference_bus_idx]
    else
        @info "No reference bus (type 3) found in this zone. Skipping θ_ref constraint."
    end

    vmin    = Dict{Int, Float64}()
    vmax    = Dict{Int, Float64}()
    vstart  = Dict{Int, Union{Nothing, Float64}}()
    vastart = Dict{Int, Union{Nothing, Float64}}()

    @inbounds for i in eachindex(bus_ids)
        b = bus_ids[i]
        vmin[b]    = bus_data["vmin"][i]
        vmax[b]    = bus_data["vmax"][i]
        vstart[b]  = isnan(bus_data["vm"][i]) ? nothing : bus_data["vm"][i]
        vastart[b] = isnan(bus_data["va"][i]) ? nothing : bus_data["va"][i]
    end

    return (
        bus_ids = bus_ids,
        reference_bus = reference_bus,
        vmin = vmin,
        vmax = vmax,
        vstart = vstart,
        initial_voltage_angle = vastart,
    )
end

function build_bus_data(data::Dict{String,Any})::BusData
    bus = data["bus"]
    n   = length(bus["bus_i"])
    params = get_bus_parameters(bus)

    return BusData(
        n,
        bus["bus_i"],
        bus["bus_type"],
        bus["pd"],
        bus["qd"],
        bus["gs"],
        bus["bs"],
        bus["area"],
        bus["vm"],
        bus["va"],
        bus["base_kv"],
        bus["zone"],
        bus["vmax"],
        bus["vmin"],

        params.reference_bus,
        params.vmin,
        params.vmax,
        params.vstart,
        params.initial_voltage_angle,
    )
end