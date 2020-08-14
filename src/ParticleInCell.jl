module ParticleInCell

export scatter_charge!, shape_1st_order, shape_2nd_order, shape_3rd_order, shape_4th_order

include("grid.jl")
include("field.jl")
include("species.jl")

abstract type IntegrationStep end

mutable struct Simulation{FT, IT}
    num_species::IT
    num_fields::IT

    species_dict::Dict{String, IT}
    fields_dict::Dict{String, IT}

    species::Vector{Species}
    fields::Vector{Field}

    timestep::FT
    time::FT
    integration_steps::Vector{IntegrationStep}

    function Simulation(timestep::FT) where FT
        num_species = 0
        num_fields = 0

        species_dict = Dict{String, Int64}()
        fields_dict = Dict{String, Int64}()

        species = Vector{Species}()
        fields = Vector{Field}()

        time = zero(FT)
        integration_steps = Vector{IntegrationStep}()

        return new{FT, Int64}(num_species, num_fields, species_dict, fields_dict, species,
                   fields, timestep, time, integration_steps)
    end
end

function add_species!(sim::Simulation{FT, IT}, species::Species{FT, IT},
                      label::String) where {FT, IT}
    push!(sim.species, species)
    sim.num_species += 1
    sim.species_dict[label] = sim.num_species
end

function add_field!(sim::Simulation{FT, IT}, field::Field{FT},
                      label::String) where {FT, IT}
    push!(sim.fields, field)
    sim.num_fields += 1
    sim.fields_dict[label] = sim.num_fields
end

function add_integration_step!(sim::Simulation, step::IntegrationStep)
    push!(sim.integration_steps, step)
end

function setup!(sim::Simulation)
    for step in sim.integration_steps
        setup!(step, sim)
    end
end
setup!(::IntegrationStep, ::Simulation) = nothing

function step!(sim::Simulation)
    for step in sim.integration_steps
        step!(step, sim)
    end
end
step!(::IntegrationStep, ::Simulation) = throw("No step method")

include("shape_functions.jl")
include("scatter.jl")
include("gather.jl")
include("particle_push.jl")
include("field_solve.jl")

end # module
