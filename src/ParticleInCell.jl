module ParticleInCell

using FFTW

export UniformGrid, Field, Species, total_cells, cell_length
export shape_1st_order, shape_2nd_order, shape_3rd_order
export ScatterChargeToGrid, GatherForcesFromGrid, SymplecticEulerPush
export ConstrainSpecies, FourierFieldSolve, FiniteDifferenceToNodes
export PrescribeField
export Simulation, add_species!, add_field!, add_integration_step!, setup!, step!

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
    num_timesteps::IT
    time::FT
    integration_steps::Vector{IntegrationStep}

    total_time::FT
    integration_times::Vector{FT}

    setup_complete::Bool

    function Simulation(timestep::FT) where FT
        num_species = 0
        num_fields = 0

        species_dict = Dict{String, Int64}()
        fields_dict = Dict{String, Int64}()

        species = Vector{Species}()
        fields = Vector{Field}()

        time = zero(FT)
        integration_steps = Vector{IntegrationStep}()

        integration_times = Vector{FT}()

        return new{FT, Int64}(num_species, num_fields, species_dict, fields_dict, species,
                              fields, timestep, 0, time, integration_steps, zero(FT), integration_times, false)
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

function add_integration_step!(sim::Simulation{FT}, step::IntegrationStep) where FT
    push!(sim.integration_steps, step)
    push!(sim.integration_times, zero(FT))
end

function setup!(sim::Simulation)
    for step in sim.integration_steps
        setup!(step, sim)
    end
    sim.setup_complete = true
end
setup!(::IntegrationStep, ::Simulation) = nothing

function step!(sim::Simulation)
    start_time = time()
    for (i, step) in enumerate(sim.integration_steps)
        step_start_time = time()
        step!(step, sim)
        sim.integration_times[i] += time() - step_start_time
    end
    sim.time += sim.timestep
    sim.num_timesteps += 1
    sim.total_time += time() - start_time
end
step!(::IntegrationStep, ::Simulation) = throw("No step method")

function print_summary(io::IO, sim::Simulation{FT, IT}) where {FT, IT}
    println(io, "Simulation{$FT, $IT}")
    println(io, "Simulation setup ", sim.setup_complete ? "has" : "has not", " been completed.")
    println(io, "Simulated $(sim.time) seconds using $(sim.num_steps) timesteps".)
#    println(io, "Total simulation wall time is $(sim.total_time) seconds.")
#    println(io)
#    println(io, "  $num_species species")
#    for s in keys(sim.species_dict)
#        println(io, "    ", s)
#    end
#    println(io)
#    println(io, "  $num_fields fields")
#    for f in keys(sim.fields_dict)
#        println(io, "    ", f)
#    end
#    println(io)
#    println(io, "  $(length(integration_steps)) integration steps")
#    for (i, s) in enumerate(sim.integration_steps)
#        println(io, "    $s")
#        println(io, "      $(sim.integration_times[i]) seconds ($(sim.integration_times[i] / sim.total_time * 100)%)")
#    end
end

include("shape_functions.jl")
include("scatter.jl")
include("gather.jl")
include("particle_push.jl")
include("field_solve.jl")
include("prescribe_field.jl")

end # module
