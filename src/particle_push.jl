mutable struct SymplecticEulerPush{IT} <: IntegrationStep
    species_name::String
    species_index::IT
end

function SymplecticEulerPush(species_name::String)

    species_index = -1

    return SymplecticEulerPush(species_name, species_index)
end

function setup!(step::SymplecticEulerPush, sim::Simulation)
    step.species_index = sim.species_dict[step.species_name]
end

function step!(step::SymplecticEulerPush{IT}, sim::Simulation{FT, IT}) where {FT, IT}
    positions = sim.species[step.species_index].positions
    velocities = sim.species[step.species_index].velocities
    forces = sim.species[step.species_index].forces
    mass = sim.species[step.species_index].macro_mass
    timestep = sim.timestep

    symplectic_euler_push!(positions, velocities, forces, mass, timestep)

    return
end

function symplectic_euler_push!(positions, velocities, forces, mass, timestep)
    velocities .+= (timestep / mass) .* forces
    positions  .+= timestep .* velocities

    return
end


mutable struct ConstrainSpecies{IT, G} <: IntegrationStep
    species_name::String
    species_index::IT
    grid::G
end

function ConstrainSpecies(species_name::String, grid::G) where {G <: AbstractGrid}
    species_index = -1

    return ConstrainSpecies(species_name, species_index, grid)
end

function setup!(step::ConstrainSpecies, sim::Simulation)
    step.species_index = sim.species_dict[step.species_name]
end

function step!(step::ConstrainSpecies, sim::Simulation)
    positions = sim.species[step.species_index].positions
    sim_length = simulation_length(step.grid)

    for i in 1:length(positions)
        if positions[i] >= sim_length
            positions[i] -= sim_length * floor(positions[i] / sim_length)
        elseif positions[i] < 0
            positions[i] -= sim_length * floor(positions[i] / sim_length)
        end
    end

    return
end
