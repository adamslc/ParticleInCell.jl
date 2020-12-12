mutable struct ParticleInFourier{FT, IT} <: IntegrationStep
    species_name::String
    species_index::IT

    num_modes::IT
    modes::Vector{FT}

    rho::Vector{Complex{FT}}
    E::Vector{Complex{FT}}

    epsilon_0::FT

    function ParticleInFourier(species_name::String, modes::Vector{FT}, epsilon_0::FT) where FT
        num_modes = length(modes)
        IT = typeof(num_modes)
        rho = Vector{Complex{FT}}(undef, num_modes)
        E = similar(rho)

        new{FT, IT}(species_name, -1, num_modes, modes, rho, E, epsilon_0)
    end
end

function setup!(step::ParticleInFourier, sim::Simulation)
    step.species_index = sim.species_dict[step.species_name]
end

# Sk(kΔx) = (sin(kΔx/2) / (kΔx/2))^2

function step!(step::ParticleInFourier, sim::Simulation)
    positions = sim.species[step.species_index].positions
    forces    = sim.species[step.species_index].forces
    charge    = sim.species[step.species_index].macro_charge

    rho   = step.rho
    E     = step.E
    modes = step.modes

    rho .= 0.
    E   .= 0.

    Sk = 1

    # Compute rhos
    for position in positions, (i, mode) in enumerate(modes)
        rho[i] += charge * Sk * exp(-1 * im * mode * position)
    end

    # Use rho to compute E
    for (i, mode) in enumerate(modes)
        E[i] = rho[i] / (im * mode * step.epsilon_0)
    end

    # Compute force on each particle
    for (i, position) in enumerate(positions), (j, mode) in enumerate(modes)
        forces[i] += real(charge * Sk * E[j] * exp(im * mode * position))
    end
end
