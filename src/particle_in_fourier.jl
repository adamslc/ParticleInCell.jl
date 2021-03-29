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

mutable struct ScatterToFourierGrid{IT, SF} <: IntegrationStep
    rho_name::String
    rho_index::IT

    species_name::String
    species_index::IT

    shape_function::SF
end

function ScatterToFourierGrid(species_name, rho_name, shape_function=x->1)
    return ScatterToFourierGrid(rho_name, -1, species_name, -1, shape_function)
end

function setup!(step::ScatterToFourierGrid, sim::Simulation)
    step.rho_index     = sim.fields_dict[step.rho_name]
    step.species_index = sim.species_dict[step.species_name]
end

function step!(step::ScatterToFourierGrid, sim::Simulation)
    rho        = sim.fields[step.rho_index]
    positions  = sim.species[step.species_index].positions
    num_cells  = total_cells(rho.grid)
    sim_length = simulation_length(rho.grid)
    charge     = sim.species[step.species_index].macro_charge

    # Reset rho to zero
    # This PIF implementation is only valid for periodic boundaries where
    # there can be no constant electric force that results from the
    # self-consistant interactions. Therefore, we keep rho[1] at zero.
    rho.grid_values .= 0

    for i in 1:round(Int, num_cells / 2)
        k1 = i + 1
        k2 = num_cells + 1 - i
        k  = 2 * pi * i / sim_length

        # FrequencyGrid stores values to be compliant with FFTW. For real
        # values, this means that indices k1 and k2 are conjugate pairs.
        for position in positions
            rho.grid_values[k1] += charge * step.shape_function(k) * exp(-1 * im * k * position)
            rho.grid_values[k2] += charge * step.shape_function(k) * exp(     im * k * position)
        end
    end
end

mutable struct PIFFieldSolve{FT, IT} <: IntegrationStep
    rho_name::String
    rho_index::IT

    elec_name::String
    elec_index::IT

    epsilon_0::FT
end

function PIFFieldSolve(rho_name, elec_name, epsilon_0)
    return PIFFieldSolve(rho_name, -1, elec_name, -1, epsilon_0)
end

function setup!(step::PIFFieldSolve, sim::Simulation)
    step.rho_index  = sim.fields_dict[step.rho_name]
    step.elec_index = sim.fields_dict[step.elec_name]
end

function step!(step::PIFFieldSolve, sim::Simulation)
    rho        = sim.fields[step.rho_index]
    elec       = sim.fields[step.elec_index]
    num_cells  = total_cells(rho.grid)
    sim_length = simulation_length(rho.grid)

    for i in 1:round(Int, num_cells / 2)
        k1   = i + 1
        k2   = num_cells + 1 - i
        mode = i / sim_length

        elec.grid_values[k1] = rho.grid_values[k1] / (im * mode * step.epsilon_0)
        elec.grid_values[k2] = rho.grid_values[k2] / (im * mode * step.epsilon_0)
    end
end

mutable struct GatherFromFourierGrid{FT, IT, SF} <: IntegrationStep
    elec_name::String
    elec_index::IT

    species_name::String
    species_index::IT

    shape_function::SF
end

function GatherFromFourierGrid(species_name, elec_name, shape_function=x->1)
    return ScatterToFourierGrid(elec_name, -1, species_name, -1, shape_function)
end

function setup!(step::GatherFromFourierGrid, sim::Simulation)
    step.elec_index     = sim.fields_dict[step.elec_name]
    step.species_index = sim.species_dict[step.species_name]
end

function step!(step::GatherFromFourierGrid, sim::Simulation)
    elec       = sim.fields[step.elec_index]
    positions  = sim.species[step.species_index].positions
    forces     = sim.species[step.species_index].forces
    num_cells  = total_cells(rho.grid)
    sim_length = simulation_length(rho.grid)
    charge     = sim.species[step.species_index].macro_charge

    for i in 1:round(Int, num_cells / 2)
        k1 = i + 1
        k2 = num_cells + 1 - i
        k  = 2 * pi * i / sim_length

        # FrequencyGrid stores values to be compliant with FFTW. For real
        # values, this means that indices k1 and k2 are conjugate pairs.
        for (j, position) in enumerate(positions)
            forces[j] .+= real(charge * elec[k1] * step.shape_function(k) * exp(-1 * im * k * position))
            forces[j] .+= real(charge * elec[k2] * step.shape_function(k) * exp(     im * k * position))
        end
    end
end
