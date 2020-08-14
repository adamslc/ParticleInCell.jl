mutable struct Species{FT, IT}
    num_macroparticles::IT
    physical_charge::FT
    physical_mass::FT

    physical_particles_per_macroparticle::IT

    macro_charge::FT
    macro_mass::FT

    positions::Vector{FT}
    velocities::Vector{FT}
    forces::Vector{FT}

    function Species(physical_charge::FT, physical_mass::FT,
                     physical_particles_per_macroparticle::IT,
                     positions, velocities, forces) where {FT, IT}
        @assert eltype(positions) == FT
        @assert eltype(velocities) == FT
        @assert eltype(forces) == FT
        @assert length(positions) == length(velocities)
        @assert length(positions) == length(forces)

        num_macroparticles = length(positions)
        macro_charge = physical_particles_per_macroparticle * physical_charge
        macro_mass = physical_particles_per_macroparticle * physical_mass

        new{FT, IT}(num_macroparticles, physical_charge, physical_mass,
                    physical_particles_per_macroparticle, macro_charge, macro_mass,
                    positions, velocities, forces)
    end
end
