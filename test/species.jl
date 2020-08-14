using ParticleInCell
using Test

@testset "Species" begin
    physical_charge = 1.
    physical_mass = 1.
    physical_particles_per_macroparticle = 10^2

    num_macroparticles = 10
    positions = zeros(num_macroparticles)
    velocities = zeros(num_macroparticles)
    forces = zeros(num_macroparticles)

    s = ParticleInCell.Species(physical_charge, physical_mass,
                                  physical_particles_per_macroparticle,
                                  positions, velocities, forces)

    @test s.num_macroparticles == num_macroparticles
    @test s.macro_charge == physical_particles_per_macroparticle * physical_charge
    @test s.macro_mass == physical_particles_per_macroparticle * physical_mass

    forces = zeros(num_macroparticles - 1)
    @test_throws AssertionError ParticleInCell.Species(physical_charge, physical_mass,
                                  physical_particles_per_macroparticle,
                                  positions, velocities, forces)
end
