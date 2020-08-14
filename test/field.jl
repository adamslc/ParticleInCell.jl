using ParticleInCell
using Test

@testset "Field" begin
    num_cells = 32
    num_guard_cells = 2
    simulation_length = 10.
    g = ParticleInCell.UniformGrid(num_cells, num_guard_cells, simulation_length)

    f = ParticleInCell.Field(g)
    @test ParticleInCell.total_cells(g) == num_cells + 2*num_guard_cells + 1
    @test length(f.grid_values) == ParticleInCell.total_cells(g)
end
