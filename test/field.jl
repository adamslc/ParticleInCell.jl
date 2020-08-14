using ParticleInCell
using Test

@testset "Field" begin
    num_cells = 32
    num_guard_cells = 2
    simulation_length = 10.
    g = UniformGrid(num_cells, num_guard_cells, simulation_length)

    f = Field(g)
    @test total_cells(g) == num_cells + 2*num_guard_cells + 1
    @test length(f.grid_values) == total_cells(g)
end
