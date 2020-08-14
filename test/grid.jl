using ParticleInCell
using Test

@testset "UniformGrid" begin
    num_cells = 32
    num_guard_cells = 2
    simulation_length = 10.
    g = UniformGrid(num_cells, num_guard_cells, simulation_length)
    @test cell_length(g) == simulation_length / num_cells

    cells = collect(g)
    @test cells[1][1] == 1
    @test cells[1][2] == -1 * num_guard_cells * cell_length(g)

    @test cells[num_guard_cells+1][2] == 0.
    @test cells[num_guard_cells+num_cells+1][2] == simulation_length
end
