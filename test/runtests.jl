using ParticleInCell
using Test

struct FakeStep <: ParticleInCell.IntegrationStep end

@testset "ParticleInCell" begin
    sim = Simulation(1.)
    add_integration_step!(sim, FakeStep())

    setup!(sim)

    @test_throws String step!(sim)

    include("grid.jl")
    include("scatter.jl")
    include("species.jl")
    include("field.jl")
    include("field_solve.jl")
end
