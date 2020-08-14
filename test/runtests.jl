using ParticleInCell
using Test

using Random
Random.seed!(1)

@testset "ParticleInCell" begin
    include("grid.jl")
    include("scatter.jl")
    include("species.jl")
    include("field.jl")
    include("field_solve.jl")
end
