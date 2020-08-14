using ParticleInCell
using Test

@testset "ParticleInCell" begin
    include("grid.jl")
    include("scatter.jl")
    include("species.jl")
    include("field.jl")
    include("field_solve.jl")
end
