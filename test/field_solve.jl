using ParticleInCell
using Test

# Helper function to test is solutions satisfy the Poisson equation
function compute_poisson_derivative(phi, num_cells, offset,
        cell_length, epsilon_0)
    rho = similar(phi)
    for i in 2:num_cells-1
        rho[i] = (phi[i-1] - 2*phi[i] + phi[i+1]) / cell_length^2
    end

    # TODO: This needs to take the offset into account
    rho[1]   = (phi[end]   - 2*phi[1]   + phi[2]) / cell_length^2
    rho[end] = (phi[end-1] - 2*phi[end] + phi[1]) / cell_length^2

    rho .*= -1 * epsilon_0

    return rho
end

function compute_ksq_inv!(ksq_inv, num_cells, sim_length, epsilon_0)
    cell_length = sim_length / num_cells

    for i in 1:round(Int, 1 + num_cells / 2)
        k = 2 * pi * i / sim_length
        grid_angle = k * cell_length / 2
        ksq_inv[i] = (cell_length / (2 * sin(grid_angle)))^2 / epsilon_0
    end

    return
end

@testset "FourierFieldSolve" begin
@testset "Internals" begin
    offset = 0
    num_cells = 32
    sim_length = 1.
    cell_length = sim_length / num_cells

    grid_pos = collect(range(0, stop=sim_length, length=num_cells+1))[1:num_cells]
    f1 = similar(grid_pos)
    f2 = similar(grid_pos)

    ft_vector = Vector{Complex{Float64}}(undef, num_cells)
    ksq_inv = Vector{Float64}(undef, round(Int, 1 + num_cells / 2))

    # First a test with natural units
    # If rho = sin((2*pi/L)*k*x), then the solution to the equation
    # d^2/dx^2(phi) = -1*rho is phi = (L/(2*pi*k))^2 sin(2*pi/L)*k*L). We check
    # that this is true for a range of k's here.
    epsilon_0 = 1.
    compute_ksq_inv!(ksq_inv, num_cells, sim_length, epsilon_0)
    for wave_num in 1:8
        f1 .= sin.(2*pi*wave_num/sim_length .* grid_pos)

        ParticleInCell.fourier_field_solve!(f1, f2, ft_vector, ksq_inv, offset, num_cells)

        computed_f1 = compute_poisson_derivative(f2, num_cells, offset,
                cell_length, 1.)

        @test all(isapprox.(f1, computed_f1, atol=1e-8))
    end

    # Now a test with an arbitrary non-unit epsilon_0
    epsilon_0 = 0.001
    compute_ksq_inv!(ksq_inv, num_cells, sim_length, epsilon_0)
    for wave_num in 1:8
        f1 .= sin.(2*pi*wave_num/sim_length .* grid_pos)

        ParticleInCell.fourier_field_solve!(f1, f2, ft_vector, ksq_inv, offset, num_cells)

        computed_f1 = compute_poisson_derivative(f2, num_cells, offset,
                cell_length, epsilon_0)

        @test all(isapprox.(f1, computed_f1, atol=1e-8))
    end
end

@testset "Integration Step" begin
    sim = Simulation(1.)

    grid = UniformGrid(num_cells=32, num_guard_cells=4, simulation_length=8.)

    f1 = Field(grid)
    add_field!(sim, f1, "f1")

    f2 = Field(grid)
    add_field!(sim, f2, "f2")

    ffs = FourierFieldSolve("f1", "f2", 1.)
    add_integration_step!(sim, ffs)

    setup!(sim)
    @test sim.fields[ffs.input_field_index]  === f1
    @test sim.fields[ffs.output_field_index] === f2

    step!(sim)
    # TODO: Need to test that f2 actually satisfies the Poisson equation
end
end
