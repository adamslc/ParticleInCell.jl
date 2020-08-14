using ParticleInCell
using Test

@testset "FourierFieldSolve" begin
    offset = 0
    num_cells = 32
    sim_length = 1.
    cell_length = sim_length / num_cells
    grid_pos = collect(range(0, stop=sim_length, length=num_cells+1))[1:num_cells]

    ft_vector = Vector{Complex{Float64}}(undef, num_cells)
    ksq_inv = Vector{Float64}(undef, round(Int, 1 + num_cells / 2))

    # Compute ksq_inv with epsilon_0 = 1
    for i in 1:round(Int, 1 + num_cells / 2)
        k = 2 * pi * i / sim_length
        grid_angle = k * cell_length / 2
        ksq_inv[i] = (cell_length / (2 * sin(grid_angle)))^2
    end

    f1 = similar(grid_pos)
    f2 = similar(f1)
    expected_sol = similar(f1)

    # If rho = sin((2*pi/L)*k*x), then the solution to the equation
    # d^2/dx^2(phi) = -1*rho is phi = (L/(2*pi*k))^2 sin(2*pi/L)*k*L). We check
    # that this is true for a range of k's here.
    # TODO: This actually computes a solution to the exact finite differenced
    # result, so I should check against that.
    for wave_num in 1:8
        f1 .= sin.(2*pi*wave_num/sim_length .* grid_pos)
        expected_sol .= (sim_length / (2*pi*wave_num))^2 .* sin.(2*pi*wave_num/sim_length .* grid_pos)

        ParticleInCell.fourier_field_solve!(f1, f2, ft_vector, ksq_inv, offset, num_cells)

        @test all(isapprox.(f2, expected_sol, atol=0.001))
    end

    # Compute ksq_inv with epsilon_0 = 0.001
    epsilon_0 = 0.001
    for i in 1:round(Int, 1 + num_cells / 2)
        k = 2 * pi * i / sim_length
        grid_angle = k * cell_length / 2
        ksq_inv[i] = (cell_length / (2 * sin(grid_angle)))^2 / epsilon_0
    end

    # Check that the field solve still correct for
    # d^2/dx^2(phi) = -1*rho/epsilon_0 where epsilon_0 = 0.001
    # TODO: See above.
    for wave_num in 1:8
        f1 .= sin.(2*pi*wave_num/sim_length .* grid_pos)
        expected_sol .= (sim_length / (2*pi*wave_num))^2 .* sin.(2*pi*wave_num/sim_length .* grid_pos) ./ epsilon_0

        ParticleInCell.fourier_field_solve!(f1, f2, ft_vector, ksq_inv, offset, num_cells)

        @test all(isapprox.(f2, expected_sol, atol=0.1))
    end
end
