mutable struct CubicSplineFieldSolve{FT, IT} <: IntegrationStep
    rho0_name::String
    rho2_name::String

    phi0_name::String
    phi2_name::String

    rho0_index::IT
    rho2_index::IT

    phi0_index::IT
    phi2_index::IT

    epsilon_0::FT

    ft_vector1::Vector{Complex{FT}}
    ft_vector2::Vector{Complex{FT}}

    sim_length::FT

    ksq_inv::Vector{FT}

    function CubicSplineFieldSolve(rho0_name::String, rho2_name::String,
                                   phi0_name::String, phi2_name::String,
                                   epsilon_0::FT) where FT
        IT = Int

        new{FT, IT}(rho0_name, rho2_name, phi0_name, phi2_name, -1, -1, -1, -1, epsilon_0)
    end
end

function setup!(step::CubicSplineFieldSolve, sim::Simulation)
    step.rho0_index = sim.fields_dict[step.rho0_name]
    step.rho2_index = sim.fields_dict[step.rho2_name]
    step.phi0_index = sim.fields_dict[step.phi0_name]
    step.phi2_index = sim.fields_dict[step.phi2_name]

    grid = sim.fields[step.rho0_index].grid

    step.ft_vector1 = Vector{Complex{Float64}}(undef, grid.num_cells)
    step.ft_vector2 = Vector{Complex{Float64}}(undef, grid.num_cells)

    step.sim_length = simulation_length(grid)

    num_cells = grid.num_cells
    cell_length = grid.cell_length
    step.ksq_inv = Vector{Float64}(undef, round(Int, 1 + num_cells / 2))
    for i in 1:round(Int, 1 + num_cells / 2)
        k = 2 * pi * i / step.sim_length
        # grid_angle can alternatively be computed as i * pi / num_cells. I
        # don't do this for intuitions sake.
        grid_angle = k * cell_length / 2
        step.ksq_inv[i] = (cell_length / (2 * sin(grid_angle)))^2 / step.epsilon_0
    end
end

function step!(step::CubicSplineFieldSolve, sim::Simulation)
    rho0 = sim.fields[step.rho0_index].grid_values
    rho2 = sim.fields[step.rho2_index].grid_values
    phi0 = sim.fields[step.phi0_index].grid_values
    phi2 = sim.fields[step.phi2_index].grid_values

    grid = sim.fields[step.rho0_index].grid
    num_cells = grid.num_cells
    num_guard_cells = grid.num_guard_cells
    offset = num_guard_cells

    # TODO: is this right?
    # We need to make rho2 charge neutral
    sum = 0.
    for i in 1:num_cells
        sum += rho2[i + offset]
    end
    rho2 .-= sum / num_cells

    # TODO: same as above
    sum = 0.
    for i in 1:num_cells
        sum += rho0[i + offset]
    end
    rho0 .-= sum / num_cells

    # Set up transformtion vectors
    for i in 1:num_cells
        step.ft_vector1[i] = rho0[i + offset]
        step.ft_vector2[i] = rho2[i + offset]
    end

    # Compute FFTs in place
    fft!(step.ft_vector1)
    fft!(step.ft_vector2)

    step.ft_vector1[1] = 0
    step.ft_vector2[1] = 0

    for i in 1:round(Int, num_cells / 2)
        k = 2 * pi * i / step.sim_length
        k1 = i + 1
        k2 = num_cells + 1 - i
        step.ft_vector1[k1] = (step.ft_vector1[k1] - k^2 * step.ft_vector2[k1]) * step.ksq_inv[i]
        k1 == k2 && break
        step.ft_vector1[k2] = (step.ft_vector1[k2] - k^2 * step.ft_vector2[k2]) * step.ksq_inv[i]
    end

    for i in 1:round(Int, num_cells / 2)
        k = 2 * pi * i / step.sim_length
        k1 = i + 1
        k2 = num_cells + 1 - i
        step.ft_vector2[k1] = -1 * k^2 * step.ft_vector1[k1]
        k1 == k2 && break
        step.ft_vector2[k2] = -1 * k^2 * step.ft_vector1[k2]
    end

    # Compute phi inverses in place
    ifft!(step.ft_vector1)
    ifft!(step.ft_vector2)

    for i in 1:num_cells
        phi0[i + offset] = real(step.ft_vector1[i])
        phi2[i + offset] = real(step.ft_vector2[i])
    end

    # Update guard/periodic cells
    for j in 1:num_guard_cells
        phi0[j] = phi0[num_cells + j]
        phi2[j] = phi2[num_cells + j]
    end
    phi0[num_guard_cells + num_cells + 1] = phi0[num_guard_cells + 1]
    phi2[num_guard_cells + num_cells + 1] = phi2[num_guard_cells + 1]
    for j in 1:num_guard_cells
        phi0[num_guard_cells + num_cells + 1 + j] = phi0[num_guard_cells + 1 + j]
        phi2[num_guard_cells + num_cells + 1 + j] = phi2[num_guard_cells + 1 + j]
    end
end
