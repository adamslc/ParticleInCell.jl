mutable struct FourierFieldSolve{FT, IT} <: IntegrationStep
    input_field_name::String
    output_field_name::String

    input_field_index::IT
    output_field_index::IT

    ksq_inv::Vector{FT}
    ft_vector::Vector{Complex{FT}}

    epsilon_0::FT

    function FourierFieldSolve(input_field_name::String, output_field_name::String,
                               epsilon_0::FT) where FT
        return new{FT, Int}(input_field_name, output_field_name, -1, -1,
                            Vector{FT}(), Vector{Complex{FT}}(), epsilon_0)
    end
end

function setup!(step::FourierFieldSolve, sim::Simulation)
    step.input_field_index = sim.fields_dict[step.input_field_name]
    step.output_field_index = sim.fields_dict[step.output_field_name]

    @assert sim.fields[step.input_field_index].grid ==
            sim.fields[step.output_field_index].grid
    @assert typeof(sim.fields[step.input_field_index].grid) <: UniformGrid

    num_cells = sim.fields[step.input_field_index].grid.num_cells
    sim_length = simulation_length(sim.fields[step.input_field_index].grid)
    cell_length = sim.fields[step.input_field_index].grid.cell_length

    step.ksq_inv = Vector{Float64}(undef, round(Int, 1 + num_cells / 2))
    step.ft_vector = Vector{Complex{Float64}}(undef, num_cells)

    for i in 1:round(Int, 1 + num_cells / 2)
        k = 2 * pi * i / sim_length
        # grid_angle can alternatively be computed as i * pi / num_cells. I
        # don't do this for intuitions sake.
        grid_angle = k * cell_length / 2
        step.ksq_inv[i] = (cell_length / (2 * sin(grid_angle)))^2 / step.epsilon_0
    end
end

function step!(step::FourierFieldSolve, sim::Simulation)
    input_values = sim.fields[step.input_field_index].grid_values
    output_values = sim.fields[step.output_field_index].grid_values
    grid = sim.fields[step.input_field_index].grid
    num_cells = sim.fields[step.input_field_index].grid.num_cells
    num_guard_cells = sim.fields[step.input_field_index].grid.num_guard_cells

    offset = num_guard_cells
    # Setup ft_vector to Fourier tramsform
    for i in 1:num_cells
        step.ft_vector[i] = input_values[i + offset]
    end

    # Compute FFT in place
    fft!(step.ft_vector)

    # Modifiy ft_vector to get phi(k)
    # This mode should be zero if the simulation cell is charge neutral. We
    # assume that it is here.
    step.ft_vector[1] = 0

    # The modes 2:num_cells (of which there are num_cells - 1) are paired in conjugate
    # pairs so rhok[2] = complex_conjugate(rhok[end]) and so on. Thus we need
    # to multiply both by the same scaling factor.
    for i in 1:round(Int, num_cells / 2)
        k1 = i + 1
        k2 = num_cells + 1 - i
        step.ft_vector[k1] = step.ksq_inv[i] * step.ft_vector[k1]
        k1 == k2 && break
        step.ft_vector[k2] = step.ksq_inv[i] * step.ft_vector[k2]
    end

    # Compute IFFT in place
    ifft!(step.ft_vector)

    # Write values to output_field
    for i in 1:num_cells
        output_values[i + offset] = real(step.ft_vector[i])
    end

    # Update guard/periodic cells
    for j in 1:grid.num_guard_cells
        output_values[j] = output_values[grid.num_cells + j]
    end
    output_values[grid.num_guard_cells + grid.num_cells + 1] = output_values[grid.num_guard_cells + 1]
    for j in 1:grid.num_guard_cells
        output_values[grid.num_guard_cells + grid.num_cells + 1 + j] = output_values[grid.num_guard_cells + 1 + j]
    end
end

function fourier_field_solve!(f1, f2, ft_vector, ksq_inv, offset, num_cells)
    # Setup ft_vector to Fourier tramsform
    for i in 1:num_cells
        ft_vector[i] = f1[i + offset]
    end

    # Compute FFT in place
    fft!(ft_vector)

    # We assume that the plasma is charge neutral, so this mode must be zero.
    ft_vector[1] = 0

    # The modes 2:num_cells (of which there are num_cells - 1) are paired in conjugate
    # pairs so rhok[2] = complex_conjugate(rhok[end]) and so on. Thus we need
    # to multiply both by the same scaling factor.
    for i in 1:round(Int, num_cells / 2)
        k1 = i + 1
        k2 = num_cells + 1 - i
        ft_vector[k1] = ksq_inv[i] * ft_vector[k1]
        k1 == k2 && break
        ft_vector[k2] = ksq_inv[i] * ft_vector[k2]
    end

    # Compute IFFT in place
    ifft!(ft_vector)

    # Write values to output_field
    for i in 1:num_cells
        f2[i + offset] = real(ft_vector[i])
    end

    ## Update guard/periodic cells
    #for j in 1:grid.num_guard_cells
    #    f2[j] = f2[grid.num_cells + j]
    #end
    #f2[grid.num_guard_cells + grid.num_cells + 1] = f2[grid.num_guard_cells + 1]
    #for j in 1:grid.num_guard_cells
    #    f2[grid.num_guard_cells + grid.num_cells + 1 + j] = f2[grid.num_guard_cells + 1 + j]
    #end
end

mutable struct FiniteDifferenceToNodes{IT} <: IntegrationStep
    input_field_name::String
    output_field_name::String

    input_field_index::IT
    output_field_index::IT
end

FiniteDifferenceToNodes(ifn::String, ofn::String) = FiniteDifferenceToNodes(ifn, ofn, -1, -1)

function setup!(step::FiniteDifferenceToNodes, sim::Simulation)
    step.input_field_index = sim.fields_dict[step.input_field_name]
    step.output_field_index = sim.fields_dict[step.output_field_name]
end

function step!(step::FiniteDifferenceToNodes, sim::Simulation)
    input_values = sim.fields[step.input_field_index].grid_values
    output_values = sim.fields[step.output_field_index].grid_values
    grid = sim.fields[step.input_field_index].grid
    num_cells = sim.fields[step.input_field_index].grid.num_cells
    num_guard_cells = sim.fields[step.input_field_index].grid.num_guard_cells

    for i in 2:length(input_values) - 1
        output_values[i] = -1 * (input_values[i + 1] - input_values[i - 1]) / (2 * cell_length(grid, i))
    end

    # Update guard/periodic cells
    for j in 1:grid.num_guard_cells
        output_values[j] = output_values[grid.num_cells + j]
    end
    output_values[grid.num_guard_cells + grid.num_cells + 1] = output_values[grid.num_guard_cells + 1]
    for j in 1:grid.num_guard_cells
        output_values[grid.num_guard_cells + grid.num_cells + 1 + j] = output_values[grid.num_guard_cells + 1 + j]
    end
end

"""
    d2_finite_diff(vec, num_cells, sim_length)

Computes a finite differenced second derivative of `vec`. It is assumed that
`vec` is periodic, and that `vec[1] == vec[end]` so that
`length(vec) == num_cells - 1`.
"""
function d2_finite_diff(vec, num_cells, sim_length)
    cell_length = sim_length / num_cells

    diff = similar(vec)
    for i in 2:num_cells
        diff[i] = (vec[i-1] - 2*vec[i] + vec[i + 1]) / cell_length^2
    end

    diff[1] = (vec[end-1] - 2*vec[1] + vec[2]) / cell_length^2
    diff[end] = diff[1]

    return diff
end
