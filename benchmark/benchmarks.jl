using ElectrostaticTest
using BenchmarkTools

const SUITE = BenchmarkGroup()

num_cells = 32
num_guard_cells = 2
simulation_length = 1.
g1 = ElectrostaticTest.UniformGrid(32, 2, 1.)
f1 = fill(0., ElectrostaticTest.total_cells(g1))

num_macros = 640
p1  = collect(range(0, stop=simulation_length, length=num_macros+1))[1:num_macros]

num_cells = 128
num_guard_cells = 2
simulation_length = 1.
g2 = ElectrostaticTest.UniformGrid(32, 2, 1.)
f2 = fill(0., ElectrostaticTest.total_cells(g1))

num_macros = 10000
p2  = collect(range(0, stop=simulation_length, length=num_macros+1))[1:num_macros]

shape_functions = [ElectrostaticTest.shape_1st_order,
                   ElectrostaticTest.shape_2nd_order]

SUITE["scatter"] = BenchmarkGroup()
for shape_function in shape_functions
    SUITE["scatter"][string(shape_function), 1] =
        @benchmarkable ElectrostaticTest.scatter_charge_to_grid!($p1, $f1, $g1, 1., $shape_function)
    SUITE["scatter"][string(shape_function), 2] =
        @benchmarkable ElectrostaticTest.scatter_charge_to_grid!($p2, $f2, $g2, 1., $shape_function)
end





offset = 0
num_cells = 32
sim_length = 1.
wave_num = 2
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

f1 = sin.(2*pi*wave_num/sim_length .* grid_pos)
f2 = similar(f1)

SUITE["field_solve"] = BenchmarkGroup()
SUITE["field_solve"]["fourier"] =
    @benchmarkable ElectrostaticTest.fourier_field_solve!(f1, f2, ft_vector, ksq_inv, offset, num_cells)