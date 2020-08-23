using ParticleInCell
using HDF5

timestep = 1e-10
sim = Simulation(timestep)

sim_length = 1.
num_cells = 32
num_guard_cells = 2
grid = UniformGrid(num_cells=num_cells, num_guard_cells=num_guard_cells, simulation_length=sim_length)

grid_multiplier = 5
num_cells_hd = num_cells * grid_multiplier
num_guard_cells_hd = num_guard_cells * grid_multiplier
grid_hd = UniformGrid(num_cells=num_cells_hd, num_guard_cells=num_guard_cells_hd, simulation_length=sim_length)

rho = Field(grid)
add_field!(sim, rho, "rho")
phi = Field(grid)
add_field!(sim, phi, "phi")
elec = Field(grid)
add_field!(sim, elec, "elec")

rho_hd = Field(grid_hd)
add_field!(sim, rho_hd, "rho_hd")

electron_charge = -1.602e-19
electron_mass = 9.109e-31
particles_per_macro = 10^10
num_macros = 6400
positions  = collect(range(0, stop=sim_length, length=num_macros+1))[1:num_macros]
velocities = 1e4 * sin.(2*pi*2/sim_length .* positions)
forces     = similar(positions)
electrons = Species(electron_charge, electron_mass, particles_per_macro,
                                      positions, velocities, forces)
add_species!(sim, electrons, "electrons")

scatter = ScatterChargeToGrid("electrons", "rho")
#scatter = ScatterChargeToGrid("electrons", "rho", shape_2nd_order)
add_integration_step!(sim, scatter)

scatter_hd = ScatterChargeToGrid("electrons", "rho_hd")
add_integration_step!(sim, scatter_hd)

epsilon_0 = 8.85e-12
field_solve = FourierFieldSolve("rho", "phi", epsilon_0)
add_integration_step!(sim, field_solve)

finite_diff = FiniteDifferenceToNodes("phi", "elec")
add_integration_step!(sim, finite_diff)

gather = GatherForcesFromGrid("electrons", "elec")
#gather = GatherForcesFromGrid("electrons", "elec", shape_2nd_order)
add_integration_step!(sim, gather)

symplectic_euler = SymplecticEulerPush("electrons")
add_integration_step!(sim, symplectic_euler)

constrain_species = ConstrainSpecies("electrons", grid)
add_integration_step!(sim, constrain_species)

setup!(sim)

num_steps = 1000
dump_steps = 10
background_charge_density = -1 * electron_charge * particles_per_macro * num_macros / num_cells
background_charge_density_hd = -1 * electron_charge * particles_per_macro * num_macros / num_cells_hd

rho.grid_values .= 0.
rho.grid_values[num_guard_cells + 1:num_guard_cells + num_cells] .= background_charge_density
electrons.forces .= 0

grid_pos = sim_length / num_cells * (collect(1:num_cells + 2*num_guard_cells + 1) .- (num_guard_cells + 1))
grid_pos_hd = sim_length / num_cells_hd * (collect(1:num_cells_hd + 2*num_guard_cells_hd + 1) .- (num_guard_cells_hd + 1))

mkpath("data")
h5open("data/dump_0.h5", "w") do file
    write(file, "positions", electrons.positions)
    write(file, "velocities", electrons.velocities)
    write(file, "forces", electrons.forces)

    write(file, "grid_pos", grid_pos)
    write(file, "grid_pos_hd", grid_pos_hd)
    write(file, "rho", rho.grid_values)
    write(file, "rho_hd", rho_hd.grid_values)
    write(file, "phi", phi.grid_values)
    write(file, "elec", elec.grid_values)
end

for i in 1:num_steps
    rho.grid_values .= 0.
    rho.grid_values[num_guard_cells + 1:num_guard_cells + num_cells] .= background_charge_density
    rho_hd.grid_values .= 0.
    rho_hd.grid_values[num_guard_cells_hd + 1:num_guard_cells_hd + num_cells_hd] .= background_charge_density_hd
    electrons.forces .= 0
    step!(sim)

    #@info "fields $i" rho.grid_values sum(rho.grid_values[num_guard_cells+1:num_cells+num_guard_cells]) elec.grid_values
    #@info "electrons $i" electrons.forces

    if i % dump_steps == 0
        dump_number = div(i, dump_steps)

        h5open("data/dump_$dump_number.h5", "w") do file
            write(file, "positions", electrons.positions)
            write(file, "velocities", electrons.velocities)
            write(file, "forces", electrons.forces)

            write(file, "grid_pos", grid_pos)
            write(file, "grid_pos_hd", grid_pos_hd)
            write(file, "rho", rho.grid_values)
            write(file, "rho_hd", rho_hd.grid_values)
            write(file, "phi", phi.grid_values)
            write(file, "elec", elec.grid_values)
        end
    end
end
