using ParticleInCell
using HDF5

function create_electron_velocities(wavenumber, velocity_amplitude, sim_length, positions; interpolated=false, num_cells=32)
    if interpolated
        grid_pos = sim_length / num_cells * (collect(1:num_cells + 1) .- 1)
        grid_vel = velocity_amplitude * sin.(2*pi*wavenumber/sim_length .* grid_pos)

        velocities = zeros(length(positions))

        for (i, _) in enumerate(velocities), (j, _) in enumerate(grid_pos)
            velocities[i] += grid_vel[j] * shape_1st_order((positions[i] - grid_pos[j]) * num_cells / sim_length)
        end

        return velocities
    else
        return velocity_amplitude * sin.(2*pi*wavenumber/sim_length .* positions)
    end
end

function add_electrons!(sim, num_macros, charge_density, sim_length, wavenumber, velocity_amplitude)
    electron_charge = -1.602e-19
    electron_mass = 9.109e-31

    particles_per_macro = round(Int, charge_density * sim_length / (electron_charge * num_macros))
    @show particles_per_macro

    positions  = collect(range(0, stop=sim_length, length=num_macros+1))[1:num_macros]
    velocities = create_electron_velocities(wavenumber, velocity_amplitude, sim_length, positions, interpolated=false)
    forces     = similar(positions)
    electrons = Species(electron_charge, electron_mass, particles_per_macro,
                                          positions, velocities, forces)
    add_species!(sim, electrons, "electrons")
    return electrons
end

function add_electrostatic_fields!(sim, num_cells, num_guard_cells, sim_length)
    grid = UniformGrid(num_cells=num_cells, num_guard_cells=num_guard_cells, simulation_length=sim_length)

    rho = Field(grid)
    add_field!(sim, rho, "rho")
    phi = Field(grid)
    add_field!(sim, phi, "phi")
    elec = Field(grid)
    add_field!(sim, elec, "elec")

    return grid, rho, phi, elec
end

function write_data(sim, data_path, data_prefix, electrons, grid_pos, grid_pos_hd, rho, rho_hd, phi, elec, dump_number, step_number)
    h5open("$data_path/$(data_prefix)_$dump_number.h5", "w") do file
        write(file, "positions", electrons.positions)
        write(file, "velocities", electrons.velocities)
        write(file, "forces", electrons.forces)

        write(file, "grid_pos", grid_pos)
        write(file, "grid_pos_hd", grid_pos_hd)
        write(file, "rho", rho.grid_values)
        write(file, "rho_hd", rho_hd.grid_values)
        write(file, "phi", phi.grid_values)
        write(file, "elec", elec.grid_values)

        write(file, "time", sim.time)
        write(file, "dump_number", dump_number)
        write(file, "step_number", step_number)
    end

    return
end

function do_simulation(;
        timestep=1e-9, num_cells=32, num_guard_cells=2, num_macros=640,
        charge_density=1e-6, wavenumber=2, velocity_amplitude=1e5,
        shape_function=shape_1st_order, prescribed=false,
        num_steps=100, dump_steps=1, data_path="data", data_prefix="dump")

    sim = Simulation(timestep)

    sim_length = 1.
    grid, rho, phi, elec = add_electrostatic_fields!(sim, num_cells, num_guard_cells, sim_length)

    # Add an additional rho grid for higher resolution diagnostics
    grid_multiplier = 5
    num_cells_hd = num_cells * grid_multiplier
    num_guard_cells_hd = num_guard_cells * grid_multiplier
    grid_hd = UniformGrid(num_cells=num_cells_hd, num_guard_cells=num_guard_cells_hd, simulation_length=sim_length)
    rho_hd = Field(grid_hd)
    add_field!(sim, rho_hd, "rho_hd")

    electrons = add_electrons!(sim, num_macros, -1 * charge_density, sim_length, wavenumber, velocity_amplitude)

    # Compute some plasma quantities
    electron_charge = -1.602e-19
    electron_mass = 9.109e-31
    particles_per_macro = round(Int, -1 * charge_density * sim_length / (electron_charge * num_macros))
    epsilon_0 = 8.85e-12
    number_density = particles_per_macro * num_macros
    plasma_freq = sqrt(number_density * electron_charge^2 / (epsilon_0 * electron_mass))
    plasma_period = 2*pi / plasma_freq
    steps_per_period = plasma_period / (timestep)
    steps_per_dump = plasma_period / (dump_steps * timestep)
    angular_wavenumber = 2 * pi * wavenumber
    elec_amp = plasma_freq * electron_mass * velocity_amplitude / abs(electron_charge)
    rho_amp = angular_wavenumber * epsilon_0 * elec_amp
    max_position_change_per_timestep = velocity_amplitude * timestep
    cell_crossings_per_timestep = max_position_change_per_timestep * num_cells / sim_length
    @show plasma_freq plasma_period steps_per_period steps_per_dump elec_amp rho_amp max_position_change_per_timestep cell_crossings_per_timestep
    println()

    scatter = ScatterChargeToGrid("electrons", "rho", shape_function)
    add_integration_step!(sim, scatter)

    scatter_hd = ScatterChargeToGrid("electrons", "rho_hd")
    add_integration_step!(sim, scatter_hd)

    st_func(t, x) = elec_amp * sin(2*pi*2/sim_length * x) * sin(2 * pi / plasma_period * t)

    if prescribed
        prescribe = PrescribeField("elec", st_func)
        add_integration_step!(sim, prescribe)
    else
        field_solve = FourierFieldSolve("rho", "phi", epsilon_0)
        add_integration_step!(sim, field_solve)

        finite_diff = FiniteDifferenceToNodes("phi", "elec")
        add_integration_step!(sim, finite_diff)
    end

    gather = GatherForcesFromGrid("electrons", "elec", shape_function)
    add_integration_step!(sim, gather)

    symplectic_euler = SymplecticEulerPush("electrons")
    add_integration_step!(sim, symplectic_euler)

    constrain_species = ConstrainSpecies("electrons", grid)
    add_integration_step!(sim, constrain_species)

    setup!(sim)

    rho.grid_values .= 0.
    rho.grid_values[num_guard_cells + 1:num_guard_cells + num_cells] .= charge_density
    electrons.forces .= 0
    rho_hd.grid_values .=0

    grid_pos = sim_length / num_cells * (collect(1:num_cells + 2*num_guard_cells + 1) .- (num_guard_cells + 1))
    grid_pos_hd = sim_length / num_cells_hd * (collect(1:num_cells_hd + 2*num_guard_cells_hd + 1) .- (num_guard_cells_hd + 1))

    mkpath(data_path)
    write_data(sim, data_path, data_prefix, electrons, grid_pos, grid_pos_hd, rho, rho_hd, phi, elec, 0, 0)

    for i in 1:num_steps
        rho.grid_values .= 0.
        rho.grid_values[num_guard_cells + 1:num_guard_cells + num_cells] .= charge_density
        rho_hd.grid_values .= 0.
        rho_hd.grid_values[num_guard_cells_hd + 1:num_guard_cells_hd + num_cells_hd] .= charge_density
        electrons.forces .= 0
        elec.grid_values .= 0.

        step!(sim)

        if i % dump_steps == 0
            dump_number = div(i, dump_steps)
            write_data(sim, data_path, data_prefix, electrons, grid_pos, grid_pos_hd, rho, rho_hd, phi, elec, dump_number, i)
        end
    end

    ParticleInCell.print_summary(stdout, sim)

    return
end

# do_simulation(; prescribed=true, velocity_amplitude=1e5, num_steps=1000, num_macros=6400, data_prefix="normal")
# do_simulation(; prescribed=true, velocity_amplitude=1e5, num_steps=1000, shape_function=shape_2nd_order, num_macros=6400, data_prefix="normal2")
# do_simulation(; prescribed=true, velocity_amplitude=1e5, num_steps=1000, shape_function=shape_3rd_order, num_macros=6400, data_prefix="normal3")
# do_simulation(; prescribed=true, velocity_amplitude=1e5, num_cells = 128, data_prefix="fine")
# do_simulation(; prescribed=true, velocity_amplitude=1e5, num_cells = 8, data_prefix="coarse")
# do_simulation(; velocity_amplitude=1e5, num_steps=1000, num_macros=6400, data_prefix="normal-sc")
# do_simulation(; velocity_amplitude=1e5, num_steps=1000, shape_function=shape_2nd_order, num_macros=6400, data_prefix="normal2-sc")
# do_simulation(; velocity_amplitude=1e5, num_steps=10000, shape_function=shape_3rd_order, num_macros=6400, data_prefix="normal3-sc")
# do_simulation(; velocity_amplitude=1e3, wavenumber=18, num_steps=1000, shape_function=shape_1st_order, num_macros=640, data_prefix="foo")
# do_simulation(; velocity_amplitude=1e3, wavenumber=18, prescribed=true, num_steps=1000, shape_function=shape_1st_order, num_macros=640, data_prefix="bar")
do_simulation(; num_steps=100, num_macros=640, data_prefix="testing")
