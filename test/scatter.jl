using ParticleInCell
using Test

@testset "Shape functions" begin
    shape_functions = [ParticleInCell.shape_1st_order,
                       ParticleInCell.shape_2nd_order]

    for shape_function in shape_functions
        for offset in 0:0.1:1
            shape_sum = 0.
            for j in -20:20
                shape_sum += shape_function(offset + j)
            end

            @test isapprox(shape_sum, 1)
        end
    end
end

@testset "ScatterChargeToGrid" begin
    num_cells = 8
    num_guard_cells = 2
    simulation_length = 1.
    grid = ParticleInCell.UniformGrid(num_cells, num_guard_cells, simulation_length)
    field = fill(0., ParticleInCell.total_cells(grid))

    num_macros = 100
    positions  = collect(range(0, stop=simulation_length, length=num_macros+1))[1:num_macros]

    shape_functions = [ParticleInCell.shape_1st_order,
                       ParticleInCell.shape_2nd_order]

    for shape_function in shape_functions
        field .= 0
        ParticleInCell.scatter_charge_to_grid!(positions, field,
            grid, 1., shape_function)

        # Check for charge conservation
        @test isapprox(sum(field[3:10]), num_macros)

        # Check that guard/periodic cells are equal to their counterparts
        # within the simulation cell
        @test field[1] == field[9]
        @test field[2] == field[10]
        @test field[3] == field[11]
        @test field[4] == field[12]
        @test field[5] == field[13]
    end




    @testset "Integration test" begin
        sim = ParticleInCell.Simulation(1.)

        num_cells = 32
        num_guard_cells = 5
        simulation_length = 1.
        g = ParticleInCell.UniformGrid(num_cells, num_guard_cells, simulation_length)

        f = ParticleInCell.Field(g)
        ParticleInCell.add_field!(sim, f, "f")

        charge = 1.
        mass = 1.
        particles_per_macro = 100
        num_macros = 320
        positions  = collect(range(0, stop=simulation_length, length=num_macros+1))[1:num_macros]
        velocities = similar(positions)
        forces     = similar(positions)
        s = ParticleInCell.Species(charge, mass, particles_per_macro, positions, velocities, forces)
        ParticleInCell.add_species!(sim, s, "s")

        sctg = ParticleInCell.ScatterChargeToGrid("s", "f")
        ParticleInCell.add_integration_step!(sim, sctg)

        ParticleInCell.setup!(sim)
        ParticleInCell.step!(sim)

        # Check charge conservation
        total_charge = num_macros * charge * particles_per_macro
        @test isapprox(sum(f.grid_values[1:num_cells]), total_charge)

        # Check guard/periodic cells
        for j in -num_guard_cells:num_guard_cells
            low_index = num_guard_cells + 1 + j
            high_index = num_cells + num_guard_cells + 1 + j

            @test f.grid_values[low_index] == f.grid_values[high_index]
        end
    end
end
