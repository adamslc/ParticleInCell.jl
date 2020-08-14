mutable struct ScatterChargeToGrid{IT, S} <: IntegrationStep
    species_name::String
    field_name::String

    species_index::IT
    field_index::IT

    shape_function::S
end

function ScatterChargeToGrid(species_name::String, field_name::String,
                             shape_function::S=shape_1st_order) where {S}
    species_index = -1
    field_index = -1

    return ScatterChargeToGrid{Int, S}(species_name, field_name, species_index,
                                field_index, shape_function)
end

function setup!(step::ScatterChargeToGrid, sim::Simulation)
    step.species_index = sim.species_dict[step.species_name]
    step.field_index = sim.fields_dict[step.field_name]
end

function step!(step::ScatterChargeToGrid, sim::Simulation)
    positions = sim.species[step.species_index].positions
    field = sim.fields[step.field_index].grid_values
    grid = sim.fields[step.field_index].grid
    charge = sim.species[step.species_index].macro_charge

    scatter_charge_to_grid!(positions, field, grid, charge, step.shape_function)
end

function scatter_charge_to_grid!(positions, field, grid, charge, shape_function)
    sleep(0.001)
    for x in positions, (index, grid_pos) in grid
        field[index] += charge * shape_function((x - grid_pos) / cell_length(grid, index))
    end

    # Correct for the guard/periodic cells
    for j in 1:grid.num_guard_cells
        field[j] += field[grid.num_cells + j]
        field[j + grid.num_cells] = field[j]
    end
    field[grid.num_guard_cells + 1] += field[grid.num_guard_cells + grid.num_cells + 1]
    field[grid.num_guard_cells + grid.num_cells + 1] = field[grid.num_guard_cells + 1]
    for j in 1:grid.num_guard_cells
        field[grid.num_guard_cells + 1 + j] += field[grid.num_guard_cells + grid.num_cells + 1 + j]
        field[grid.num_guard_cells + grid.num_cells + 1 + j] = field[grid.num_guard_cells + 1 + j]
    end
end
