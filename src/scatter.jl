mutable struct ScatterChargeToGrid{FT, IT, S} <: IntegrationStep
    species_name::String
    field_name::String

    species_index::IT
    field_index::IT

    shape_function::S
    charge::FT

    temp_field::Vector{FT}

    function ScatterChargeToGrid(species_name::String, field_name::String,
                                 shape_function::S=shape_1st_order, charge::FT=0.) where {FT, S}
        species_index = -1
        field_index = -1

        return new{FT, Int, S}(species_name, field_name, species_index,
                                    field_index, shape_function, charge)
    end
end

function setup!(step::ScatterChargeToGrid{FT}, sim::Simulation) where {FT}
    step.species_index = sim.species_dict[step.species_name]
    step.field_index = sim.fields_dict[step.field_name]

    if step.charge == zero(FT)
        step.charge = sim.species[step.species_index].macro_charge
    end

    step.temp_field = similar(sim.fields[step.field_index].grid_values)
end

function step!(step::ScatterChargeToGrid, sim::Simulation)
    positions = sim.species[step.species_index].positions
    field = sim.fields[step.field_index].grid_values
    grid = sim.fields[step.field_index].grid
    charge = step.charge

    step.temp_field .= 0

    scatter_charge_to_grid!(positions, step.temp_field, grid, charge, step.shape_function)

    field .+= step.temp_field
end

function scatter_charge_to_grid!(positions, field, grid, charge, shape_function)
    for x in positions, (index, grid_pos) in grid
        field[index] += charge * shape_function((grid_pos - x) / cell_length(grid, index)) / cell_length(grid, index)
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
