mutable struct GatherForcesFromGrid{IT, S} <: IntegrationStep
    species_name::String
    field_name::String

    species_index::IT
    field_index::IT

    shape_function::S
end

function GatherForcesFromGrid(species_name::String, field_name::String,
                             shape_function::S=shape_1st_order) where {S}
    species_index = -1
    field_index = -1

    return GatherForcesFromGrid{Int, S}(species_name, field_name, species_index,
                                field_index, shape_function)
end

function setup!(step::GatherForcesFromGrid, sim::Simulation)
    step.species_index = sim.species_dict[step.species_name]
    step.field_index = sim.fields_dict[step.field_name]
end


function step!(step::GatherForcesFromGrid, sim::Simulation)
    positions = sim.species[step.species_index].positions
    forces = sim.species[step.species_index].forces
    field = sim.fields[step.field_index].grid_values
    grid = sim.fields[step.field_index].grid
    charge = sim.species[step.species_index].macro_charge

    for (i, x) in enumerate(positions), (grid_index, grid_pos) in grid
        forces[i] += charge * field[grid_index] * step.shape_function((x - grid_pos) / cell_length(grid, grid_index))
    end

    return
end
