mutable struct SmoothingFilter{FT, IT} <: IntegrationStep
    field_name::String
    field_index::IT

    temp_vector::Vector{FT}

    function SmoothingFilter(field_name)
        FT = Float64
        IT = Int

        new{FT, IT}(field_name, -1)
    end
end

function setup!(step::SmoothingFilter, sim::Simulation)
    step.field_index = sim.fields_dict[step.field_name]

    step.temp_vector = similar(sim.fields[step.field_index].grid_values)
end

function step!(step::SmoothingFilter, sim::Simulation)
    temp_vector = step.temp_vector
    field = sim.fields[step.field_index].grid_values
    grid = sim.fields[step.field_index].grid

    for i in 2:length(temp_vector) - 1
        temp_vector[i] = (field[i-1] + field[i] + field[i+1]) / 4
    end

    field .= temp_vector

    # Update guard/periodic cells
    for j in 1:grid.num_guard_cells
        field[j] = field[grid.num_cells + j]
    end
    field[grid.num_guard_cells + grid.num_cells + 1] = field[grid.num_guard_cells + 1]
    for j in 1:grid.num_guard_cells
        field[grid.num_guard_cells + grid.num_cells + 1 + j] = field[grid.num_guard_cells + 1 + j]
    end
end
