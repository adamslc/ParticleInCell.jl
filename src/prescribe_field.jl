"""
    PrescribeField(field_name::String, spacetime_function)

Adds the value of `spacetime_function` evaluated at each grid point to the
field with `field_name`.
"""
mutable struct PrescribeField{IT, STF} <: IntegrationStep
    field_name::String
    field_index::IT

    spacetime_function::STF

    function PrescribeField(field_name::String, spacetime_function::STF) where STF
        return new{Int, STF}(field_name, -1, spacetime_function)
    end
end

function setup!(step::PrescribeField, sim::Simulation)
    step.field_index = sim.fields_dict[step.field_name]
end

# TODO: When the code is generalized to 2 and 3 dimensions, this should be a
# generated function, so that each component of each node is individually
# passed to the spacetime_function. Some benchmarking has revealed that
# splatting incurs an unacceptable performance penalty.
function step!(step::PrescribeField, sim::Simulation)
    field = sim.fields[step.field_index]
    field_values = field.grid_values
    grid = field.grid

    for (index, grid_pos) in grid
        field_values[index] += step.spacetime_function(sim.time, grid_pos)
    end
end
