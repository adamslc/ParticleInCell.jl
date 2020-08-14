"""
Subtypes of `AbstractGrid` represent a one-dimensional, periodic spatial grid. By convention,
the grid is though to begin with cell 0 at the left edge of the simulation cell, and end with
cell `num_cells` at the right edge of the simulation, where cell 0 and cell `num_cells`
actually refer to the same physical location.

Subtypes should implement the following interface:
`cell_length(grid, index)`: returns the length of the cell at `index`.
`total_cells(grid)`: returns the total number of cells defined on the grid, including duplicate guard cells.
`iterate(grid)`: returns an iterator that yields a tuple of the previous two iterators.
"""
abstract type AbstractGrid{FT, IT} end

"""
    UniformGrid(num_cells, num_guard_cells, simulation_length, offset=0)

Stores information about a uniform spatial grid.
"""
struct UniformGrid{FT, IT} <: AbstractGrid{FT, IT}
    num_cells::IT
    num_guard_cells::IT

    simulation_length::FT
    cell_length::FT
    offset::FT

    function UniformGrid(num_cells::IT, num_guard_cells::IT,
                         simulation_length::FT, offset::FT=zero(FT)) where {FT, IT}
        cell_length = simulation_length / num_cells

        return new{FT, IT}(num_cells, num_guard_cells, simulation_length, cell_length, offset)
    end
end

function UniformGrid(; num_cells::IT, num_guard_cells=one(IT), simulation_length) where IT
    return UniformGrid(num_cells, num_guard_cells, simulation_length)
end

simulation_length(g::UniformGrid) = g.simulation_length
cell_length(g::UniformGrid, i=0) = g.cell_length
total_cells(g::UniformGrid) = g.num_cells + 2*g.num_guard_cells + 1

function Base.iterate(g::UniformGrid, state=(1, -g.num_guard_cells))
    if state[1] == g.num_cells + 2*g.num_guard_cells + 2
        return nothing
    end

    return ((state[1], state[2]*g.cell_length), (state[1]+1, state[2]+1))
end

Base.length(g::UniformGrid) = g.num_cells + 2*g.num_guard_cells + 1
