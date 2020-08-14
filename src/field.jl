"""
    Field(grid)

Stores the values of a field defined at the grid points defined by `grid`.
"""
mutable struct Field{FT, G <: AbstractGrid}
    grid::G
    grid_values::Vector{FT}

    function Field(grid::G) where {G <: AbstractGrid{FT}} where FT
        grid_values = Vector{FT}(undef, total_cells(grid))

        return new{FT, G}(grid, grid_values)
    end
end
