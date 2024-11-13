using JuAFEM
using LinearAlgebra

# Function to generate a rectangular mesh with Q4 elements
# Function to generate a rectangular mesh with Q4 elements
function generate_rectangular_mesh(chord_length, span_length, nx, ny)
    # Define node positions
    dx = chord_length / nx
    dy = span_length / ny
    nodes = Vector{Vector{Float64}}((nx + 1) * (ny + 1))
    
    # Generate node positions in x (chord-wise), y (span-wise), and z=0 (vertical)
    count = 1
    for j in 0:ny
        for i in 0:nx
            nodes[count] = [i * dx, j * dy, 0.0]  # z is initially zero
            count += 1
        end
    end
    
    # Define Q4 elements with connectivity (node indices)
    elements = []
    for j in 0:ny-1
        for i in 0:nx-1
            n1 = i + j * (nx + 1) + 1
            n2 = n1 + 1
            n3 = n1 + nx + 2
            n4 = n1 + nx + 1
            push!(elements, (n1, n2, n3, n4))
        end
    end
    
    # Create the grid from nodes and elements
    grid = generate_grid(elements, nodes)
    return grid
end


# Plate dimensions
chord_length = 1.0
span_length = 0.5
nx = 10
ny = 5

mesh = generate_rectangular_mesh(chord_length, span_length, nx, ny)

println("Passing package extraction!")