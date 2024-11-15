using LinearAlgebra
using PyPlot

# Material property struct
struct Ply
    E1::Float64
    E2::Float64
    G12::Float64
    v12::Float64
    thickness::Float64
    angle::Float64
end



# Laminate stiffness matries
function compute_laminate_matricess(plies)
    A = zeros(3, 3)    # Extensional stiffness
    B = zeros(3, 3)    # Coupling stiffness
    D = zeros(3, 3)    # Bending stiffness
    # Total laminate thickness
    z = -sum(ply.thickness for ply in plies) / 2

    for ply in plies
        Q_local = Q_matrix(ply)
        Q_global = transform_Q(Q_local, ply.angle)
        z_next = z + ply.thickness

        # Populate stiffness matrices
        A += Q_global * ply.thickness
        B += Q_global * (z_next^2 - z^2) / 2
        D += Q_global * (z_next^3 - z^3) / 3
        z = z_next
    end
    return A, B, D
end

# Stiffness matrix Q for each ply
function Q_matrix(ply::Ply)
    v21 = ply.v12 * ply.E2 / ply.E1
    Q11 = ply.E1 / (1 - ply.v12 * v21)
    Q22 = ply.E2 / (1 - ply.v12 * v21)
    Q12 = ply.v12 * ply.E2 / (1 - ply.v12 * v21)
    Q66 = ply.G12

    return [
        Q11  Q12   0;
        Q12  Q22   0;
         0    0   Q66
    ]
end

# Transform Q to global coordinates for a given ply angle
function transform_Q(Q_local, angle)

    θ = deg2rad(angle)
    c = cos(θ)
    s = sin(θ)

    # Transformation matrix
    T = [c^2   s^2    2*c*s;
         s^2   c^2   -2*c*s;
        -c*s   c*s    c^2-s^2]

    return T' * Q_local * T
end

# Element stiffness matrix
function element_stiffness_matrix(A, B, D, node_coords)
    Ke = zeros(12, 12)
    gauss_points = [-1/sqrt(3), 1/sqrt(3)]
    weights = [1.0, 1.0]

    # Cycle thorough each gauss point
    for ξ in gauss_points
        for η in gauss_points
            N, dN_dξ, dN_dη = shape_functions_Q4(ξ, η)
            J = jacobian(node_coords, dN_dξ, dN_dη)
            detJ = det(J)
            J_inv = inv(J)
            B_matrix = strain_displacement_matrix_Q4(dN_dξ, dN_dη, J_inv)

            # Simplified combination: TODO: Check for actual implemtation
            stiffness_matrix = A + B + D

            # Element stiffness matrix
            Ke += B_matrix' * stiffness_matrix * B_matrix * detJ * weights[1] * weights[2]
        end
    end
    return Ke
end











function plot_mesh(nodes, elements)
    # Initialize a new figure
    figure()
    title("2D Mesh")
    xlabel("x")
    ylabel("y")
    axis("equal")
    
    # Loop over each element and plot its edges
    for element in elements
        x_coords = []
        y_coords = []
        
        # Collect x and y coordinates for each node in the element
        for node_idx in element
            x, y, _ = nodes[node_idx]  # Ignore z-coordinate for 2D plot
            push!(x_coords, x)
            push!(y_coords, y)
        end

        # Close the quadrilateral by appending the first node's coordinates again
        push!(x_coords, x_coords[1])
        push!(y_coords, y_coords[1])

        # Plot the element as a closed shape
        plot(x_coords, y_coords, "b-", lw=1.5, alpha=0.7)
    end

    # Display the plot
    savefig("mesh.png", dpi=300)
end

# Ply definition with stacking order
plies = [
    Ply(150e9, 10e9, 5e9, 0.3, 0.125, 0),
    Ply(150e9, 10e9, 5e9, 0.3, 0.125, 45),
    Ply(150e9, 10e9, 5e9, 0.3, 0.125, -45),
    Ply(150e9, 10e9, 5e9, 0.3, 0.125, 90)
]

# Generate mesh
nx, ny = 10, 5               # Elements in x and y directions
chord_length, span_length = 1.0, 0.5
dx = chord_length / nx
dy = span_length / ny

nodes = []
for j in 0:ny
    for i in 0:nx
        push!(nodes, (i * dx, j * dy, 0.0))  # (x, y, z=0)
    end
end

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


plot_mesh(nodes, elements)