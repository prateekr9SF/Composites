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


function apply_uniform_force(f, nodes, chord_length, total_magnitude)
    # Find nodes at the free end (x = chord_length)
    free_nodes = [n for n in 1:length(nodes) if nodes[n][1] ≈ chord_length]

    # Compute force per node
    num_free_nodes = length(free_nodes)

    if num_free_nodes == 0
        error("No nodes found at the free end to apply forces.")
    end

    force_per_node = total_magnitude / num_free_nodes

    # Apply forces in the z-direction
    for n in free_nodes
        f[(n - 1) * 3 + 3] += force_per_node  # Add force to z-direction DOF
    end

    println("Applied a total force of $total_magnitude N in z-direction, distributed uniformly.")

    return f
end


# Laminate stiffness matries
function compute_laminate_matrices(plies)
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

    # Construct the block stiffness matrix
    C = [A B;
         B D]

    # Cycle thorough each gauss point
    for ξ in gauss_points
        for η in gauss_points
            # Get shape function for his gauss point
            N, dN_dξ, dN_dη = shape_functions_Q4(ξ, η)

            # Get jacobian for this gauss point
            J = jacobian(node_coords, dN_dξ, dN_dη)

            detJ = det(J)
            J_inv = inv(J)

            # Get strain-displacement matrix
            B_matrix = strain_displacement_matrix_Q4(dN_dξ, dN_dη, J_inv)

            # Extend B to include bending components
            B_extended = vcat(B_matrix, B_matrix)  # Duplicate rows for bending

            # Compute stiffness matrix contribution at this Gauss point
            Ke += B_extended' * C * B_extended * detJ * weights[1] * weights[2]

            # Simplified combination: TODO: Check for actual implemtation
            # Done: Look above
            #stiffness_matrix = A + B + D

            # Element stiffness matrix
            #Ke += B_matrix' * stiffness_matrix * B_matrix * detJ * weights[1] * weights[2]
        end
    end
    return Ke
end

# Shape functions and their derivatives for Q4 element
function shape_functions_Q4(ξ, η)
    N = [
        (1 - ξ) * (1 - η) / 4;
        (1 + ξ) * (1 - η) / 4;
        (1 + ξ) * (1 + η) / 4;
        (1 - ξ) * (1 + η) / 4
    ]
    dN_dξ = [
        -(1 - η) / 4, (1 - η) / 4, (1 + η) / 4, -(1 + η) / 4
    ]

    dN_dη = [
        -(1 - ξ) / 4, -(1 + ξ) / 4, (1 + ξ) / 4, (1 - ξ) / 4
    ]

    return N, dN_dξ, dN_dη
end

# Jacobian for Q4 element
function jacobian(node_coords, dN_dξ, dN_dη)
    J = zeros(2, 2)
    for i in 1:4
        x, y = node_coords[i][1], node_coords[i][2]
        J[1, 1] += dN_dξ[i] * x
        J[1, 2] += dN_dξ[i] * y
        J[2, 1] += dN_dη[i] * x
        J[2, 2] += dN_dη[i] * y
    end
    return J
end

# Strain-displacement matrix B
function strain_displacement_matrix_Q4(dN_dξ, dN_dη, J_inv)
    B = zeros(3, 12)
    for i in 1:4
        dN_dx = J_inv[1, 1] * dN_dξ[i] + J_inv[1, 2] * dN_dη[i]
        dN_dy = J_inv[2, 1] * dN_dξ[i] + J_inv[2, 2] * dN_dη[i]
        B[1, 3*(i-1) + 1] = dN_dx
        B[2, 3*(i-1) + 2] = dN_dy
        B[3, 3*(i-1) + 1] = dN_dy
        B[3, 3*(i-1) + 2] = dN_dx
    end
    return B
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

#### FEA CALCULATION ####

# Assemble global stiffness matrix
num_nodes = length(nodes)
num_dofs = 3 * num_nodes

K = zeros(num_dofs, num_dofs)

# Stiffness matrix for the entire laminate
A, B, D = compute_laminate_matrices(plies)

for elem in elements
    node_coords = [nodes[n] for n in elem]
    Ke = element_stiffness_matrix(A, B, D, node_coords)

    # Get indices and assemble global stiffness matrix
    for i in 1:4
        for j in 1:4
            for d1 in 1:3
                for d2 in 1:3
                    global_i = (elem[i] - 1) * 3 + d1
                    global_j = (elem[j] - 1) * 3 + d2
                    K[global_i, global_j] += Ke[(i-1)*3 + d1, (j-1)*3 + d2]
                end
            end
        end
    end
end

# Boundary conditions
f = zeros(num_dofs)


# Apply uniform force on the free end
total_force = 1.0  # Total force in N
f = apply_uniform_force(f, nodes, chord_length, total_force)

fixed_nodes = [n for n in 1:num_nodes if nodes[n][1] ≈ 0.0]


# Reduce global stiffness matrix
for n in fixed_nodes
    # Get dof indices
    dofs = (n - 1) * 3 .+ (1:3)

    # Loop over indices and enforece zero dispacement
    for dof in dofs
        K[dof, :] .= 0.0
        K[:, dof] .= 0.0
        K[dof, dof] = 1.0
        f[dof] = 0.0
    end
end

# Solve system
u = K \ f





