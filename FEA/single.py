import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

# Material properties
E = 210e9  # Young's modulus (Pa)
nu = 0.3   # Poisson's ratio
t = 0.01   # Plate thickness (m)

# Geometry and nodes
L = 1.0  # Plate length (m)
H = 1.0  # Plate height (m)

# Node coordinates (for two elements)
nodes = np.array([
    [0.0, 0.0],  # Node 1
    [L, 0.0],    # Node 2
    [L, H / 2],  # Node 3
    [0.0, H / 2],  # Node 4
    [L, H],  # Node 5
    [0.0, H]   # Node 6
])

# Connectivity (two elements)
elements = np.array([
    [0, 1, 2, 3],  # Element 1
    [3, 2, 4, 5]   # Element 2
])

# Degrees of freedom per node
dof_per_node = 3  # w, theta_x, theta_y

# Total degrees of freedom
num_nodes = nodes.shape[0]
num_dof = num_nodes * dof_per_node

# Shape function derivatives in natural coordinates (xi, eta)
def shape_function_derivatives(xi, eta):
    dN_dxi = np.array([
        [-0.25 * (1 - eta), 0.25 * (1 - eta), 0.25 * (1 + eta), -0.25 * (1 + eta)],
        [-0.25 * (1 - xi), -0.25 * (1 + xi), 0.25 * (1 + xi), 0.25 * (1 - xi)]
    ])
    return dN_dxi

# Compute element stiffness matrix
def compute_element_stiffness(xy, t, E, nu):
    D_b = (E * t**3) / (12 * (1 - nu**2)) * np.array([
        [1, nu, 0],
        [nu, 1, 0],
        [0, 0, (1 - nu) / 2]
    ])
    G = E / (2 * (1 + nu))
    kappa = 5 / 6
    D_s = kappa * G * t * np.eye(2)

    gauss_points = [(-np.sqrt(1/3), -np.sqrt(1/3)), 
                    (np.sqrt(1/3), -np.sqrt(1/3)),
                    (np.sqrt(1/3), np.sqrt(1/3)),
                    (-np.sqrt(1/3), np.sqrt(1/3))]
    weights = [1, 1, 1, 1]

    ke = np.zeros((12, 12))  # 4 nodes x 3 DOFs per node

    for gp, w in zip(gauss_points, weights):
        xi, eta = gp
        dN_dxi = shape_function_derivatives(xi, eta)
        J = np.dot(dN_dxi, xy)
        det_J = np.linalg.det(J)
        inv_J = np.linalg.inv(J)
        dN_dx = np.dot(inv_J, dN_dxi)

        B_b = np.zeros((3, 12))
        for i in range(4):
            B_b[0, i * 3 + 1] = dN_dx[0, i]
            B_b[1, i * 3 + 2] = dN_dx[1, i]
            B_b[2, i * 3 + 1] = dN_dx[1, i]
            B_b[2, i * 3 + 2] = dN_dx[0, i]

        B_s = np.zeros((2, 12))
        for i in range(4):
            B_s[0, i * 3] = dN_dx[0, i]
            B_s[1, i * 3] = dN_dx[1, i]
            B_s[0, i * 3 + 1] = dN_dx[1, i]
            B_s[1, i * 3 + 2] = dN_dx[0, i]

        ke += det_J * w * (np.dot(B_b.T, np.dot(D_b, B_b)) + np.dot(B_s.T, np.dot(D_s, B_s)))

    return ke

# Assemble global stiffness matrix
K = np.zeros((num_dof, num_dof))
for elem in elements:
    xy = nodes[elem, :]
    ke = compute_element_stiffness(xy, t, E, nu)

    element_dofs = []
    for node in elem:
        element_dofs.extend([node * dof_per_node, node * dof_per_node + 1, node * dof_per_node + 2])

    for i in range(12):
        for j in range(12):
            K[element_dofs[i], element_dofs[j]] += ke[i, j]

# Define the load vector (transverse load applied at Node 5)
load = np.zeros(num_dof)
load[12] = -1000  # Apply -1000 N transverse force at Node 5 (w DOF)

# Boundary conditions (fix all DOFs at Node 1 and Node 6)
fixed_dofs = [0, 1, 2, 15, 16, 17]  # Fix w, theta_x, theta_y at Node 1 and Node 6
free_dofs = np.setdiff1d(np.arange(num_dof), fixed_dofs)

# Reduce stiffness matrix and load vector
K_reduced = K[np.ix_(free_dofs, free_dofs)]
load_reduced = load[free_dofs]

# Solve for displacements
displacements = np.zeros(num_dof)
displacements[free_dofs] = np.linalg.solve(K_reduced, load_reduced)

# Output displacements
print("Nodal displacements (w, theta_x, theta_y):")
for i in range(num_nodes):
    print(f"Node {i + 1}: w = {displacements[i * dof_per_node]:.6e} m, "
          f"theta_x = {displacements[i * dof_per_node + 1]:.6e} rad, "
          f"theta_y = {displacements[i * dof_per_node + 2]:.6e} rad")

# Plot displacement contours
def plot_displacement_contours(nodes, displacements, dof_per_node, elements, title="Displacement Contours"):
    w_displacements = displacements[0::dof_per_node]
    triangulation = tri.Triangulation(nodes[:, 0], nodes[:, 1], elements[:, [0, 1, 2, 0, 2, 3]].reshape(-1, 3))
    plt.figure(figsize=(8, 6))
    contour = plt.tricontourf(triangulation, w_displacements, levels=12, cmap='viridis')
    plt.colorbar(contour, label='Transverse Displacement (w) [m]')
    for elem in elements:
        plt.plot(nodes[elem, 0], nodes[elem, 1], color='black', linestyle='--', linewidth=0.8)
    plt.title(title)
    plt.xlabel('X Coordinate [m]')
    plt.ylabel('Y Coordinate [m]')
    plt.axis('equal')
    plt.grid(True)
    plt.show()

plot_displacement_contours(nodes, displacements, dof_per_node, elements, title="Transverse Displacement Contours")
