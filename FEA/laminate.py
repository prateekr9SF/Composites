import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

# Plate geometry
Lx = 1.0  # Chord length (m)
Ly = 2.0  # Span length (m)
nx = 4    # Number of elements along the chord (x-direction)
ny = 4    # Number of elements along the span (y-direction)

# Generate nodes
x_coords = np.linspace(0, Lx, nx + 1)  # Chord-wise direction
y_coords = np.linspace(0, Ly, ny + 1)  # Span-wise direction
nodes = np.array([[x, y] for y in y_coords for x in x_coords])

# Define element connectivity
elements = []
for j in range(ny):
    for i in range(nx):
        n1 = j * (nx + 1) + i
        n2 = n1 + 1
        n3 = n1 + (nx + 1) + 1
        n4 = n1 + (nx + 1)
        elements.append([n1, n2, n3, n4])
elements = np.array(elements)

# Degrees of freedom per node
dof_per_node = 3  # w, theta_x, theta_y

# Total degrees of freedom
num_nodes = nodes.shape[0]
num_dof = num_nodes * dof_per_node

# Define ply properties
plies = [
    {"E1": 140e9, "E2": 10e9, "G12": 5e9, "nu12": 0.3, "thickness": 0.002, "angle": 0},
    {"E1": 140e9, "E2": 10e9, "G12": 5e9, "nu12": 0.3, "thickness": 0.002, "angle": 90},
    {"E1": 140e9, "E2": 10e9, "G12": 5e9, "nu12": 0.3, "thickness": 0.002, "angle": 45},
]

# Helper functions
def compute_ply_stiffness(ply):
    """Compute the transformed stiffness matrix Q* for a ply."""
    E1, E2, G12, nu12 = ply["E1"], ply["E2"], ply["G12"], ply["nu12"]
    angle = np.radians(ply["angle"])
    
    # Local stiffness matrix
    Q = np.array([
        [E1 / (1 - nu12 * E2 / E1), nu12 * E2 / (1 - nu12 * E2 / E1), 0],
        [nu12 * E2 / (1 - nu12 * E2 / E1), E2 / (1 - nu12 * E2 / E1), 0],
        [0, 0, G12]
    ])
    
    # Transformation matrix
    c = np.cos(angle)
    s = np.sin(angle)
    T = np.array([
        [c**2, s**2, 2 * c * s],
        [s**2, c**2, -2 * c * s],
        [-c * s, c * s, c**2 - s**2]
    ])
    
    return T @ Q @ T.T

def compute_abd_matrix(plies):
    """Compute the A, B, D stiffness matrices for the laminate."""
    A = np.zeros((3, 3))
    B = np.zeros((3, 3))
    D = np.zeros((3, 3))
    z_prev = -sum(ply["thickness"] for ply in plies) / 2  # Bottom surface

    for ply in plies:
        z_next = z_prev + ply["thickness"]  # Top surface of the current ply
        Q_star = compute_ply_stiffness(ply)
        A += Q_star * (z_next - z_prev)
        B += 0.5 * Q_star * (z_next**2 - z_prev**2)
        D += (1 / 3) * Q_star * (z_next**3 - z_prev**3)
        z_prev = z_next  # Move to the next layer

    return A, B, D

def compute_element_stiffness(xy, D_b, t, kappa, G13, G23):
    """Compute the element stiffness matrix for a composite laminate."""
    gauss_points = [(-np.sqrt(1/3), -np.sqrt(1/3)), 
                    (np.sqrt(1/3), -np.sqrt(1/3)),
                    (np.sqrt(1/3), np.sqrt(1/3)),
                    (-np.sqrt(1/3), np.sqrt(1/3))]
    weights = [1, 1, 1, 1]

    ke = np.zeros((12, 12))  # Element stiffness matrix

    for gp, w in zip(gauss_points, weights):
        xi, eta = gp

        # Derivatives of shape functions in natural coordinates
        dN_dxi = np.array([
            [-0.25 * (1 - eta), 0.25 * (1 - eta), 0.25 * (1 + eta), -0.25 * (1 + eta)],
            [-0.25 * (1 - xi), -0.25 * (1 + xi), 0.25 * (1 + xi), 0.25 * (1 - xi)]
        ])

        # Jacobian and its inverse
        J = np.dot(dN_dxi, xy)
        det_J = np.linalg.det(J)
        inv_J = np.linalg.inv(J)

        # Derivatives of shape functions in physical coordinates
        dN_dx = np.dot(inv_J, dN_dxi)

        # Bending strain-displacement matrix
        B_b = np.zeros((3, 12))
        for i in range(4):
            B_b[0, i * 3 + 1] = dN_dx[0, i]
            B_b[1, i * 3 + 2] = dN_dx[1, i]
            B_b[2, i * 3 + 1] = dN_dx[1, i]
            B_b[2, i * 3 + 2] = dN_dx[0, i]

        # Shear strain-displacement matrix
        B_s = np.zeros((2, 12))
        for i in range(4):
            B_s[0, i * 3] = dN_dx[0, i]
            B_s[1, i * 3] = dN_dx[1, i]
            B_s[0, i * 3 + 1] = dN_dx[1, i]
            B_s[1, i * 3 + 2] = dN_dx[0, i]

        # Bending stiffness contribution
        ke += det_J * w * (np.dot(B_b.T, np.dot(D_b, B_b)))

        # Shear stiffness contribution
        D_s = kappa * t * np.diag([G13, G23])  # Shear stiffness matrix
        ke += det_J * w * (np.dot(B_s.T, np.dot(D_s, B_s)))

    return ke

# Compute the ABD matrix
A, B, D = compute_abd_matrix(plies)

# Global stiffness matrix assembly
K = np.zeros((num_dof, num_dof))
for elem in elements:
    xy = nodes[elem, :]
    ke = compute_element_stiffness(xy, D, t=sum(ply["thickness"] for ply in plies),
                                   kappa=5/6, G13=5e9, G23=5e9)

    element_dofs = []
    for node in elem:
        element_dofs.extend([node * dof_per_node, node * dof_per_node + 1, node * dof_per_node + 2])

    for i in range(12):
        for j in range(12):
            K[element_dofs[i], element_dofs[j]] += ke[i, j]

# Define the load vector (distributed load along Y = Ly)
load = np.zeros(num_dof)
trailing_edge_nodes = [i for i, node in enumerate(nodes) if np.isclose(node[1], Ly)]
total_load = 1.0
load_per_node = total_load / len(trailing_edge_nodes)
for node in trailing_edge_nodes:
    load[node * dof_per_node] += load_per_node

# Boundary conditions (fix all DOFs at nodes along Y = 0)
fixed_dofs = []
for i, node in enumerate(nodes):
    if np.isclose(node[1], 0.0):  # Fix nodes where Y = 0
        fixed_dofs.extend([i * dof_per_node, i * dof_per_node + 1, i * dof_per_node + 2])

free_dofs = np.setdiff1d(np.arange(num_dof), fixed_dofs)

# Reduce stiffness matrix and load vector
K_reduced = K[np.ix_(free_dofs, free_dofs)]
load_reduced = load[free_dofs]

# Solve for displacements
displacements = np.zeros(num_dof)
displacements[free_dofs] = np.linalg.solve(K_reduced, load_reduced)

# Output displacements
#print("Nodal displacements (w, theta_x, theta_y):")
#for i in range(num_nodes):
#    print(f"Node {i + 1}: w = {displacements[i * dof_per_node]:.6e} m, "
#          f"theta_x = {displacements[i * dof_per_node + 1]:.6e} rad, "
#          f"theta_y = {displacements[i * dof_per_node + 2]:.6e} rad")


def plot_displacement_contours(nodes, displacements, dof_per_node, elements, title="Displacement Contours"):
    """
    Plot the displacement contours for the plate with an overlay of the FEA mesh.
    
    Parameters:
        nodes : ndarray
            Nodal coordinates of the plate (Nx2 array where N is the number of nodes).
        displacements : ndarray
            Global displacement vector (size: num_dof).
        dof_per_node : int
            Number of degrees of freedom per node (e.g., 3 for Reissner-Mindlin plate elements).
        elements : ndarray
            Element connectivity matrix (Ex4 array for quadrilateral elements).
        title : str
            Title of the plot.
    """
    import matplotlib.pyplot as plt
    import matplotlib.tri as tri

    # Extract transverse displacements (w) for each node
    w_displacements = displacements[0::dof_per_node]

    # Create triangulation for plotting
    triangulation = tri.Triangulation(nodes[:, 0], nodes[:, 1], elements[:, [0, 1, 2, 0, 2, 3]].reshape(-1, 3))

    # Plot the contour with the FEA mesh overlay
    plt.figure(figsize=(10, 8), dpi=300)
    contour = plt.tricontourf(triangulation, w_displacements, levels=12, cmap='RdBu')
    plt.colorbar(contour, label='Transverse Displacement (w) [m]')

    # Overlay the FEA mesh
    for elem in elements:
        plt.plot(nodes[elem, 0], nodes[elem, 1], color='black', linestyle='-', linewidth=0.5)

    # Formatting
    plt.title(title, fontsize=14)
    plt.xlabel('Chord-wise Direction (X) [m]', fontsize=12)
    plt.ylabel('Span-wise Direction (Y) [m]', fontsize=12)
    plt.axis('equal')
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    
    # Remove the grid
    plt.grid(False)

    # Show plot
    plt.tight_layout()
    plt.show()



plot_displacement_contours(nodes, displacements, dof_per_node, elements, title="Transverse Displacement Contours")
