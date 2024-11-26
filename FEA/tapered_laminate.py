import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Plate geometry with tapered planform
Lx_root = 1.0  # Chord length at the root (m)
Lx_tip = 0.5   # Chord length at the tip (m)
Ly = 2.0       # Span length (m)
nx = 10        # Number of elements along the chord (x-direction)
ny = 20        # Number of elements along the span (y-direction)

# Tapered thickness
t_min = 0.0001  # Minimum thickness (tip) in meters
t_max = 0.001   # Maximum thickness (root) in meters

# Functions for tapering
def compute_local_chord(y, Ly, Lx_root, Lx_tip):
    return Lx_tip + (Lx_root - Lx_tip) * (Ly - y) / Ly


def generateTaperedMesh(Lx_root, Lx_tip, Ly, nx, ny):
    y_coords = np.linspace(0, Ly, ny + 1)
    nodes = []
    elements = []
    for j, y in enumerate(y_coords):
        Lx_local = compute_local_chord(y, Ly, Lx_root, Lx_tip)
        x_coords = np.linspace(0, Lx_local, nx + 1)
        for x in x_coords:
            nodes.append([x, y])

    nodes = np.array(nodes)

    # Generate element connectivity
    for j in range(ny):
        for i in range(nx):
            n1 = j * (nx + 1) + i
            n2 = n1 + 1
            n3 = n1 + (nx + 1) + 1
            n4 = n1 + (nx + 1)
            elements.append([n1, n2, n3, n4])

    elements = np.array(elements)
    return nodes, elements


nodes, elements = generateTaperedMesh(Lx_root, Lx_tip, Ly, nx, ny)

def compute_tapered_thickness(x, y, Lx, Ly, t_min, t_max):
    return t_min + (t_max - t_min) * (Ly - y) / Ly

# Define symmetric laminate
base_ply = {
    "E1": 140e9,
    "E2": 10e9,
    "G12": 5e9,
    "nu12": 0.3,
    "thickness": 0.000125,
    "density": 1600,
}

stacking_sequence = [0, 45, -45, 90]

def define_symmetric_laminate(base_ply, stacking_sequence):
    laminate = []
    for angle in stacking_sequence:
        ply = base_ply.copy()
        ply["angle"] = angle
        laminate.append(ply)
    symmetric_part = laminate[::-1]
    laminate += symmetric_part
    return laminate


plies = define_symmetric_laminate(base_ply, stacking_sequence)

# Compute ABD matrix
def compute_ply_stiffness(ply):
    E1, E2, G12, nu12 = ply["E1"], ply["E2"], ply["G12"], ply["nu12"]
    angle = np.radians(ply["angle"])
    Q = np.array([
        [E1 / (1 - nu12 * E2 / E1), nu12 * E2 / (1 - nu12 * E2 / E1), 0],
        [nu12 * E2 / (1 - nu12 * E2 / E1), E2 / (1 - nu12 * E2 / E1), 0],
        [0, 0, G12]
    ])
    c, s = np.cos(angle), np.sin(angle)
    T = np.array([
        [c**2, s**2, 2 * c * s],
        [s**2, c**2, -2 * c * s],
        [-c * s, c * s, c**2 - s**2]
    ])
    return T @ Q @ T.T


def compute_abd_matrix_tapered(plies, nodes, elements, thickness_func):
    A = np.zeros((3, 3))
    B = np.zeros((3, 3))
    D = np.zeros((3, 3))

    for elem in elements:
        xy = nodes[elem, :]
        avg_thickness = np.mean([thickness_func(x, y, Lx_root, Ly, t_min, t_max) for x, y in xy])
        z_prev = -sum(ply["thickness"] for ply in plies) / 2 * avg_thickness
        for ply in plies:
            z_next = z_prev + ply["thickness"] * avg_thickness
            Q_star = compute_ply_stiffness(ply)
            A += Q_star * (z_next - z_prev)
            B += 0.5 * Q_star * (z_next**2 - z_prev**2)
            D += (1 / 3) * Q_star * (z_next**3 - z_prev**3)
            z_prev = z_next

    return A, B, D


A, B, D = compute_abd_matrix_tapered(plies, nodes, elements, compute_tapered_thickness)

# Global stiffness matrix
dof_per_node = 3
num_nodes = len(nodes)
num_dof = num_nodes * dof_per_node
K = np.zeros((num_dof, num_dof))

def compute_element_stiffness(xy, D):
    """
    Compute element stiffness matrix.
    """
    ke = np.zeros((12, 12))
    for i in range(4):
        ke[i * 3:(i + 1) * 3, i * 3:(i + 1) * 3] = D
    return ke

for elem in elements:
    xy = nodes[elem, :]
    ke = compute_element_stiffness(xy, D)
    element_dofs = []
    for node in elem:
        element_dofs.extend([node * dof_per_node, node * dof_per_node + 1, node * dof_per_node + 2])
    for i in range(12):
        for j in range(12):
            K[element_dofs[i], element_dofs[j]] += ke[i, j]

# Apply boundary conditions
fixed_dofs = []
for i, node in enumerate(nodes):
    if np.isclose(node[1], 0.0):
        fixed_dofs.extend([i * dof_per_node, i * dof_per_node + 1, i * dof_per_node + 2])

free_dofs = np.setdiff1d(np.arange(num_dof), fixed_dofs)

# Example load
load = np.zeros(num_dof)
load[-dof_per_node] = 1000.0  # Apply load at the last node in the w direction

# Solve system
K_reduced = K[np.ix_(free_dofs, free_dofs)]
load_reduced = load[free_dofs]
displacements = np.zeros(num_dof)
displacements[free_dofs] = np.linalg.solve(K_reduced, load_reduced)

# Plot results
def plot_displacement_contours(nodes, displacements, dof_per_node, elements):
    w_displacements = displacements[0::dof_per_node]
    plt.tricontourf(nodes[:, 0], nodes[:, 1], w_displacements, cmap="RdBu")
    plt.colorbar(label="Transverse Displacement (m)")
    plt.title("Displacement Contours")
    plt.xlabel("X (m)")
    plt.ylabel("Y (m)")
    plt.axis("equal")
    plt.show()


plot_displacement_contours(nodes, displacements, dof_per_node, elements)
