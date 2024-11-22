import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri


# Plate geometry
Lx = 1.0  # Chord length (m)
Ly = 2.0  # Span length (m)
nx = 4    # Number of elements along the chord (x-direction)
ny = 4    # Number of elements along the span (y-direction)


def generateMesh(Lx, Ly, nx, ny):
    # Uniformly spaced x and y coordinates
    x_coords = np.linspace(0, Lx, nx + 1)
    y_coords = np.linspace(0, Ly, ny + 1)
    
    # Generate nodes
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
    
    return x_coords, y_coords, nodes, elements
    #end
    
x_coords, y_coords, nodes, elements = generateMesh(Lx, Ly, nx, ny)


# Degrees of freedom per node
dof_per_node = 3  # w, theta_x, theta_y (default for RM theory)

# Define ply properties
plies = [
    {"E1": 140e9, "E2": 10e9, "G12": 5e9, "nu12": 0.3, "thickness": 0.002, "angle": 0},
    {"E1": 140e9, "E2": 10e9, "G12": 5e9, "nu12": 0.3, "thickness": 0.002, "angle": 90},
    {"E1": 140e9, "E2": 10e9, "G12": 5e9, "nu12": 0.3, "thickness": 0.002, "angle": 45},
]


# Select theory
use_reissner_mindlin = True  # Set to False to use Classical Laminate Theory


# Helper functions
def compute_ply_stiffness(ply):
    """Compute the transformed stiffness matrix Q* for a ply."""
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

def compute_abd_matrix(plies):
    """Compute the A, B, D stiffness matrices for the laminate."""
    A = np.zeros((3, 3))
    B = np.zeros((3, 3))
    D = np.zeros((3, 3))
    z_prev = -sum(ply["thickness"] for ply in plies) / 2
    for ply in plies:
        z_next = z_prev + ply["thickness"]
        Q_star = compute_ply_stiffness(ply)
        A += Q_star * (z_next - z_prev)
        B += 0.5 * Q_star * (z_next**2 - z_prev**2)
        D += (1 / 3) * Q_star * (z_next**3 - z_prev**3)
        z_prev = z_next
    return A, B, D

def compute_element_stiffness(xy, D_b, t, kappa, G13, G23, use_reissner_mindlin):
    """Compute the element stiffness matrix."""
    gauss_points = [(-np.sqrt(1/3), -np.sqrt(1/3)), (np.sqrt(1/3), -np.sqrt(1/3)),
                    (np.sqrt(1/3), np.sqrt(1/3)), (-np.sqrt(1/3), np.sqrt(1/3))]
    weights = [1, 1, 1, 1]
    ke = np.zeros((12, 12))

    for gp, w in zip(gauss_points, weights):
        xi, eta = gp
        dN_dxi = np.array([
            [-0.25 * (1 - eta), 0.25 * (1 - eta), 0.25 * (1 + eta), -0.25 * (1 + eta)],
            [-0.25 * (1 - xi), -0.25 * (1 + xi), 0.25 * (1 + xi), 0.25 * (1 - xi)]
        ])
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

        ke += det_J * w * (np.dot(B_b.T, np.dot(D_b, B_b)))

        if use_reissner_mindlin:
            B_s = np.zeros((2, 12))
            for i in range(4):
                B_s[0, i * 3] = dN_dx[0, i]
                B_s[1, i * 3] = dN_dx[1, i]
            D_s = kappa * t * np.diag([G13, G23])
            ke += det_J * w * (np.dot(B_s.T, np.dot(D_s, B_s)))

    return ke

# Compute ABD matrix
A, B, D = compute_abd_matrix(plies)

# Global stiffness matrix
num_nodes = nodes.shape[0]
num_dof = num_nodes * dof_per_node
K = np.zeros((num_dof, num_dof))



