import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from output import *
from load import *


# Plate geometry
Lx = 1.0  # Chord length (m)
Ly = 2.0  # Span length (m)
nx = 10    # Number of elements along the chord (x-direction)
ny = 20    # Number of elements along the span (y-direction)

# Generate mesh
def generateMesh(Lx, Ly, nx, ny):
    x_coords = np.linspace(0, Lx, nx + 1)
    y_coords = np.linspace(0, Ly, ny + 1)
    nodes = np.array([[x, y] for y in y_coords for x in x_coords])
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

x_coords, y_coords, nodes, elements = generateMesh(Lx, Ly, nx, ny)

# Degrees of freedom per node
dof_per_node = 3  # w, theta_x, theta_y

# Define symmetric laminate
def define_symmetric_laminate(base_ply, stacking_sequence):
    laminate = []
    for angle in stacking_sequence:
        ply = base_ply.copy()
        ply["angle"] = angle
        laminate.append(ply)
    symmetric_part = laminate[::-1]
    laminate += symmetric_part
    return laminate

base_ply = {
    "E1": 140e9, "E2": 10e9, "G12": 5e9,
    "nu12": 0.3, "thickness": 0.000125, "density": 1600
}
stacking_sequence = [0, 45, -45, 90]
plies = define_symmetric_laminate(base_ply, stacking_sequence)

# Thickness taper function
def thickness_taper(y, y_min, y_max, t_min, t_max):
    return t_min + (t_max - t_min) * (y - y_min) / (y_max - y_min)

# Compute ply stiffness
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

# Compute A, B, D matrices with taper
def compute_abd_matrix_with_taper(plies, y_coords, y_min, y_max, t_min, t_max):
    A = np.zeros((3, 3))
    B = np.zeros((3, 3))
    D = np.zeros((3, 3))
    for y in y_coords:
        tapered_thickness = thickness_taper(y, y_min, y_max, t_min, t_max)
        z_prev = -tapered_thickness / 2
        for ply in plies:
            z_next = z_prev + ply["thickness"]
            Q_star = compute_ply_stiffness(ply)
            A += Q_star * (z_next - z_prev)
            B += 0.5 * Q_star * (z_next**2 - z_prev**2)
            D += (1 / 3) * Q_star * (z_next**3 - z_prev**3)
            z_prev = z_next
    return A, B, D

y_min, y_max = 0.0, Ly
t_min, t_max = 0.005, 0.005  # Thickness varies from 1 mm to 5 mm
A, B, D = compute_abd_matrix_with_taper(plies, y_coords, y_min, y_max, t_min, t_max)

def apply_weight_load(nodes, elements, dof_per_node, load_vector, plies, gravity=9.81):
    """
    Compute and apply the correct weight loads for each node, accounting for element areas.

    Parameters:
        nodes : ndarray
            Array of nodal coordinates (Nx2, where N is the number of nodes).
        elements : ndarray
            Element connectivity array (Ex4, where E is the number of elements).
        dof_per_node : int
            Degrees of freedom per node (e.g., 3 for Reissner-Mindlin elements).
        load_vector : ndarray
            Global load vector to be updated with weight loads (size: num_dof).
        plies : list of dict
            List of ply properties, where each dict includes 'thickness' and 'density'.
        gravity : float
            Acceleration due to gravity (m/s^2). Default is 9.81.

    Returns:
        load_vector : ndarray
            Updated global load vector with weight loads applied.
    """
    # Compute total laminate thickness and average density
    total_thickness = sum(ply["thickness"] for ply in plies)
    average_density = sum(ply["density"] * ply["thickness"] for ply in plies) / total_thickness

    # Weight per unit area
    weight_per_unit_area = total_thickness * average_density * gravity

    # Initialize node weights to zero
    node_weights = np.zeros(nodes.shape[0])

    # Loop through elements to distribute weight
    for elem in elements:
        # Extract node coordinates for the current element
        x_coords = nodes[elem, 0]
        y_coords = nodes[elem, 1]

        # Compute the area of the quadrilateral element using the shoelace formula
        element_area = 0.5 * abs(
            x_coords[0] * y_coords[1] + x_coords[1] * y_coords[2] +
            x_coords[2] * y_coords[3] + x_coords[3] * y_coords[0] -
            (y_coords[0] * x_coords[1] + y_coords[1] * x_coords[2] +
             y_coords[2] * x_coords[3] + y_coords[3] * x_coords[0])
        )

        # Total weight for this element
        element_weight = weight_per_unit_area * element_area
        
        # Distribute weight equally among the 4 nodes of the element
        for i in range(4):
            node_weights[elem[i]] += element_weight / 4

    # Apply the weight loads to the global load vector (w DOF only)
    for i, weight in enumerate(node_weights):
        load_vector[i * dof_per_node] -= weight

    return load_vector




def compute_element_stiffness(xy, D, t, kappa, G13, G23, use_reissner_mindlin):
    """
    Compute the element stiffness matrix for a quadrilateral element.

    Parameters:
        xy : ndarray
            Coordinates of the element nodes (4x2 array).
        D : ndarray
            Bending stiffness matrix (3x3 array).
        t : float
            Element thickness (m).
        kappa : float
            Shear correction factor.
        G13, G23 : float
            Shear moduli in the 1-3 and 2-3 planes (Pa).
        use_reissner_mindlin : bool
            Whether to include shear deformation effects (Reissner-Mindlin theory).

    Returns:
        ke : ndarray
            Element stiffness matrix (12x12 array).
    """
    # Gaussian quadrature points and weights for 2D integration
    gauss_points = [(-np.sqrt(1/3), -np.sqrt(1/3)), (np.sqrt(1/3), -np.sqrt(1/3)),
                    (np.sqrt(1/3), np.sqrt(1/3)), (-np.sqrt(1/3), np.sqrt(1/3))]
    weights = [1, 1, 1, 1]
    ke = np.zeros((12, 12))

    for gp, w in zip(gauss_points, weights):
        xi, eta = gp
        # Shape function derivatives in natural coordinates
        dN_dxi = np.array([
            [-0.25 * (1 - eta), 0.25 * (1 - eta), 0.25 * (1 + eta), -0.25 * (1 + eta)],
            [-0.25 * (1 - xi), -0.25 * (1 + xi), 0.25 * (1 + xi), 0.25 * (1 - xi)]
        ])
        # Jacobian matrix
        J = np.dot(dN_dxi, xy)
        det_J = np.linalg.det(J)
        inv_J = np.linalg.inv(J)
        dN_dx = np.dot(inv_J, dN_dxi)

        # Bending strain-displacement matrix
        B_b = np.zeros((3, 12))
        for i in range(4):
            B_b[0, i * 3 + 1] = dN_dx[0, i]
            B_b[1, i * 3 + 2] = dN_dx[1, i]
            B_b[2, i * 3 + 1] = dN_dx[1, i]
            B_b[2, i * 3 + 2] = dN_dx[0, i]

        # Integrate bending contribution
        ke += det_J * w * (np.dot(B_b.T, np.dot(D, B_b)))

        # Shear strain-displacement matrix (if Reissner-Mindlin is used)
        if use_reissner_mindlin:
            B_s = np.zeros((2, 12))
            for i in range(4):
                B_s[0, i * 3] = dN_dx[0, i]
                B_s[1, i * 3] = dN_dx[1, i]
            D_s = kappa * t * np.diag([G13, G23])
            ke += det_J * w * (np.dot(B_s.T, np.dot(D_s, B_s)))

    return ke


# Global stiffness matrix
num_nodes = nodes.shape[0]
num_dof = num_nodes * dof_per_node
K = np.zeros((num_dof, num_dof))
for elem in elements:
    xy = nodes[elem, :]
    t_elem = thickness_taper(np.mean(xy[:, 1]), y_min, y_max, t_min, t_max)
    ke = compute_element_stiffness(xy, D, t=t_elem, kappa=5/6, G13=5e9, G23=5e9, use_reissner_mindlin=True)
    element_dofs = []
    for node in elem:
        element_dofs.extend([node * dof_per_node, node * dof_per_node + 1, node * dof_per_node + 2])
    for i in range(12):
        for j in range(12):
            K[element_dofs[i], element_dofs[j]] += ke[i, j]

# Initialize load vector
gravity_load = np.zeros(num_dof)  # Gravity load vector
pressure_load = np.zeros(num_dof)  # Pressure load vector


# Apply weight loads
# Apply corrected weight loads to the global load vector
gravity_load = apply_weight_load(nodes, elements, dof_per_node, gravity_load, plies)


# Apply elliptical pressure load
pressure_max = 5000.0  # Pascals

# Apply linearly varying load
load_max = 5000.0  # Maximum load (e.g., N/m^2)
pressure_load = apply_spanwise_varying_load(nodes, dof_per_node, pressure_load, load_max, Ly)

# Combine gravity and pressure loads
load = gravity_load + pressure_load

# Apply boundary conditions
fixed_dofs = []
for i, node in enumerate(nodes):
    if np.isclose(node[1], 0.0):  # Fix trailing edge
        fixed_dofs.extend([i * dof_per_node, i * dof_per_node + 1, i * dof_per_node + 2])
free_dofs = np.setdiff1d(np.arange(num_dof), fixed_dofs)

# Solve the system
K_reduced = K[np.ix_(free_dofs, free_dofs)]
load_reduced = load[free_dofs]
displacements = np.zeros(num_dof)
displacements[free_dofs] = np.linalg.solve(K_reduced, load_reduced)


# Post-processing....
plot_displaced_body(nodes, displacements, dof_per_node, elements, scale_factor=1.0)

# Plot the weight load distribution with FEA mesh
plot_weight_loads(nodes, elements, gravity_load, dof_per_node)

plot_pressure_loads(nodes, elements, pressure_load, dof_per_node)