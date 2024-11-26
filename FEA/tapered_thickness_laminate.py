import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


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
t_min, t_max = 0.001, 0.005  # Thickness varies from 1 mm to 5 mm
A, B, D = compute_abd_matrix_with_taper(plies, y_coords, y_min, y_max, t_min, t_max)

# Apply weight loads
def apply_ply_weight(nodes, dof_per_node, load_vector, plies, gravity=9.81):
    total_thickness = sum(ply["thickness"] for ply in plies)
    average_density = sum(ply["density"] * ply["thickness"] for ply in plies) / total_thickness
    weight_per_unit_area = average_density * total_thickness * gravity
    num_nodes = nodes.shape[0]
    force_per_node = weight_per_unit_area * (Lx * Ly) / num_nodes
    for i in range(num_nodes):
        load_vector[i * dof_per_node] += force_per_node
    return load_vector

# Apply elliptical pressure loads
def apply_elliptic_pressure(nodes, dof_per_node, load_vector, pressure_max, a, b):
    for i, node in enumerate(nodes):
        x, y = node
        if (x**2 / a**2 + y**2 / b**2) <= 1.0:
            pressure = pressure_max * (1 - x**2 / a**2 - y**2 / b**2)
            load_vector[i * dof_per_node] += pressure
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
load = np.zeros(num_dof)

# Apply weight loads
load = apply_ply_weight(nodes, dof_per_node, load, plies)

# Apply elliptical pressure load
pressure_max = 5000.0  # Pascals
a = 0.5  # Semi-major axis (m)
b = 1.0  # Semi-minor axis (m)
load = apply_elliptic_pressure(nodes, dof_per_node, load, pressure_max, a, b)

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

def plot_displaced_body(nodes, displacements, dof_per_node, elements, scale_factor=1.0):
    """
    Plot the displaced body of the plate in 3D without axes and with a horizontal colorbar.
    
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
        scale_factor : float
            Factor to scale the displacements for better visualization.
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    # Extract transverse displacements (w) for each node
    w_displacements = displacements[0::dof_per_node] * scale_factor
    
    #print(w_displacements)

    # Create the 3D figure
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_facecolor('white')  # Set background color to white
    
    # Make the background enclosed by the axes transparent
    ax.patch.set_alpha(0)  # Set patch background transparency to 0
    
    # Plot the plate surface (displaced shape)
    for elem in elements:
        # Get nodal coordinates and displacements for the current element
        x = nodes[elem, 0]
        y = nodes[elem, 1]
        z = w_displacements[elem]  # Scaled displacements
        
        # Create the polygon for the element
        verts = [list(zip(x, y, z))]
        poly = Poly3DCollection(verts, alpha=0.8, edgecolor='k')
        poly.set_facecolor(plt.cm.viridis((z - z.min()) / (z.max() - z.min())))  # Map displacement to color
        ax.add_collection3d(poly)

    # Set axis limits explicitly to ensure the plot is not cut
    x_min, x_max = nodes[:, 0].min(), nodes[:, 0].max()
    y_min, y_max = nodes[:, 1].min(), nodes[:, 1].max()
    z_min, z_max = w_displacements.min(), w_displacements.max()
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_zlim(z_min, z_max)
    
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)
    
    # Set axis labels
    ax.set_xlabel('Chord (m)', fontsize=14, fontname="Times New Roman", labelpad=20)
    ax.set_ylabel('Span (m)', fontsize=14, fontname="Times New Roman", labelpad=20)
    
    # Turn off the z-axis
    ax.zaxis.set_visible(False)
    ax.zaxis.set_tick_params(labelsize=0, colors=(1, 1, 1, 0))  # Hide ticks
    ax.zaxis.line.set_color((1, 1, 1, 0))  # Hide the axis line
    ax.zaxis.label.set_color((1, 1, 1, 0))  # Hide the label

    
    #ax.set_zlabel('Transverse Displacement (w) [m]', fontsize=24, fontname="Times New Roman", labelpad=20)


    # Hide the axes and grid
    ax.grid(False)
    ax.axis("on")
    
    # Add colorbar horizontally
    mappable = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=w_displacements.min(), vmax=w_displacements.max()))
    mappable.set_array(w_displacements)
    cbar = fig.colorbar(mappable, ax=ax, orientation='horizontal', pad=0.1, fraction=0.02)
    cbar.set_label('Transverse Displacement (w) [m]', fontsize=12, fontname = "Times New Roman")

    # Set the aspect ratio and view angle
    ax.view_init(elev=30, azim=135)  # Isometric view
    ax.set_box_aspect([1, 2, 0.5])  # Adjust aspect ratio (X:Y:Z)
    
    F = plt.gcf()
    Size = F.get_size_inches()
    F.set_size_inches(Size[0]*1.5, Size[1]*1.5, forward=True) # Set forward to True to resize window along with plot in figure.
#F.set_size_inches(Size[0]*3.5, Size[1]*1.5, forward=True) # Set forward to True to resize window along with plot in figure.
    
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300
    
    # Show plot
    plt.tight_layout()
    plt.show()



plot_displaced_body(nodes, displacements, dof_per_node, elements, scale_factor=1.0)
