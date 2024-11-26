import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri


# Plate geometry
Lx = 1.0  # Chord length (m)
Ly = 2.0  # Span length (m)
nx = 10    # Number of elements along the chord (x-direction)
ny = 20    # Number of elements along the span (y-direction)


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

plies = [
    {"E1": 140e9, "E2": 10e9, "G12": 5e9, "nu12": 0.3, "thickness": 0.002, "angle": 0, "density": 1600},
    {"E1": 140e9, "E2": 10e9, "G12": 5e9, "nu12": 0.3, "thickness": 0.002, "angle": 90, "density": 1600},
    {"E1": 140e9, "E2": 10e9, "G12": 5e9, "nu12": 0.3, "thickness": 0.002, "angle": 45, "density": 1600},
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

def apply_elliptic_pressure(nodes, dof_per_node, load_vector, pressure_max, a, b):
    """
    Apply an elliptic pressure load on the top surface of the plate.

    Parameters:
        nodes : ndarray
            Nodal coordinates of the plate (Nx2 array where N is the number of nodes).
        dof_per_node : int
            Number of degrees of freedom per node (e.g., 3 for Reissner-Mindlin plate elements).
        load_vector : ndarray
            Global load vector to update (size: num_dof).
        pressure_max : float
            Maximum pressure at the center of the ellipse.
        a : float
            Semi-major axis of the ellipse (aligned with the x-direction).
        b : float
            Semi-minor axis of the ellipse (aligned with the y-direction).

    Returns:
        load_vector : ndarray
            Updated global load vector with the elliptic pressure applied.
    """
    # Loop over nodes to identify those on the top surface
    for i, node in enumerate(nodes):
        x, y = node
        # Elliptic pressure distribution
        if (x**2 / a**2 + y**2 / b**2) <= 1.0:  # Check if node is within the ellipse
            pressure = pressure_max * (1 - x**2 / a**2 - y**2 / b**2)  # Elliptic pressure formula
            load_vector[i * dof_per_node] += pressure  # Update transverse displacement DOF (w)
    
    return load_vector

def apply_ply_weight(nodes, dof_per_node, load_vector, plies, gravity=9.81):
    """
    Compute the weight of the composite ply and apply it as a distributed force on each node.

    Parameters:
        nodes : ndarray
            Nodal coordinates of the plate (Nx2 array where N is the number of nodes).
        dof_per_node : int
            Number of degrees of freedom per node (e.g., 3 for Reissner-Mindlin plate elements).
        load_vector : ndarray
            Global load vector to update (size: num_dof).
        plies : list of dict
            Ply properties, where each dict includes 'thickness' and 'density'.
        gravity : float
            Acceleration due to gravity (m/s^2). Default is 9.81 m/s^2.

    Returns:
        load_vector : ndarray
            Updated global load vector with the ply weight applied.
    """
    # Compute the total thickness and density of the laminate
    total_thickness = sum(ply["thickness"] for ply in plies)
    average_density = sum(ply["density"] * ply["thickness"] for ply in plies) / total_thickness

    # Compute the weight of the composite laminate per unit area
    weight_per_unit_area = average_density * total_thickness * gravity

    # Distribute the weight across all nodes
    num_nodes = nodes.shape[0]
    force_per_node = weight_per_unit_area * (Lx * Ly) / num_nodes

    # Apply the force in the transverse displacement direction (w) for each node
    for i in range(num_nodes):
        load_vector[i * dof_per_node] += force_per_node  # Add force to the w DOF
    
    return load_vector


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

for elem in elements:
    xy = nodes[elem, :]
    ke = compute_element_stiffness(xy, D, t=sum(ply["thickness"] for ply in plies),
                                   kappa=5/6, G13=5e9, G23=5e9,
                                   use_reissner_mindlin=use_reissner_mindlin)
    element_dofs = []
    for node in elem:
        element_dofs.extend([node * dof_per_node, node * dof_per_node + 1, node * dof_per_node + 2])
    for i in range(12):
        for j in range(12):
            K[element_dofs[i], element_dofs[j]] += ke[i, j]

# Apply boundary conditions and load
load = np.zeros(num_dof)
trailing_edge_nodes = [i for i, node in enumerate(nodes) if np.isclose(node[1], Ly)]

# Define parameters for the elliptic pressure
pressure_max = 5000.0  # Maximum pressure in Pascals
a = 0.5  # Semi-major axis (m)
b = 1.0  # Semi-minor axis (m)

#-------------------LOAD DEFINITION----------------------#

# Define gravity and material density for each ply
gravity = 9.81  # m/s^2
for ply in plies:
    ply["density"] = 1600  # Density in kg/m^3 for composite material (example value)
    
# Apply ply weight as distributed force
load = apply_ply_weight(nodes, dof_per_node, load, plies, gravity)

# Apply the elliptic pressure load
load = apply_elliptic_pressure(nodes, dof_per_node, load, pressure_max, a, b)

#--------------------------------------------------------#

fixed_dofs = []
for i, node in enumerate(nodes):
    if np.isclose(node[1], 0.0):
        fixed_dofs.extend([i * dof_per_node, i * dof_per_node + 1, i * dof_per_node + 2])

free_dofs = np.setdiff1d(np.arange(num_dof), fixed_dofs)
K_reduced = K[np.ix_(free_dofs, free_dofs)]
load_reduced = load[free_dofs]

# Solve for displacements
displacements = np.zeros(num_dof)
displacements[free_dofs] = np.linalg.solve(K_reduced, load_reduced)

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
    
    
    # Explicitly define axis limits
    plt.xlim(0, 1.0)  # Set x-axis limits (example: 0 to 1.0 meters)
    plt.ylim(0, 2.0)  # Set y-axis limits (example: 0 to 2.0 meters)

    
    # Remove the grid
    plt.grid(False)
    
    F = plt.gcf()
    Size = F.get_size_inches()
    F.set_size_inches(Size[0]*1.5, Size[1]*1.5, forward=True) # Set forward to True to resize window along with plot in figure.

    # Show plot
    plt.tight_layout()
    
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300
    #plt.savefig('Plots/MKIII_CLCD.png')
    plt.show()
    
    
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

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

#plot_displacement_contours(nodes, displacements, dof_per_node, elements, title="Displacement Contours")