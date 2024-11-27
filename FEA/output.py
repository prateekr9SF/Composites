import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import seaborn as sns
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
        scale_factor : float
            Factor to scale the displacements for better visualization.
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    import os

    # Create output directory
    os.makedirs('Plots', exist_ok=True)

    # Extract transverse displacements (w) for each node
    w_displacements = displacements[0::dof_per_node] * scale_factor

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
        # Normalize z globally for color mapping
        poly.set_facecolor(plt.cm.viridis((z - w_displacements.min()) / (w_displacements.max() - w_displacements.min())))
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

    # Add colorbar horizontally
    mappable = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=w_displacements.min(), vmax=w_displacements.max()))
    mappable.set_array(w_displacements)
    cbar = fig.colorbar(mappable, ax=ax, orientation='horizontal', pad=0.1, fraction=0.02)
    cbar.set_label('Transverse Displacement (w) [m]', fontsize=12, fontname="Times New Roman")

    # Hide the axes and grid
    ax.grid(False)
    ax.axis("on")
    
    # Set the view angle
    ax.view_init(elev=30, azim=135)
    ax.set_box_aspect([1, 2, 0.5])  # Adjust aspect ratio (X:Y:Z)
    
    F = plt.gcf()
    Size = F.get_size_inches()
    F.set_size_inches(Size[0]*1.5, Size[1]*1.5, forward=True) # Set forward to True to resize window along with plot in figure.

    
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300

    # Save the figure
    plt.tight_layout()
    plt.savefig('Plots/CFRP_Plate_Displaced.png', dpi=300)

    # Show plot
    #plt.show()

def plot_displacement_contours(nodes, displacements, dof_per_node, elements, scale_factor=1.0):
    """
    Plot a 2D surface plot of the displacement contours.

    Parameters:
        nodes : ndarray
            Nodal coordinates of the plate (Nx2 array where N is the number of nodes).
        displacements : ndarray
            Global displacement vector (size: num_dof).
        dof_per_node : int
            Number of degrees of freedom per node (e.g., 3 for Reissner-Mindlin plate elements).
        elements : ndarray
            Element connectivity matrix (Ex4 array for quadrilateral elements).
        scale_factor : float
            Factor to scale the displacements for better visualization.
    """
    import matplotlib.pyplot as plt
    import matplotlib.tri as tri
    
    # Convert quadrilateral elements to triangles
    triangles = []
    for elem in elements:
        # Split the quadrilateral into two triangles
        triangles.append([elem[0], elem[1], elem[2]])  # Triangle 1
        triangles.append([elem[0], elem[2], elem[3]])  # Triangle 2
    triangles = np.array(triangles)

    # Extract transverse displacements (w) for each node
    w_displacements = displacements[0::dof_per_node] * scale_factor

    # Create a triangular grid for plotting
    triangulation = tri.Triangulation(nodes[:, 0], nodes[:, 1], triangles)

    # Create the 2D contour plot
    plt.figure()
    contour = plt.tricontourf(triangulation, w_displacements, levels=10, cmap='viridis')

    # Add a colorbar
    cbar = plt.colorbar(contour, orientation='horizontal', pad=0.1, fraction = 0.03)
    cbar.set_label('Transverse Displacement (m)', fontsize=12, fontname="Times New Roman")
    cbar.ax.tick_params(labelsize=10)

    # Set plot titles and labels
    plt.xlabel('Chord (m)', fontsize=12, fontname="Times New Roman")
    plt.ylabel('Span (m)', fontsize=12, fontname="Times New Roman")
    plt.axis("equal")

    # Adjust axis limits
    #plt.xlim(nodes[:, 0].min(), nodes[:, 0].max())
    #plt.ylim(nodes[:, 1].min(), nodes[:, 1].max())

    plt.tight_layout()
    F = plt.gcf()
    Size = F.get_size_inches() 
    F.set_size_inches(Size[0]*1.5, Size[1]*1.5, forward=True) # Set forward to True to resize window along with plot in figure.

    plt.axis("off")
    plt.grid(False)
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300
    plt.savefig('Plots/CFRP_Plate_Weight.png')

    # Show the plot
    plt.show()


def plot_weight_loads(nodes, elements, load_vector, dof_per_node):
    """
    Plot the 2D surface contour of weight loads applied to each mesh node and overlay the FEA mesh.

    Parameters:
        nodes : ndarray
            Array of nodal coordinates (Nx2, where N is the number of nodes).
        elements : ndarray
            Element connectivity array (Ex4, where E is the number of elements).
        load_vector : ndarray
            Global load vector containing the weight loads (size: num_dof).
        dof_per_node : int
            Degrees of freedom per node (e.g., 3 for Reissner-Mindlin elements).
        title : str
            Title of the plot.
    """
    import matplotlib.pyplot as plt
    import matplotlib.tri as tri

    # Extract the transverse load (weight) at each node (w-direction)
    weight_loads = load_vector[0::dof_per_node]

    # Create a triangulation for the mesh
    triangulation = tri.Triangulation(nodes[:, 0], nodes[:, 1])
    
    sns_cmap = sns.color_palette("crest", as_cmap=True)

    # Plot the contour
    plt.figure()
    contour = plt.tricontourf(triangulation, weight_loads, levels=12, cmap=sns_cmap)
    #contour = plt.tricontourf(triangulation, weight_loads, levels=12)
    cbar = plt.colorbar(contour, orientation='horizontal', pad=0.1, fraction = 0.03)
    cbar.set_label("Weight (N)", fontsize=12, fontname="Times New Roman")
    
    plt.xlabel("Chord (m)", fontname="Times New Roman", fontsize=14)
    plt.ylabel("Span (m)",  fontname="Times New Roman", fontsize=14)
    plt.axis("equal")

    # Overlay the FEA mesh
    for elem in elements:
        x_coords = nodes[elem, 0]
        y_coords = nodes[elem, 1]
        plt.plot(np.append(x_coords, x_coords[0]), np.append(y_coords, y_coords[0]), 'k-', linewidth=0.5)

    plt.tight_layout()
    F = plt.gcf()
    Size = F.get_size_inches() 
    F.set_size_inches(Size[0]*1.5, Size[1]*1.5, forward=True) # Set forward to True to resize window along with plot in figure.

    plt.axis("off")
    plt.grid(False)
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300
    plt.savefig('Plots/CFRP_Plate_Weight.png')



def plot_pressure_loads(nodes, elements, load_vector, dof_per_node):
    """
    Plot the 2D surface contour of weight loads applied to each mesh node and overlay the FEA mesh.

    Parameters:
        nodes : ndarray
            Array of nodal coordinates (Nx2, where N is the number of nodes).
        elements : ndarray
            Element connectivity array (Ex4, where E is the number of elements).
        load_vector : ndarray
            Global load vector containing the pressure loads (size: num_dof).
        dof_per_node : int
            Degrees of freedom per node (e.g., 3 for Reissner-Mindlin elements).
        title : str
            Title of the plot.
    """
    import matplotlib.pyplot as plt
    import matplotlib.tri as tri

    # Extract the transverse load (weight) at each node (w-direction)
    weight_loads = load_vector[0::dof_per_node]

    # Create a triangulation for the mesh
    triangulation = tri.Triangulation(nodes[:, 0], nodes[:, 1])
    
    sns_cmap = sns.color_palette("flare", as_cmap=True)
    

    # Plot the contour
    plt.figure()
    contour = plt.tricontourf(triangulation, weight_loads, levels=12, cmap=sns_cmap)
    #contour = plt.tricontourf(triangulation, weight_loads, levels=12)
    cbar = plt.colorbar(contour, orientation='horizontal', pad=0.1, fraction = 0.03)
    cbar.set_label("Tractiion (N)", fontsize=12, fontname="Times New Roman")
    
    plt.xlabel("Chord (m)", fontname="Times New Roman", fontsize=14)
    plt.ylabel("Span (m)",  fontname="Times New Roman", fontsize=14)
    plt.axis("equal")

    # Overlay the FEA mesh
    for elem in elements:
        x_coords = nodes[elem, 0]
        y_coords = nodes[elem, 1]
        plt.plot(np.append(x_coords, x_coords[0]), np.append(y_coords, y_coords[0]), 'k-', linewidth=0.5)

    plt.tight_layout()
    F = plt.gcf()
    Size = F.get_size_inches() 
    F.set_size_inches(Size[0]*1.5, Size[1]*1.5, forward=True) # Set forward to True to resize window along with plot in figure.

    plt.axis("off")
    plt.grid(False)
    
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300
    plt.savefig('Plots/CFRP_Plate_Traction.png')
    #plt.show()
