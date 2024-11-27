import numpy as np


def apply_spanwise_varying_load(nodes, dof_per_node, load_vector, load_max, Ly):
    """
    Apply a load that varies linearly along the spanwise direction (y) and is constant across the chordwise direction (x).

    Parameters:
        nodes : ndarray
            Array of nodal coordinates (Nx2, where N is the number of nodes).
        dof_per_node : int
            Degrees of freedom per node (e.g., 3 for Reissner-Mindlin elements).
        load_vector : ndarray
            Global load vector to update.
        load_max : float
            Maximum load value at the root (y = 0).
        Ly : float
            Span length (y-direction).

    Returns:
        load_vector : ndarray
            Updated global load vector with the spanwise varying load applied.
    """
    for i, node in enumerate(nodes):
        x, y = node

        # Linearly varying load along the spanwise direction (y)
        load = load_max * (1 - y / Ly)

        # Apply the load to the w DOF of the node
        load_vector[i * dof_per_node] += load

    return load_vector
