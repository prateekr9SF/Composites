
# Function to compute stress distributions at each point along the plate
function compute_stress_distribution(lamina_properties, thickness, stacking_sequence, mid_plane_strain, curvature, x_vals)
    n = length(stacking_sequence)
    h = thickness / n  # Thickness of each lamina
    z_vals = range(-thickness/2, thickness/2, length=n+1)  # z-values through the laminate thickness
    
    # Store stress distributions for each lamina at each x
    stress_distributions = []
    
    for x in x_vals
        stress_at_x = []
        
        z_prev = -thickness / 2  # Start from the bottom of the laminate
        for i in 1:n
            # Get the z value for this ply
            z = (z_vals[i] + z_vals[i+1]) / 2  # Midpoint of ply thickness
            
            # Compute the strain at this z
            strain_z = mid_plane_strain .+ z .* curvature
            
            # Get the transformed stiffness matrix Q for this ply's orientation
            Q_bar = transform_Q(lamina_properties, stacking_sequence[i])
            
            # Compute the stress using sigma = Q * epsilon
            stress_z = Q_bar * strain_z
            
            # Store the stress for this ply
            push!(stress_at_x, stress_z)
        end
        
        # Store the stress for this location x
        push!(stress_distributions, stress_at_x)
    end
    
    return stress_distributions
end


# Function to compute stress distributions at each point along the plate
function compute_stress_distribution_all(lamina_properties, thickness, stacking_sequence, mid_plane_strain, curvature, x_vals)
    n = length(stacking_sequence)
    h = thickness / n  # Thickness of each lamina
    z_vals = range(-thickness/2, thickness/2, length=n+1)  # z-values through the laminate thickness
    
    # Store stress distributions for each lamina at each x
    stress_distributions = []
    
    for x in x_vals
        stress_at_x = []
        
        z_prev = -thickness / 2  # Start from the bottom of the laminate
        for i in 1:n
            # Get the z value for this ply
            z = (z_vals[i] + z_vals[i+1]) / 2  # Midpoint of ply thickness
            
            # Compute the strain at this z
            strain_z = mid_plane_strain .+ z .* curvature
            
            # Get the transformed stiffness matrix Q for this ply's orientation
            Q_bar = transform_Q(lamina_properties, stacking_sequence[i])
            
            # Compute the stress using sigma = Q * epsilon
            stress_z = Q_bar * strain_z
            
            # Store the stress for this ply (all components: σₓₓ, σᵧᵧ, σₓᵧ)
            push!(stress_at_x, stress_z)
        end
        
        # Store the stress for this location x
        push!(stress_distributions, stress_at_x)
    end
    
    return stress_distributions
end


# Helper function to create a meshgrid
function meshgrid(x, y)
    X = repeat(x', length(y), 1)
    Y = repeat(y, 1, length(x))
    return X, Y
end


# Function to create contour plots of stress distributions (σₓₓ, σᵧᵧ, σₓᵧ) as subplots
function plot_all_stress_contours(x_vals, stress_distributions, thickness, stacking_sequence)
    n = length(stacking_sequence)
    h = thickness / n  # Thickness of each lamina
    z_vals = range(-thickness/2, thickness/2, length=n+1)  # z-values through the laminate thickness

    # Create a mesh grid for x (length) and z (through-thickness) values
    X, Z = meshgrid(x_vals, z_vals)

    # Initialize figure for subplots
    figure()

    # Stress components to plot (σₓₓ, σᵧᵧ, σₓᵧ)
    components = ["xx", "yy", "xy"]
    titles = ["σₓₓ - Longitudinal Stress", "σᵧᵧ - Transverse Stress", "σₓᵧ - Shear Stress"]

    for k in 1:3
        # Initialize a matrix to store interpolated stress values at each (x, z) point
        stress_matrix = zeros(length(z_vals), length(x_vals))

        # Populate the stress_matrix by interpolating stress values along the z-axis for each x
        component_idx = k  # Choose which stress component to plot
        for i in 1:length(x_vals)
            for j in 1:n  # For each lamina
                # Get the stress for the current lamina at this x position for the selected component
                stress_z = stress_distributions[i][j][component_idx]

                # Fill in the corresponding z range for this lamina
                stress_matrix[j, i] = stress_z
            end
        end

        # Plot each stress component in its own subplot
        subplot(3, 1, k)  # Create a 3-row, 1-column grid of subplots
        contourf(X, Z, stress_matrix, cmap="viridis")
        colorbar()

        # Label the subplot
        xlabel("Length along plate (m)")
        ylabel("Through-thickness position (m)")
        title(titles[k])
    end
    
    # Adjust layout and save the combined plot
    tight_layout()
    savefig("stress_contour_plots.png")
    println("Combined stress contour plots saved as stress_contour_plots.png")
end
