# Function to compute stress distributions at each point along the plate for a tapered plate
function compute_stress_distribution_tapered(lamina_properties, thickness_max, thickness_min, a, stacking_sequence, mid_plane_strain, curvature, x_vals)
    n = length(stacking_sequence)
    stress_distributions = []

    for x in x_vals
        # Compute local thickness at x
        thickness_x = thickness_at_x(x, thickness_max, thickness_min, a)

        # Compute the ABD matrix for this location based on local thickness
        ABD_matrix_x = compute_abd_matrix_at_x(lamina_properties, thickness_x, stacking_sequence)

        # Get mid-plane strain and curvature at this location (mid_plane_strain and curvature are constant)
        strain_curvature_x = mid_plane_strain .+ curvature .* x
        
        stress_at_x = []
        z_prev = -thickness_x / 2  # Start from the bottom of the laminate

        for i in 1:n
            # Get the z value for this ply
            z = z_prev + (thickness_x / n)

            # Compute strain at this z position (mid-plane strain + curvature * z)
            strain_z = strain_curvature_x .+ z .* curvature

            # Get the transformed stiffness matrix Q for this ply's orientation
            Q_bar = transform_Q(lamina_properties, stacking_sequence[i])

            # Compute stress using sigma = Q * epsilon
            stress_z = Q_bar * strain_z

            # Store the stress for this ply
            push!(stress_at_x, stress_z)

            z_prev = z
        end

        # Store the stress at this location x
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
# Converts stress values to MPa before plotting and adjusts for varying thickness due to taper
function plot_all_stress_contours(x_vals, stress_distributions, thickness_max, thickness_min, a, stacking_sequence)
    n = length(stacking_sequence)
    
    # Initialize figure for subplots with DPI set to 300
    figure(dpi=300)

    # Stress components to plot (σₓₓ, σᵧᵧ, σₓᵧ)
    components = ["xx", "yy", "xy"]
    titles = ["σₓₓ - Longitudinal Stress (MPa)", "σᵧᵧ - Transverse Stress (MPa)", "σₓᵧ - Shear Stress (MPa)"]

    # Predefine arrays to store z-values (through-thickness) for every x location
    z_vals_matrix = []
    z_length = 0

    # Precompute the maximum length of z-values (based on maximum thickness)
    for x in x_vals
        thickness_x = thickness_at_x(x, thickness_max, thickness_min, a)
        z_vals = range(-thickness_x/2, thickness_x/2, length=n+1)
        if length(z_vals) > z_length
            z_length = length(z_vals)
        end
    end

    # Create the mesh grid for z-values (through-thickness) and span-wise x-values
    z_grid = zeros(z_length, length(x_vals))
    for i in 1:length(x_vals)
        thickness_x = thickness_at_x(x_vals[i], thickness_max, thickness_min, a)
        z_vals = range(-thickness_x/2, thickness_x/2, length=n+1)
        z_grid[:, i] .= collect(z_vals)  # Assign z-values for each spanwise x location
    end

    for k in 1:3
        # Initialize a matrix to store interpolated stress values at each (x, z) point
        stress_matrix = zeros(n, length(x_vals))

        # Loop through spanwise locations (x) and laminae (through thickness)
        for i in 1:length(x_vals)
            # Get the local z-values for this x (due to taper)
            z_vals = z_grid[:, i]

            # Populate the stress_matrix by interpolating stress values along the z-axis for each x
            component_idx = k  # Choose which stress component to plot
            for j in 1:n  # For each lamina
                # Get the stress for the current lamina at this x position for the selected component
                stress_z = stress_distributions[i][j][component_idx]

                # Convert stress from Pa to MPa (1 MPa = 10^6 Pa)
                stress_z_mpa = stress_z / 1e6

                # Store the stress value at the corresponding (x, z) location
                stress_matrix[j, i] = stress_z_mpa
            end
        end

        # Determine min and max values for contours
        min_stress = minimum(stress_matrix)
        max_stress = maximum(stress_matrix)

        # Set contour levels (10 equally spaced levels between min and max)
        levels = range(min_stress, stop=max_stress, length=10)

        # Plot each stress component in its own subplot
        subplot(3, 1, k)  # Create a 3-row, 1-column grid of subplots
        contourf(x_vals, z_grid[:, 1], stress_matrix, levels=levels, cmap="bwr")

        # Add the colorbar and limit it to 4 ticks
        cbar = colorbar()
        cbar.set_ticks(range(min_stress, stop=max_stress, length=4))  # 4 tick values

        # Label the subplot
        xlabel("Span (m)")
        ylabel("T (m)")
        title(titles[k])
    end
    
    # Adjust layout and save the combined plot with DPI set to 300
    tight_layout()
    savefig("stress_contour_plots_tapered.png", dpi=300)
end
