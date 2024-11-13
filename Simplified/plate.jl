using LinearAlgebra
using PyPlot

<<<<<<< HEAD:Simplified/plate.jl

include("stress.jl")
=======
#include("stress.jl")

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


>>>>>>> 1b09fbcfeeb66be909b054d4bbc985a0527b99c3:plate.jl

# Material properties for a single lamina of CFRP
function get_CFRP_stiffness()

    # Material properties
    E1 = 150e9  # Pa (fiber direction modulus)
    E2 = 10e9   # Pa (transverse direction modulus)
    G12 = 5e9   # Pa (in-plane shear modulus)
    nu12 = 0.3  # Poisson's ratio

    # Compute stiffness elements
    nu21 = nu12 * E2/ E1
    Q11 = E1 / (1 - nu12 * nu21)
    Q12 = nu12 * E2 / (1 - nu12 * nu21)
    Q22 = E2 / (1 - nu12 * nu21)
    Q66 = G12

    # Assemble matrix stiffness Q
    Q = [Q11 Q12 0; Q12 Q22 0; 0 0 Q66]

    return Q

end

# Rotation matix to transform stiffness for off-axis lamina
function get_transformation_matrix(theta)
    theta_rad = deg2rad(theta)
    c = cos(theta_rad)
    s = sin(theta_rad)
    
    # Assemble matrix T
    T = [c^2 s^2 2*c*s; s^2 c^2 -2*c*s; -c*s c*s c^2 - s^2]
    return T
end

# Transform stiffness marix to laminate coordinates
function transform_Q(Q, theta)
    T = get_transformation_matrix(theta)
    T_inv = inv(T)
    Q_Trans = T_inv * Q * T_inv
    return Q_Trans
end

# Compute the ABD matrix for a laminate
function compute_abd_matrix(lamina_properties, thickness, stacking_sequence)
    n = length(stacking_sequence)
    h = thickness / n  # Thickness of each lamina
    A = zeros(3, 3)
    B = zeros(3, 3)
    D = zeros(3, 3)
    z_prev = -thickness / 2  # Mid-plane offset
    for theta in stacking_sequence
        z_next = z_prev + h
        Q_bar = transform_Q(lamina_properties, theta)
        A .+= Q_bar * (z_next - z_prev)
        B .+= Q_bar * (z_next^2 - z_prev^2) / 2
        D .+= Q_bar * (z_next^3 - z_prev^3) / 3
        z_prev = z_next
    end
    return [A B; B D]
end

# Compute the ABD matrix for a laminate with varying thickness
function compute_abd_matrix_at_x(lamina_properties, thickness_at_x, stacking_sequence)
    n = length(stacking_sequence)
    h = thickness_at_x / n  # Thickness of each lamina at position x
    A = zeros(3, 3)
    B = zeros(3, 3)
    D = zeros(3, 3)
    z_prev = -thickness_at_x / 2  # Mid-plane offset
    for theta in stacking_sequence
        z_next = z_prev + h
        Q_bar = transform_Q(lamina_properties, theta)
        A .+= Q_bar * (z_next - z_prev)
        B .+= Q_bar * (z_next^2 - z_prev^2) / 2
        D .+= Q_bar * (z_next^3 - z_prev^3) / 3
        z_prev = z_next
    end
    return [A B; B D]
end


# Calculate deflection along the length of the tapered plate with elliptical load
function deflection_at_x(x, a, q_max, lamina_properties, stacking_sequence, thickness_max, thickness_min, density)
    # Get the local thickness at this x
    thickness_x = thickness_at_x(x, thickness_max, thickness_min, a)

    # Get the ABD matrix at this x (depends on local thickness)
    ABD_matrix_x = compute_abd_matrix_at_x(lamina_properties, thickness_x, stacking_sequence)

    # Extract D11 for the bending stiffness at this x
    D11_x = ABD_matrix_x[4, 1]  # Bending stiffness for x-direction

    # Get local load at this spanwise station
    q_x = q_total(x, a, q_max, thickness_max, thickness_min, density)

    # Compute the deflection at x considering the local thickness
    return (q_x * x^2 / (24 * D11_x)) * (6*a - x)
end

# Total load at any point x for elliptical lift distribution w/ tapered thickness
function q_total(x, a, q_max, thickness_max, thickness_min, density)
    # Compute thickness at x for the taper
    thickness_x = thickness_at_x(x, thickness_max, thickness_min, a)
    
    # Self-weight at this x
    q_self_weight_x = density * 9.81 * thickness_x
    
    # Total load (elliptical + self-weight)
    return q_max * sqrt(1 - (x/a)^2) + q_self_weight_x
end


# Function to compute bending moment considering varying self-weight and elliptical lift distribution
function bending_moment_tapered(x_vals, q_max, thickness_max, thickness_min, density, a)
    Mx_total = 0.0  # Initialize total bending moment at the fixed end (x = 0)

    for i in 1:length(x_vals)-1
        # Get x positions (current and next) for integration
        x = x_vals[i]
        x_next = x_vals[i + 1]
        
        # Compute thickness at x and x_next
        thickness_x = thickness_at_x(x, thickness_max, thickness_min, a)
        thickness_x_next = thickness_at_x(x_next, thickness_max, thickness_min, a)

        # Self-weight at x and x_next
        q_self_weight_x = density * 9.81 * thickness_x
        q_self_weight_x_next = density * 9.81 * thickness_x_next

        # Elliptical load at x and x_next
        q_elliptical_x = q_max * sqrt(1 - (x/a)^2)
        q_elliptical_x_next = q_max * sqrt(1 - (x_next/a)^2)

        # Total load (self-weight + elliptical) at x and x_next
        total_load_x = q_self_weight_x + q_elliptical_x
        total_load_x_next = q_self_weight_x_next + q_elliptical_x_next

        # Trapezoidal integration for bending moment (area under the load curve)
        Mx_total += 0.5 * (total_load_x + total_load_x_next) * (x_next - x) * (x_next^2 + x^2) / 2
    end

    return Mx_total
end

<<<<<<< HEAD:Simplified/plate.jl
=======


# Function to compute elliptical lift distribution and body weight distribution
>>>>>>> 1b09fbcfeeb66be909b054d4bbc985a0527b99c3:plate.jl
function plot_force_distributions(a, q_max, q_self_weight)
 
    # Set font to Times New Roman
    rc("font", family="Times New Roman")
    
    # Define span-wise locations along the plate
    x_vals = range(0, a, length=100)
    
    # Calculate the elliptical load and the body weight load at each point
    q_elliptical_vals = [q_max * sqrt(1 - (x/a)^2) for x in x_vals]
    q_body_weight_vals = [q_self_weight for _ in x_vals]  # Constant across the length
    
    # Plot the force distributions
    figure()
    plot(x_vals, q_elliptical_vals, label="Elliptical Load", color="b", linewidth=2)
    plot(x_vals, q_body_weight_vals, label="Body Weight Load", color="r", linestyle="--", linewidth=2)
    
    # Label the plot
    xlabel("Length along plate (m)")
    ylabel("Force Distribution (N/m²)")
    legend()
    
    # Save the plot as a PNG file with a DPI of 300
    savefig("Plots/force_distributions.png", dpi=300)
end

# Function to define thickness as a function of x (for a tapered plate)
function thickness_at_x(x, thickness_max, thickness_min, a)
    return thickness_max - (thickness_max - thickness_min) * (x / a)
end

# Calculate the mid-plane strain and curvature from applied moments
function compute_strain_curvature(A, B, D, M)
    # Assemble ABD matrix
    ABD_matrix = [A B; B D]
    
    # In-plane forces are zero for pure bending
    N = [0, 0, 0]

    # Combine in-plane forces and moments
    NM = vcat(N, M)

    # Solve for strain and curvature
    strain_curvature = ABD_matrix \ NM
    
    # Split into mid-plane strain and curvature components
    mid_plane_strain = strain_curvature[1:3]
    curvature = strain_curvature[4:6]
    
    return mid_plane_strain, curvature
end

<<<<<<< HEAD:Simplified/plate.jl
# Total thickness of laminate and stacking sequence
thickness = 0.003
=======

# Total thickness of laminate (max and min for tapered plate)
thickness_max = 0.003  # Max thickness at the root (x = 0)
thickness_min = 0.001  # Min thickness at the tip (x = a)

>>>>>>> 1b09fbcfeeb66be909b054d4bbc985a0527b99c3:plate.jl
stacking_sequence = [0, 45, -45, 90] .* 7

# Get stiffness matrix for a single lamina
lamina_stiffness = get_CFRP_stiffness()

# Compute the ABD matrix for the laminate
#ABD_matrix = compute_abd_matrix(lamina_stiffness, thickness, stacking_sequence)

# Plate dimensions (assuming a rectangular cantilevered plate)
a = 1.0  # meters (length of the cantilevered plate)

# Discretize the plate length
x_vals = range(0, a, length=100)  # 100 points along the length of the plate

# Apply uniform pressure (N/m^2)
q_max = 500  # N/m^2 (max load at the wing root)
density = 1600 # kg/m^3 for CFRP


# Compute the bending moment at the root, considering the varying thickness and self-weight
Mx = bending_moment_tapered(x_vals, q_max, thickness_max, thickness_min, density, a)

# No moment in the transverse direction for cantilevered loading
My = 0            

# External moments applied (Mx, My)
M = [Mx, My, 0]  # No twisting moment Mxy


# Compute the ABD matrix for the laminate at the root (thickest section)
ABD_matrix_root = compute_abd_matrix_at_x(lamina_stiffness, thickness_max, stacking_sequence)

# Extract A, B, and D matrices
A = ABD_matrix_root[1:3, 1:3]
B = ABD_matrix_root[1:3, 4:6]
D = ABD_matrix_root[4:6, 4:6]

# Compute mid-plane strain and curvature at the root
mid_plane_strain, curvature = compute_strain_curvature(A, B, D, M)



# Compute deflections considering the taper
w_vals = [deflection_at_x(x, a, q_max, lamina_stiffness, stacking_sequence, thickness_max, thickness_min, density) for x in x_vals]


stress_distribution = compute_stress_distribution_tapered(lamina_stiffness, thickness_max, thickness_min, a, stacking_sequence, mid_plane_strain, curvature, x_vals)

# Plot the stress distributions using the function from the included file

plot_all_stress_contours(x_vals, stress_distribution, thickness_max, thickness_min, a, stacking_sequence)

<<<<<<< HEAD:Simplified/plate.jl
plot_force_distributions(a, q_max, q_self_weight)
=======
#plot_stress_distributions(x_vals, stress_distributions, thickness, stacking_sequence)


>>>>>>> 1b09fbcfeeb66be909b054d4bbc985a0527b99c3:plate.jl

