using LinearAlgebra
using PyPlot

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

# Calculate deflection along the length of the plate with elliptical load
function deflection_at_x(x, a, q_max, D11, q_self_weight)
    # Get local load at this spanwise station
    q_x = q_total(x, a, q_max, q_self_weight)
    return (q_x * x^2 / (24 * D11)) * (6a - x)
end

# Total load at any point x for elliptical lift distribution
function q_total(x, a, q_max, q_self_weight)
    # Pressure load + body load
    return q_max * sqrt(1 - (x/a)^2) + q_self_weight
end

# Calculate the bending moment due to elliptical lift distribution at the cantilevered ene
function bending_moment_elliptical(a, q_max, q_self_weight)
    # Elliptical moment integral for maximum bending moment at the fixed end
    Mx = (q_max * a^2 / 4) + (q_self_weight * a^2 / 2)  # Combine self-weight and elliptical moment
    return Mx
end

# Total thickness of laminate and stacking sequence
thickness = 0.003
stacking_sequence = [0, 45, -45, 90] .* 7

# Get stiffness matrix for a single lamina
lamina_stiffness = get_CFRP_stiffness()

# Compute the ABD matrix for the laminate
ABD_matrix = compute_abd_matrix(lamina_stiffness, thickness, stacking_sequence)

# Plate dimensions (assuming a rectangular cantilevered plate)
a = 1.0  # meters (length of the cantilevered plate)

# Apply uniform pressure (N/m^2)
q_max = 500  # N/m^2 (max load at the wing root)
density = 1600 # kg/m^3 for CFRP

q_self_weight = density * 9.81 * thickness  # Load due to self-weight

# Calculate bending moments due to uniform load at the cantilevered end
Mx = bending_moment_elliptical(a, q_max, q_self_weight) 

# No moment in the transverse direction for cantilevered loading
My = 0            

# External moments applied (Mx, My)
M = [Mx, My, 0]  # No twisting moment Mxy


# In-plane forces are zero for bending due to uniform load
N = [0, 0, 0]

# Combine forces and moments
NM = vcat(N, M)  # Combine force and moment vectors
strain_curvature = ABD_matrix \ NM

D11 = ABD_matrix[4, 1]  # Bending stiffness for x-direction

# Extract mid-plane strains and curvatures
mid_plane_strain = strain_curvature[1:3]
curvature = strain_curvature[4:6]

# Discretize the plate length
x_vals = range(0, a, length=100)  # 100 points along the length of the plate
w_vals = [deflection_at_x(x, a, q_max, D11, q_self_weight) for x in x_vals]

# Plot deflection along the plate length using PyPlot
figure()
plot(x_vals, w_vals, color="b", linewidth=2)
xlabel("Length along plate (m)")
ylabel("Deflection (m)")
title("Deflection of Cantilevered Plate under Uniform Load")
# Save the plot as a PNG file
savefig("cantilevered_plate_deflection.png")

