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


# Total thickness of laminate and stacking sequence
thickness = 0.003
stacking_sequence = [0, 45, -45, 90] .* 7

# Get stiffness matrix for a single lamina
lamina_stiffness = get_CFRP_stiffness()

# Compute the ABD matrix for the laminate
ABD_matrix = compute_abd_matrix(lamina_stiffness, thickness, stacking_sequence)

# Display the ABD matrix
println("ABD Matrix for the Carbon Fiber Reinforced Plastic (CFRP) laminate:\n")
println(ABD_matrix)




println("Passing")
