using LinearAlgebra
using PyPlot

# Material properties for a single lamina of CFRP
function get_lamina_CFRP(E1, E2, G12, nu12)
    nu21 = nu12 * E2/ E1
    Q11 = E1 / (1 - nu12 * nu21)
    Q12 = nu12 * E2 / (1 - nu12 * nu21)
    Q22 = E2 / (1 - nu12 * nu21)
    Q66 = G12

    # Assemble the matrix
    Q = [Q11 Q12 0; Q12 Q22 0; 0 0 Q66]

    return Q

end

println("Passing")
