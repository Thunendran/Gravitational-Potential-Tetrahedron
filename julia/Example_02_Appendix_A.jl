# 1. Include the module file.
include("GP_TetrahedronBigFloat_Optimized.jl")

# 2. Bring the module's exported functions
using .GP_TetrahedronBigFloat_Optimized
using Printf # For formatted printing
using StaticArrays
#
# ----- New Functions for Gradient Calculation -----
#

"""
    compute_gradient(phi_func, x, y, z, h)

Computes the gradient of a scalar function `phi_func` at point (x, y, z)
using the central difference method. The result is the negative gradient,
representing the force field.
"""
function compute_gradient(phi_func, x::BigFloat, y::BigFloat, z::BigFloat, h::BigFloat)
    gx = (phi_func(x + h, y, z) - phi_func(x - h, y, z)) / (BigFloat(2.0) * h)
    gy = (phi_func(x, y + h, z) - phi_func(x, y - h, z)) / (BigFloat(2.0) * h)
    gz = (phi_func(x, y, z + h) - phi_func(x, y, z - h)) / (BigFloat(2.0) * h)
    return -gx, -gy, -gz
end

#
# ----- Main Execution Logic -----
#

function main()
    # Set the precision for BigFloat calculations.
    setprecision(256)

    # Constants and vertex setup
    G = BigFloat(1.0)
    sigma = BigFloat(1.0)
    h_grad = BigFloat(1e-20) # A small step for gradient calculation
    epsilon = BigFloat(1e-10)
    eps_symbol = "ε"

    # Define tetrahedron vertices as BigFloat vectors
    #vA_vec = BigFloat.([1.0, 1.0, 1.0])
    #vB_vec = BigFloat.([2.0, 1.0, 1.0])
    #vC_vec = BigFloat.([1.0, 2.0, 1.0])
    #vD_vec = BigFloat.([1.0, 1.0, 2.0])

    vA_vec = BigFloat.([0.0, 0.0, 0.0])
    vB_vec = BigFloat.([1.0, 0.0, 0.0])
    vC_vec = BigFloat.([0.0, 1.0, 0.0])
    vD_vec = BigFloat.([0.0, 0.0, 1.0])

    # 
    svA = SVector{3}(vA_vec)
    svB = SVector{3}(vB_vec)
    svC = SVector{3}(vC_vec)
    svD = SVector{3}(vD_vec)

    # Update the wrapper to use the SVectors
    # --- build the tetra object once (from the optimized module) ---
    tet = GP_TetrahedronBigFloat_Optimized.TetrahedronBigFloat(svA, svB, svC, svD)

    # --- call the optimized potential function (no other changes) ---
    potential_phi(x, y, z) = GP_TetrahedronBigFloat_Optimized.compute_potential(
        tet, SVector{3}(x, y, z); G=G, sigma=sigma
    )


    # Define labeled test points (original and perturbed)
    labeled_test_points = [
    ("Along AB 1", [BigFloat.([2, 0, 0]), BigFloat.([2 + epsilon, 0, 0])], ["2, 0, 0", "2+$(eps_symbol), 0, 0"]),
    ("Along AB 2", [BigFloat.([-1, 0, 0]), BigFloat.([-1 + epsilon, 0, 0])], ["-1, 0, 0", "-1+$(eps_symbol), 0, 0"]),
    
    ("Along BC 1", [BigFloat.([2, -1, 0]), BigFloat.([2 + epsilon, -1 - epsilon, 0])], ["2, -1, 0", "2+$(eps_symbol), -1-$(eps_symbol), 0"]),
    ("Along BC 2", [BigFloat.([-1, 2, 0]), BigFloat.([-1 - epsilon, 2 + epsilon, 0])], ["-1, 2, 0", "-1-$(eps_symbol), 2+$(eps_symbol), 0"]),
    
    ("Along AC 1", [BigFloat.([0, -1, 0]), BigFloat.([0, -1 - epsilon, 0])], ["0, -1, 0", "0, -1-$(eps_symbol), 0"]),
    ("Along AC 2", [BigFloat.([0, 2, 0]), BigFloat.([0, 2 + epsilon, 0])], ["0, 2, 0", "0, 2+$(eps_symbol), 0"]),
    
    ("Along BD 1", [BigFloat.([2, 0, -1]), BigFloat.([2 + epsilon, 0, -1 - epsilon])], ["2, 0, -1", "2+$(eps_symbol), 0, -1-$(eps_symbol)"]),
    ("Along BD 2", [BigFloat.([-1, 0, 2]), BigFloat.([-1 - epsilon, 0, 2 + epsilon])], ["-1, 0, 2", "-1-$(eps_symbol), 0, 2+$(eps_symbol)"]),
    
    ("Along CD 1", [BigFloat.([0, 2, -1]), BigFloat.([0, 2 + epsilon, -1 - epsilon])], ["0, 2, -1", "0, 2+$(eps_symbol), -1-$(eps_symbol)"]),
    ("Along CD 2", [BigFloat.([0, -1, 2]), BigFloat.([0, -1 - epsilon, 2 + epsilon])], ["0, -1, 2", "0, -1-$(eps_symbol), 2+$(eps_symbol)"]),
    
    ("Along AD 1", [BigFloat.([0, 0, -1]), BigFloat.([0, 0, -1 + epsilon])], ["0, 0, -1", "0, 0, -1+$(eps_symbol)"]),
    ("Along AD 2", [BigFloat.([0, 0, 2]), BigFloat.([0, 0, 2 + epsilon])], ["0, 0, 2", "0, 0, 2+$(eps_symbol)"]),
    
    ("Outer ABC", [BigFloat.([1/3, 1/3, -0.5]), BigFloat.([1/3, 1/3, -0.5 - epsilon])], ["1/3, 1/3, -0.5", "1/3, 1/3, -0.5-$(eps_symbol)"]),
    ("Outer ABD", [BigFloat.([1/3, -0.5, 1/3]), BigFloat.([1/3, -0.5 - epsilon, 1/3])], ["1/3, -0.5, 1/3", "1/3, -0.5-$(eps_symbol), 1/3"]),
    ("Outer ACD", [BigFloat.([-0.5, 1/3, 1/3]), BigFloat.([-0.5 - epsilon, 1/3, 1/3])], ["-0.5, 1/3, 1/3", "-0.5-$(eps_symbol), 1/3, 1/3"]),
    
    ("Centroid", [BigFloat.([1/4, 1/4, 1/4]), BigFloat.([1/4 + epsilon, 1/4 + epsilon, 1/4 + epsilon])], ["0.25, 0.25, 0.25", "0.25+$(eps_symbol), 0.25+$(eps_symbol), 0.25+$(eps_symbol)"]),
    
    ("Vertex A", [BigFloat.([0, 0, 0]), BigFloat.([0-epsilon, 0-epsilon, 0-epsilon])], ["0, 0, 0", "$(eps_symbol), $(eps_symbol), $(eps_symbol)"]),
    ("Vertex B", [BigFloat.([1, 0, 0]), BigFloat.([1 + epsilon, epsilon, epsilon]), BigFloat.([1 - 3*epsilon, epsilon, epsilon])], ["1, 0, 0", "1+$(eps_symbol), $(eps_symbol), $(eps_symbol)", "1-$('3'eps_symbol), $(eps_symbol), $(eps_symbol)"])
    ]


    # Prepare to collect results for file output
    output_lines = []
    header = @sprintf("%-4s %-15s %-30s %-45s %-45s", "ID", "Location", "Test Point (x, y, z)", "Potential", "Gravitational Acceleration")
    push!(output_lines, header)
    push!(output_lines, "-"^140)
    println(header)
    println("-"^140)


    for (idx, (location, points, labels)) in enumerate(labeled_test_points)
        for (i, (P, label)) in enumerate(zip(points, labels))
            show_id = i == 1 ? string(idx) : ""
            show_loc = i == 1 ? location : ""

            # Calculate potential and gradient
            U = potential_phi(P...)
            gx, gy, gz = compute_gradient(potential_phi, P..., h_grad)

            # Format for printing
            potential_str = @sprintf("%.15f", U)
            gx_str = @sprintf("g_x = %.15f", gx)
            gy_str = @sprintf("g_y = %.15f", gy)
            gz_str = @sprintf("g_z = %.15f", gz)

            # Print to console
            println(@sprintf("%-4s %-15s %-30s %-45s %-45s", show_id, show_loc, label, potential_str, gx_str))
            println(@sprintf("%-4s %-15s %-30s %-45s %-45s", "", "", "", "", gy_str))
            println(@sprintf("%-4s %-15s %-30s %-45s %-45s", "", "", "", "", gz_str))

            # Store for file
            push!(output_lines, @sprintf("%-4s %-15s %-30s %-45s %-45s", show_id, show_loc, label, potential_str, gx_str))
            push!(output_lines, @sprintf("%-4s %-15s %-30s %-45s %-45s", "", "", "", "", gy_str))
            push!(output_lines, @sprintf("%-4s %-15s %-30s %-45s %-45s", "", "", "", "", gz_str))
        end
    end

    # Save the results to a .txt file
    open("gravitational_results_julia.txt", "w") do f
        for line in output_lines
            write(f, line * "\n")
        end
    end
    println("\nResults have been saved to gravitational_results_julia.txt")

end

# Run the main function
main()