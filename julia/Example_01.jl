# =============================================================================
# Example 01 (BigFloat): Gravitational Potential at a Random Point
# -----------------------------------------------------------------------------
# Description:
#     This example demonstrates the computation of the gravitational potential
#     at a single observation point using the BigFloat version of the Julia
#     implementation of the singularity-free analytical formulation for a
#     homogeneous tetrahedron.
#
#     The vertices (A, B, C, D) define a unit tetrahedron, and the potential
#     is evaluated at a test point P located outside the tetrahedron.
#
# Reference:
#     Periyandy, T., & Bevis, M. (2025).
#     "The Gravitational Potential Inside, On and Outside of a Homogeneous Tetrahedron." 
#     Supplementary Material:
#     Julia Implementation – High-Precision (BigFloat) Version
#
# Authors:
#     Thunendran Periyandy  (corresponding author)
#     Michael Bevis
#
# Date: November 2025
# =============================================================================

using StaticArrays
using Printf

# --- Include the BigFloat module file ---
include("GP_TetrahedronBigFloat_Optimized.jl")
using .GP_TetrahedronBigFloat_Optimized

# --- Set BigFloat precision: 50 significant digits (base 10) ---
setprecision(BigFloat, 50)

# --- Define vertices and evaluation point (as BigFloat vectors) ---
vA = @SVector BigFloat[1.0, 1.0, 1.0]
vB = @SVector BigFloat[2.0, 1.0, 1.0]
vC = @SVector BigFloat[1.0, 2.0, 1.0]
vD = @SVector BigFloat[1.0, 1.0, 2.0]
p_eval = @SVector BigFloat[0.0, 0.0, 0.0]

# --- Create tetrahedron object (with module prefix) ---
tet = GP_TetrahedronBigFloat_Optimized.TetrahedronBigFloat(vA, vB, vC, vD)

# --- Compute gravitational potential ---
potential = GP_TetrahedronBigFloat_Optimized.compute_potential(
    tet, p_eval; G=BigFloat(1), sigma=BigFloat(1)
)

# --- Print result with 50 significant digits ---
@printf("Gravitational Potential (50 digits):\n%.50e\n", potential)
