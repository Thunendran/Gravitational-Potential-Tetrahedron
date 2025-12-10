# =============================================================================
# Example 05: Benchmarking Float64 vs BigFloat Implementations
# -----------------------------------------------------------------------------
# Description:
#     This script benchmarks the computational performance and numerical
#     precision of the analytical gravitational potential calculation for a
#     homogeneous tetrahedron using two implementations:
#
#         (1) G_TetrahedronPotential.jl          → Machine-precision (Float64)
#         (2) GP_TetrahedronBigFloat_Optimized.jl → High-precision (BigFloat)
#
#     It computes potentials for a large set of evaluation points (e.g., 100k)
#     and reports both total runtime and per-point average time for each case.
#
#     The comparison highlights:
#         • The computational cost of high-precision arithmetic.
#         • Consistency of potential results across precisions.
#
# Reference:
#     Periyandy, T. & Bevis, M. (2025)
#     “The Gravitational Potential Inside, On and Outside of a Homogeneous
#      Tetrahedron.” Supplementary Material, Appendix A & B.
#
# Authors:
#     Thunendran Periyandy (corresponding author)
#     Michael Bevis
#
# Version: 1.0
# Date:    November 2025
# =============================================================================

using DelimitedFiles
using LinearAlgebra
using StaticArrays
using Printf

# -----------------------------------------------------------------------------
# Include the computational modules
# -----------------------------------------------------------------------------
include("G_TetrahedronPotential.jl")             # Float64 optimized version
include("GP_TetrahedronBigFloat_Optimized.jl")   # BigFloat high-precision version

# -----------------------------------------------------------------------------
# Load tetrahedron geometry and evaluation points
# -----------------------------------------------------------------------------
# File paths:
#   ../data/vertices.txt  → contains 4 rows (A, B, C, D)
#   ../data/points_100k.txt → contains N rows of (x, y, z)
# -----------------------------------------------------------------------------
A = readdlm("../data/vertices.txt")
points = readdlm("../data/points_100k.txt")

# Extract vertex coordinates
A1, B, C, D = A[1,:], A[2,:], A[3,:], A[4,:]

println("───────────────────────────────────────────────")
println("Running Julia Float64 benchmark (Optimized)...")

# =============================================================================
#  FLOAT64 BENCHMARK
# =============================================================================
# Description:
#   Evaluates potential using double precision (machine precision) arithmetic.
#   This version is optimized for speed and used as the baseline reference.
# =============================================================================

# Construct the Float64 tetrahedron object
tetra = G_TetrahedronPotential.TetrahedronPotential(A1, B, C, D)

# Warm-up pass (to trigger compilation before timing)
G_TetrahedronPotential.compute_potentials(tetra, points[1:2, :])

# Benchmark execution
t_start = time()
U_f64 = G_TetrahedronPotential.compute_potentials(tetra, points)
t_elapsed = time() - t_start

# Save results and timing
writedlm("julia_float_results.txt", U_f64)
open("julia_float_time.txt", "w") do f
    @printf(f, "%.6f", t_elapsed)
end

println("Total runtime (Float64): ", round(t_elapsed, digits=6), " s")
@printf("Average time per point: %.3e s/pt\n", t_elapsed / size(points, 1))

# =============================================================================
#  BIGFLOAT BENCHMARK
# =============================================================================
# Description:
#   Evaluates potential using arbitrary-precision BigFloat arithmetic.
#   Provides high-precision validation for analytical correctness and
#   continuity testing across boundaries.
# =============================================================================
println("\nRunning Julia BigFloat benchmark (Optimized)...")

# Convert tetrahedron vertices to SVector{3,BigFloat}
A_b = SVector{3,BigFloat}(A[1,1], A[1,2], A[1,3])
B_b = SVector{3,BigFloat}(A[2,1], A[2,2], A[2,3])
C_b = SVector{3,BigFloat}(A[3,1], A[3,2], A[3,3])
D_b = SVector{3,BigFloat}(A[4,1], A[4,2], A[4,3])

# Convert evaluation points to BigFloat matrix
points_b = Matrix{BigFloat}(points)

# Create BigFloat tetrahedron object
tetra_big = GP_TetrahedronBigFloat_Optimized.TetrahedronBigFloat(A_b, B_b, C_b, D_b)

# Warm-up compilation
GP_TetrahedronBigFloat_Optimized.compute_potentials(tetra_big, points_b[1:2, :])

# Benchmark execution
t_start = time()
U_big = GP_TetrahedronBigFloat_Optimized.compute_potentials(tetra_big, points_b)
t_elapsed = time() - t_start

# Save results and timing
writedlm("julia_bigfloat_results.txt", U_big)
open("julia_bigfloat_time.txt", "w") do f
    @printf(f, "%.6f", t_elapsed)
end

println("Total runtime (BigFloat): ", round(t_elapsed, digits=6), " s")
@printf("Average time per point: %.3e s/pt\n", t_elapsed / size(points_b, 1))
println("───────────────────────────────────────────────")

# =============================================================================
# End of Benchmark Script
# =============================================================================
