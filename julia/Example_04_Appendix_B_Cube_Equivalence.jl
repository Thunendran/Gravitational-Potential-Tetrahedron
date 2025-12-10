# =============================================================================
# Example 04 Appendix B: Cube Equivalence — Superposition of Five Tetrahedra
# -----------------------------------------------------------------------------
# Description:
#     This Julia script demonstrates that the combined gravitational acceleration
#     field from five tetrahedra forming a cube is equivalent to that of a
#     homogeneous cuboid, as shown in Appendix B.
#
#     The implementation uses high-precision BigFloat arithmetic and the
#     analytical tetrahedron potential from the module:
#     `GP_TetrahedronBigFloat_Optimized.jl`.
#
# Reference:
#     Periyandy, T., & Bevis, M. (2025).
#     "The Gravitational Potential Inside, On and Outside of a Homogeneous
#      Tetrahedron." Appendix B
#
# Authors:
#     Thunendran Periyandy (corresponding author)
#     Michael Bevis
#
# Date: November 2025
# =============================================================================

using LinearAlgebra
using StaticArrays
using PyPlot
using PyCall
using Printf

# --- Include the BigFloat tetrahedron module ---
include("GP_TetrahedronBigFloat_Optimized.jl")
using .GP_TetrahedronBigFloat_Optimized

# --- Constants ---
setprecision(256)
const G = BigFloat(1.0)
const σ = BigFloat(1.0)
const h_grad = BigFloat(1e-10)

# -----------------------------------------------------------------------------
# Function: numerical gravitational acceleration via finite difference
# -----------------------------------------------------------------------------
function gravitational_acceleration_numerical(x::BigFloat, y::BigFloat, z::BigFloat,
                                              verts::Matrix{BigFloat})
    A, B, C, D = eachrow(verts)
    tet = GP_TetrahedronBigFloat_Optimized.TetrahedronBigFloat(
        SVector{3,BigFloat}(A),
        SVector{3,BigFloat}(B),
        SVector{3,BigFloat}(C),
        SVector{3,BigFloat}(D)
    )

    φ(P::SVector{3,BigFloat}) =
        GP_TetrahedronBigFloat_Optimized.compute_potential(tet, P; G=G, sigma=σ)

    f(x, y, z) = φ(SVector{3,BigFloat}(x, y, z))

    gx = (f(x + h_grad, y, z) - f(x - h_grad, y, z)) / (2 * h_grad)
    gy = (f(x, y + h_grad, z) - f(x, y - h_grad, z)) / (2 * h_grad)
    gz = (f(x, y, z + h_grad) - f(x, y, z - h_grad)) / (2 * h_grad)

    return gx, gy, gz
end

# -----------------------------------------------------------------------------
# Function: intersection polygon of tetrahedron with z = const plane
# -----------------------------------------------------------------------------
function intersection_polygon(verts::Matrix{BigFloat}; z_plane=BigFloat(0))
    edges = [(1,2), (1,3), (1,4), (2,3), (2,4), (3,4)]
    points = SVector{3,BigFloat}[]

    for (i, j) in edges
        v1 = SVector{3,BigFloat}(verts[i, :])
        v2 = SVector{3,BigFloat}(verts[j, :])
        z1, z2 = v1[3], v2[3]

        if (z1 - z_plane) * (z2 - z_plane) < 0
            t = (z_plane - z1) / (z2 - z1)
            p = v1 .+ t .* (v2 .- v1)
            push!(points, p)
        end
    end

    return points
end

# -----------------------------------------------------------------------------
# Function: compute field (U,V) over grid for one tetrahedron
# -----------------------------------------------------------------------------
function compute_field(vertices::Matrix{BigFloat}, X::Matrix{BigFloat}, Y::Matrix{BigFloat}, z_fixed::BigFloat)
    U = zeros(Float64, size(X))
    V = zeros(Float64, size(Y))
    for i in eachindex(X)
        gx, gy, gz = gravitational_acceleration_numerical(X[i], Y[i], z_fixed, vertices)
        U[i] = Float64(gx)
        V[i] = Float64(gy)
    end
    return U, V
end

# -----------------------------------------------------------------------------
# Define the five tetrahedra forming a cube (BigFloat)
# -----------------------------------------------------------------------------
tetrahedrons = [
    BigFloat.([[-1 -1  1]; [ 1 -1 -1]; [-1  1 -1]; [-1 -1 -1]]),
    BigFloat.([[-1 -1  1]; [ 1  1  1]; [-1  1 -1]; [-1  1  1]]),
    BigFloat.([[ 1 -1 -1]; [ 1  1 -1]; [-1  1 -1]; [ 1  1  1]]),
    BigFloat.([[-1 -1  1]; [ 1  1  1]; [ 1 -1 -1]; [ 1 -1  1]]),
    BigFloat.([[-1 -1  1]; [ 1 -1 -1]; [-1  1 -1]; [ 1  1  1]])
]

# -----------------------------------------------------------------------------
# Define grid in z = 0 plane
# -----------------------------------------------------------------------------
x_vals = BigFloat.(range(-2, 2, length=25))
y_vals = BigFloat.(range(-2, 2, length=25))
X = [x for y in y_vals, x in x_vals]
Y = [y for y in y_vals, x in x_vals]
z_fixed = BigFloat(0)

# -----------------------------------------------------------------------------
# Compute field for each tetrahedron
# -----------------------------------------------------------------------------
println("Computing gravitational acceleration fields for 5 tetrahedra...")
fields = []
polygons = []

for verts in tetrahedrons
    U, V = compute_field(verts, X, Y, z_fixed)
    push!(fields, (U, V))
    push!(polygons, intersection_polygon(verts; z_plane=z_fixed))
end

# -----------------------------------------------------------------------------
# Compute total (superposed) field
# -----------------------------------------------------------------------------
U_total = zeros(Float64, size(X))
V_total = zeros(Float64, size(Y))
for (U, V) in fields
    U_total .+= U
    V_total .+= V
end

# -----------------------------------------------------------------------------
# Plotting setup (robust flattening of PyPlot axes)
# -----------------------------------------------------------------------------

fig, ax_array = subplots(2, 3, figsize=(18, 12))

function collect_all_axes(obj)
    result = PyObject[]
    # Case 1: it's already a single Axes
    if pybuiltin("hasattr")(obj, "quiver") == true
        push!(result, obj)
        return result
    end
    # Case 2: Python iterable (list, tuple, ndarray)
    if pybuiltin("hasattr")(obj, "__iter__") == true
        for sub in obj
            append!(result, collect_all_axes(sub))
        end
    end
    return result
end

axes = collect_all_axes(ax_array)

# Ensure we have exactly 6 axes
if length(axes) < 6
    error("Expected 6 subplots but only got $(length(axes)) — check PyPlot version.")
end


# -----------------------------------------------------------------------------
# Plot each tetrahedron’s field
# -----------------------------------------------------------------------------
titles = [@sprintf("Tetrahedron %d", i) for i in 1:5]
push!(titles, "Combined Field")

for idx in 1:5
    ax = axes[idx]
    U, V = fields[idx]
    ax.quiver(Float64.(X), Float64.(Y), U, V,
              color="red", scale=60, pivot="middle", width=0.003)

    poly = polygons[idx]
    if length(poly) ≥ 3
        xcoords = [Float64(p[1]) for p in poly]
        ycoords = [Float64(p[2]) for p in poly]
        ax.plot(vcat(xcoords, xcoords[1]), vcat(ycoords, ycoords[1]),
                color="black", linewidth=1)
    end

    ax.set_title(titles[idx])
    ax.set_xlim([-2, 2])
    ax.set_ylim([-2, 2])
    ax.set_aspect("equal")
    ax.grid(false)
end

# -----------------------------------------------------------------------------
# Plot combined field (6th panel)
# -----------------------------------------------------------------------------
ax = axes[6]
ax.quiver(Float64.(X), Float64.(Y), U_total, V_total,
          color="red", scale=50, pivot="middle", width=0.003)

for poly in polygons
    if length(poly) ≥ 3
        xcoords = [Float64(p[1]) for p in poly]
        ycoords = [Float64(p[2]) for p in poly]
        ax.plot(vcat(xcoords, xcoords[1]), vcat(ycoords, ycoords[1]),
                color="black", linewidth=1)
    end
end

ax.set_title(titles[6])
ax.set_xlim([-2, 2])
ax.set_ylim([-2, 2])
ax.set_aspect("equal")
ax.grid(false)

# -----------------------------------------------------------------------------
# Final figure formatting
# -----------------------------------------------------------------------------
fig.suptitle("Equivalence of the Gravitational Acceleration of a Cuboid and the Superposition of Five Tetrahedral Accelerations in the Plane z = 0",
             fontsize=20, fontweight="bold")
tight_layout()
savefig("five_tetrahedrons_fields_equivalence_cube_julia.png", dpi=400, bbox_inches="tight")
show()

println("Figure saved as five_tetrahedrons_fields_equivalence_cube_julia.png")
