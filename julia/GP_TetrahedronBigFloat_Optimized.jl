# =============================================================================
#  Module: G_TetrahedronBigFloat_Optimized
# -----------------------------------------------------------------------------
#  Description:
#      High-precision (BigFloat) Julia implementation of the analytical
#      formulation for computing the gravitational potential of a homogeneous
#      tetrahedron. This version is designed for numerical validation,
#      precision benchmarking, and continuity testing across tetrahedral
#      boundaries.
#
#      The module defines:
#          • TetrahedronBigFloat    –  precomputes all geometric relations
#                                      using BigFloat arithmetic.
#          • compute_potential      –  computes the potential at a given point.
#          • compute_potentials     –  evaluates potentials at multiple points
#                                      using multithreading.
#
#      The BigFloat environment enables evaluation at arbitrary precision
#      (set via `setprecision()`), allowing exact replication of the analytical
#      tests described in Appendix A (continuity of potential and acceleration).
#
#  Numerical Scope:
#      Uses arbitrary precision via `BigFloat`. Recommended for validation
#      against analytical limits, precision studies, and convergence tests.
#
#  Reference:
#      Periyandy, T., & Bevis, M. (2025).
#      "The Gravitational Potential Inside, On and Outside of a Homogeneous Tetrahedron."
#
#      Supplementary Material:
#      Julia Implementation — High-Precision (BigFloat) Optimized Version
#
#  Authors:
#      Thunendran Periyandy  (corresponding author)
#          Division of Geodetic Science, School of Earth Sciences,
#          The Ohio State University, USA
#          Faculty of Geomatics, Sabaragamuwa University of Sri Lanka
#          Email: thunendran@gmail.com
#
#      Michael Bevis
#          Division of Geodetic Science, School of Earth Sciences,
#          The Ohio State University, USA
#          Email: mbevis@gmail.com
#
#  Version: 1.0
#  Date:    November 2025
# =============================================================================





module GP_TetrahedronBigFloat_Optimized

using LinearAlgebra
using StaticArrays
using Base.Threads

# ============================================================
# Struct Definition
# ============================================================
"""
    struct TetrahedronBigFloat

Stores precomputed geometric data for a homogeneous tetrahedron
with vertices (A, B, C, D) using BigFloat arithmetic.
"""
struct TetrahedronBigFloat
    verts::NTuple{4, SVector{3,BigFloat}}
    centroid::SVector{3,BigFloat}
    edges::Matrix{SVector{3,BigFloat}}   # 4×4 edge vectors
    lengths::Matrix{BigFloat}            # 4×4 edge lengths
end

# ============================================================
# Helper Functions
# ============================================================

# Preserve tetrahedron orientation (ensures positive volume)
@inline function preserve_orientation_tetrahedron(A, B, C, D)
    V = dot(cross(B - A, C - A), D - A)
    return V < 0 ? (A, C, B, D) : (A, B, C, D)
end

# Special logarithmic term
@inline function special_log_term(ra::BigFloat, rb::BigFloat, rab::BigFloat)
    if rab <= 0
        return BigFloat(0)
    end
    num = ra + rb + rab
    den = ra + rb - rab
    return den <= 0 ? BigFloat(0) : log(num / den) / rab
end

# Special arctangent term
@inline function special_atan_term(rval::BigFloat, detA::BigFloat, detB::BigFloat,
                                dotv::BigFloat, diffz::BigFloat, num::BigFloat)
    return rval > 0 ? atan(num, -((detA * detB) + diffz^2 * dotv) / rval) : BigFloat(0)
end

# ============================================================
# Constructor
# ============================================================
"""
    TetrahedronBigFloat(A, B, C, D)

Creates a tetrahedron object with precomputed edge geometry.
"""
function TetrahedronBigFloat(A::SVector{3,BigFloat},
                             B::SVector{3,BigFloat},
                             C::SVector{3,BigFloat},
                             D::SVector{3,BigFloat})
    A, B, C, D = preserve_orientation_tetrahedron(A, B, C, D)
    verts = (A, B, C, D)
    centroid = (A + B + C + D) / BigFloat(4)

    edges = Matrix{SVector{3,BigFloat}}(undef, 4, 4)
    lengths = Matrix{BigFloat}(undef, 4, 4)
    for i in 1:4, j in 1:4
        if i != j
            v = verts[j] - verts[i]
            edges[i,j] = v
            lengths[i,j] = norm(v)
        else
            edges[i,j] = SVector{3,BigFloat}(0,0,0)
            lengths[i,j] = BigFloat(0)
        end
    end
    return TetrahedronBigFloat(verts, centroid, edges, lengths)
end

# ============================================================
# Compute Potential at One Point
# ============================================================
function compute_potential(tet::TetrahedronBigFloat, P::SVector{3,BigFloat};
                           G=BigFloat(1), sigma=BigFloat(1))
    (A, B, C, D) = tet.verts
    faces = ((1,3,2), (1,2,4), (2,3,4), (3,1,4))
    total_U = BigFloat(0)

    for (i1, i2, i3) in faces
        v1, v2, v3 = tet.verts[i1], tet.verts[i2], tet.verts[i3]
        n_raw = cross(v2 - v1, v3 - v1)
        if dot(n_raw, (v1 + v2 + v3)/BigFloat(3) - tet.centroid) < 0
            n_raw = -n_raw
        end
        norm_val = norm(n_raw)
        n_unit = norm_val > 0 ? n_raw / norm_val : SVector{3,BigFloat}(0,0,0)

        diff_z = dot(v1 - P, n_unit)
        numerator = diff_z * norm_val

        S = [cross(tet.verts[i] - P, tet.verts[j] - P) for i in 1:4, j in 1:4]
        det12 = dot(S[i1,i2], n_unit)
        det23 = dot(S[i2,i3], n_unit)
        det31 = dot(S[i3,i1], n_unit)

        r = (norm(v1 - P), norm(v2 - P), norm(v3 - P))
        edges = (tet.lengths[i1,i2], tet.lengths[i2,i3], tet.lengths[i3,i1])

        L12 = special_log_term(r[1], r[2], edges[1])
        L23 = special_log_term(r[2], r[3], edges[2])
        L31 = special_log_term(r[3], r[1], edges[3])

        dot1 = dot(tet.edges[i1,i2], tet.edges[i3,i1])
        dot2 = dot(tet.edges[i2,i3], tet.edges[i1,i2])
        dot3 = dot(tet.edges[i3,i1], tet.edges[i2,i3])

        S1 = special_atan_term(r[1], det31, det12, dot1, diff_z, numerator)
        S2 = special_atan_term(r[2], det12, det23, dot2, diff_z, numerator)
        S3 = special_atan_term(r[3], det23, det31, dot3, diff_z, numerator)

        term1 = diff_z * (det12 * L12 + det23 * L23 + det31 * L31)
        term2 = (diff_z^2) * (S1 + S2 + S3 - sign(diff_z) * BigFloat(π))
        total_U += 0.5 * (term1 - term2)
    end
    return G * sigma * total_U
end

# ============================================================
# Compute Potentials for Multiple Points (Threaded)
# ============================================================
function compute_potentials(tet::TetrahedronBigFloat, points::Matrix{BigFloat};
                            G=BigFloat(1), sigma=BigFloat(1))
    N = size(points, 1)
    U = Vector{BigFloat}(undef, N)
    @threads for i in 1:N
        P = SVector{3,BigFloat}(points[i,1], points[i,2], points[i,3])
        U[i] = compute_potential(tet, P; G=G, sigma=sigma)
    end
    return U
end

end # module
