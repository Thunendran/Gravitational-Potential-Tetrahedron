# =============================================================================
#  Module: G_TetrahedronPotential
# -----------------------------------------------------------------------------
#  Description:
#      Julia implementation of the singularity-free analytical formulation
#      for computing the gravitational potential of a homogeneous tetrahedron.
#
#      This module defines the `TetrahedronPotential` type and provides
#      two main functions:
#          - `compute_potential`  : Computes potential at a single point.
#          - `compute_potentials` : Computes potentials at multiple points using multithreading for performance.
#
#      The formulation ensures numerical stability and analytical continuity
#      across the interior, boundary, and exterior of the tetrahedron by
#      resolving singularities in logarithmic and arctangent terms.
#
#  Reference:
#      Periyandy, T., & Bevis, M. (2025).
#      "The Gravitational Potential Inside, On and Outside of a Homogeneous Tetrahedron."
#
#      Supplementary Material:
#      Julia Implementation (High-Precision and Threaded Version)
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

module G_TetrahedronPotential

using LinearAlgebra, StaticArrays, SpecialFunctions, Base.Threads

export TetrahedronPotential, compute_potential, compute_potentials

# ==========================================================
# Type Definition
# ==========================================================
struct TetrahedronPotential{T}
    A::SVector{3,T}
    B::SVector{3,T}
    C::SVector{3,T}
    D::SVector{3,T}
    G::T
    σ::T
end

# ----------------------------------------------------------
# Constructor (orientation preserved)
# ----------------------------------------------------------
function TetrahedronPotential(A::AbstractVector{T}, B::AbstractVector{T},
                              C::AbstractVector{T}, D::AbstractVector{T};
                              G::T = one(T), σ::T = one(T)) where {T<:Real}
    V = dot(cross(B - A, C - A), D - A)
    if V < 0
        B, C = C, B
    end
    return TetrahedronPotential(SVector{3,T}(A), SVector{3,T}(B),
                                SVector{3,T}(C), SVector{3,T}(D), G, σ)
end

# ----------------------------------------------------------
# Inline helpers
# ----------------------------------------------------------
@inline safe_log_term(ra, rb, rab) =
    (rab > 0 && ra + rb > rab) ? log((ra + rb + rab) / (ra + rb - rab)) / rab : zero(ra)

@inline safe_atan_term(rval, detA, detB, dotv, diffz, num) =
    rval > 0 ? atan(num, -((detA * detB) + diffz^2 * dotv) / rval) : zero(rval)

# ----------------------------------------------------------
# Potential at single point
# ----------------------------------------------------------
function compute_potential(tet::TetrahedronPotential{T}, P::SVector{3,T}) where {T}
    A, B, C, D, G, σ = tet.A, tet.B, tet.C, tet.D, tet.G, tet.σ
    centroid = (A + B + C + D) / 4
    verts = (A, B, C, D)
    faces = ((1, 3, 2), (1, 2, 4), (2, 3, 4), (3, 1, 4))

    # --- Vertex distances ---
    vertex_r = SVector(norm(verts[1] - P), norm(verts[2] - P),
                       norm(verts[3] - P), norm(verts[4] - P))

    # --- Edge lengths and vectors ---
    edge_len = Matrix{T}(undef, 4, 4)
    edge_vec = Matrix{SVector{3,T}}(undef, 4, 4)
    for i in 1:4, j in 1:4
        edge_vec[i, j] = verts[j] - verts[i]
        edge_len[i, j] = norm(edge_vec[i, j])
    end

    # --- Cross products S_ij ---
    S = Matrix{SVector{3,T}}(undef, 4, 4)
    for i in 1:4, j in 1:4
        S[i, j] = cross(verts[i] - P, verts[j] - P)
    end

    totalU = zero(T)
    @inbounds for (n1, n2, n3) in faces
        v1, v2, v3 = verts[n1], verts[n2], verts[n3]
        nraw = cross(v2 - v1, v3 - v1)
        if dot(nraw, (v1 + v2 + v3) / 3 - centroid) < 0
            nraw = -nraw
        end
        normn = norm(nraw)
        nunit = normn > 0 ? nraw / normn : zeros(T, 3)
        diffz = dot(v1 - P, nunit)
        num = diffz * normn

        det12 = dot(S[n1, n2], nunit)
        det23 = dot(S[n2, n3], nunit)
        det31 = dot(S[n3, n1], nunit)
        r1, r2, r3 = vertex_r[n1], vertex_r[n2], vertex_r[n3]
        e12, e23, e31 = edge_len[n1, n2], edge_len[n2, n3], edge_len[n3, n1]

        # --- Log terms ---
        L12 = safe_log_term(r1, r2, e12)
        L23 = safe_log_term(r2, r3, e23)
        L31 = safe_log_term(r3, r1, e31)
        term1 = diffz * (det12 * L12 + det23 * L23 + det31 * L31)

        # --- Arctangent terms ---
        dot1 = dot(edge_vec[n1, n2], edge_vec[n3, n1])
        dot2 = dot(edge_vec[n2, n3], edge_vec[n1, n2])
        dot3 = dot(edge_vec[n3, n1], edge_vec[n2, n3])
        S1 = safe_atan_term(r1, det31, det12, dot1, diffz, num)
        S2 = safe_atan_term(r2, det12, det23, dot2, diffz, num)
        S3 = safe_atan_term(r3, det23, det31, dot3, diffz, num)
        term2 = (diffz^2) * (S1 + S2 + S3 - sign(diffz) * π)

        totalU += 0.5 * (term1 - term2)
    end

    return G * σ * totalU
end

# ----------------------------------------------------------
# Vectorized potential (multithreaded)
# ----------------------------------------------------------
function compute_potentials(tet::TetrahedronPotential{T}, points::Matrix{T}) where {T}
    N = size(points, 1)
    U = Vector{T}(undef, N)
    @threads for i in 1:N
        P = SVector{3,T}(points[i, 1], points[i, 2], points[i, 3])
        U[i] = compute_potential(tet, P)
    end
    return U
end

end # module
