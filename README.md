# Gravitational Potential of a Homogeneous Tetrahedron  
### Python • MATLAB • Julia Implementations  
### Periyandy & Bevis (2025)

---

## Overview

This repository provides the **first complete, singularity-free analytical algorithm** for computing the **gravitational potential of a homogeneous tetrahedron**:

- Valid **inside**, **on**, and **outside** the tetrahedron  
- Free of all removable singularities (faces, edges, vertices)  
- Fully analytical formulation using stabilized log/atan expressions  
- Orientation-independent implementation  
- Verified using **continuity tests**, **acceleration direction tests**, and **Laplacian tests**

These codes accompany the peer-reviewed publication:

**Periyandy, T. & Bevis, M. (2025).  
_The gravitational potential inside, on, and outside of a homogeneous tetrahedron._  
Journal of Geodesy.**  

DOI: _[insert once assigned]_

---

## Physical Model

### Tetrahedron Geometry  
The gravitational potential is evaluated for a tetrahedron defined by four vertices:

- **A**, **B**, **C**, **D** ∈ ℝ³  
- Vertices may be arbitrarily oriented  
- The algorithm automatically enforces consistent outward face normals  

Given any observation point **P(x, y, z)**, the method computes the gravitational potential using a **closed-form analytical expression** valid for:

- Interior points  
- Boundary points (face, edge, vertex)  
- Exterior points  

### Singularity-Free Analytical Approach  
The formulation removes numerical singularities by:

- Handling **co-planarity** exactly (P on a face)  
- Neutralizing **edge cases** where denominator terms vanish  
- Ensuring **vertex coincidence** does not create undefined atan/log terms  
- Using stable forms of all logarithmic and arctangent components  
- Introducing analytic cancellation of divergent parts  

The resulting gravitational potential is **continuous**, **finite**, and **stable** across the entire spatial domain.

---

## Key Features

- Valid across **all spatial regions**  
- No special casing needed for interior vs exterior points  
- Correct behavior at faces, edges, and vertices  
- Fully analytical expression (no numerical regularization)  
- Orientation-safe implementation using signed volume  
- High-precision support via Julia **BigFloat**  
- Cross-checked consistency across Python, MATLAB, and Julia  
- Includes reproducible **tests**, **examples**, and **benchmark results**

---

## Scientific Validation

### Continuity and Perturbation Tests  
The gravitational potential and acceleration fields were evaluated at:

- All faces  
- All edges  
- All vertices  
- The centroid  
- Exterior points aligned with face planes and edge extensions  

Each point was re-evaluated under **ε-level perturbations** (ε = 10⁻¹⁰).  
Results confirm:

- No numerical instability  
- Smooth transitions across boundaries  
- Correct removal of removable singularities  
- Fully continuous potential field  

---

### Acceleration Directionality Test  
A cube is decomposed into five tetrahedra (four right-angled, one regular).  
The gravitational acceleration is computed for:

- Each tetrahedron individually  
- The superposed sum of all five tetrahedra  

The combined field matches the analytical cuboid gravitational field, confirming:

- Correct vector direction  
- Correct magnitude behavior  
- Accurate representation of interior vs exterior forces  

---

### Laplacian Test  

Poisson’s equation requires:

- **Inside the tetrahedron:** ∇²U = −4πGρ  
- **Outside the tetrahedron:** ∇²U = 0  

A grid-based second-order finite-difference Laplacian was applied to a plane slicing the tetrahedron.  
Results show:

- Constant negative Laplacian inside  
- Zero Laplacian outside  
- Sharp transition along boundary contours  
- High consistency with theoretical Poisson behavior  

This confirms the correctness of the analytical potential and the numerical implementation.


## Repository Structure

```text
Gravitational-Potential-Tetrahedron/
│
├── benchmark/
│   ├── benchmark_summary.txt
│   └── figures/
│
├── data/
│   ├── vertices_example.txt
│   ├── evaluation_points.txt
│   └── additional_test_points.txt
│
├── julia/
│   ├── src/
│   │   ├── GP_TetrahedronBigFloat_Optimized.jl
│   │   ├── G_TetrahedronPotential_Optimized.jl
│   │   └── supporting_modules.jl
│   └── examples/
│
├── matlab/
│   ├── source/
│   └── examples/
│
├── python/
│   ├── source/
│   └── examples/
│
└── README.md   <-- (this file)



