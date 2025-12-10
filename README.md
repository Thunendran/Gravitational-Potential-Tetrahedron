# Gravitational Potential of a Homogeneous Tetrahedron  
### Python • MATLAB • Julia Implementations  
### Periyandy & Bevis (2025)

---

## Overview

This repository provides the first complete, singularity-free analytical algorithm for computing the **gravitational potential of a homogeneous tetrahedron**:

- Valid **inside**, **on**, and **outside** the tetrahedron  
- Removes all removable singularities (faces, edges, vertices)  
- Fully analytical formulation with stable log/atan expressions  
- Orientation-independent implementation  
- Verified using **continuity**, **acceleration**, and **Laplacian** tests

These codes accompany the peer-reviewed publication:

**Periyandy, T. & Bevis, M. (2025).  
_The gravitational potential inside, on, and outside of a homogeneous tetrahedron._  
Journal of Geodesy.**  
DOI: _[]_

---

## Physical Model

### Tetrahedron Geometry  
The potential is computed for a tetrahedron defined by four vertices:

- **A**, **B**, **C**, **D** ∈ ℝ³  
- Arbitrary orientation  
- Outward face normals enforced automatically  

The formulation evaluates the gravitational potential at **any field point** P(x,y,z):

- Interior  
- Faces  
- Edges  
- Vertices  
- Exterior  

### Singularity-Free Formulation  
This method analytically resolves geometrical degeneracies:

- Field point on a face (coplanar)  
- Field point on an edge  
- Field point at a vertex  

All logarithmic and arctangent components are stabilized, ensuring a **continuous potential field** across all regions.

---

## Key Features

- Valid across **all spatial regions**  
- No removable singularities  
- Stable log and arctan evaluations  
- Orientation-safe using signed volume  
- Vectorized Julia, MATLAB, and Python implementations  
- Arbitrary-precision Julia BigFloat version  
- Includes examples, tests, and benchmarking  

---

## Scientific Validation

### Continuity Tests  
Evaluations at faces, edges, and vertices confirm:

- Smooth limiting behavior  
- No jumps or numerical instabilities  
- Proper cancellation of singularities  

### Acceleration Directionality Test  
A cube decomposed into five tetrahedra reproduces the analytical cuboid gravitational field:

- Correct acceleration direction  
- Correct magnitudes  
- Accurate interior/exterior behavior  

### Laplacian Test  
Poisson’s equation requires:

- **Inside:** ∇²U = −4πGρ  
- **Outside:** ∇²U = 0  

Finite-difference tests confirm:

- Constant interior Laplacian  
- Zero exterior Laplacian  
- Correct boundary transition  

---

## Repository Structure

```text
Gravitational-Potential-Tetrahedron/
│
├── benchmark/
│   └── benchmark_summary.txt
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
```
---
## Primary Citation

Periyandy, T. & Bevis, M. (2025).  
**The gravitational potential inside, on, and outside of a homogeneous tetrahedron.**  
*Journal of Geodesy*.  
DOI: _[]_

If you use these codes in research, **please cite this work**.

---

## Contact

For questions or collaborations:
**periyandy.1@osu.edu**  
**thunendran@gmail.com**  

**Division of Geodetic Science**, School of Earth Sciences, The Ohio State University, Columbus, OH 43210, USA  
**Faculty of Geomatics**, Sabaragamuwa University of Sri Lanka, Sri Lanka  

---

## About

Open-source MATLAB, Python, and Julia codes for computing the gravitational potential of a homogeneous tetrahedron.

---

## License

Distributed under the **MIT License**.

