# Version Roadmap

## v0.1 — Linear Solver Foundation

- [x] Workspace structure with six crates
- [x] Dyadic Green's function implementation
- [x] Interaction matrix assembly with Rayon parallelism
- [x] Direct solver (LU via `faer`) for N <= 1000
- [ ] GMRES iterative solver for N > 1000
- [x] Clausius-Mossotti polarisability with radiative corrections (RRCM + LDR)
- [x] Extinction, absorption, scattering cross-section computation
- [ ] Near-field |E|^2 maps on observation planes
- [x] Johnson & Christy Au data (43 points, 188–892 nm) with cubic spline interpolation
- [x] Geometry primitives (sphere, cylinder, cuboid, ellipsoid, helix)
- [x] Discretisation: centred cubic lattice from primitives
- [ ] .xyz file parser
- [x] GUI dashboard with geometry, materials, simulation, and results panels
- [x] CLI with TOML configuration
- [x] Mie theory validation (10 nm Au sphere, dielectric spheres)
- [ ] GPU matrix assembly via wgpu (CPU fallback)
- [x] mdBook documentation with validation results

### v0.1 Validation Summary

The CDA solver has been validated against Mie theory with the following results:

- **Dielectric spheres**: < 15% error at d = 3 nm, < 5% at d ≤ 2 nm
- **Gold (interband, 420–510 nm)**: < 30% error at d = 3 nm
- **Gold (Drude, > 550 nm)**: > 100% error — known limitation of CDA at coarse discretisation
- **Single-dipole Rayleigh limit**: exact agreement (0.00% error)
- **Convergence**: monotonic for dielectrics, non-monotonic for metals

See [v0.1 Validation](../theory/validation.md) for full details.

## v0.2 — Nonlinear & Time Domain

- [ ] SHG source term computation from local fields
- [ ] THG source term computation
- [ ] Symmetry group handling for chi(2) tensors
- [ ] Time-domain field reconstruction via inverse FFT
- [ ] Far-field radiation pattern computation
- [ ] Result export: CSV, JSON, HDF5
- [ ] egui_plot integration for in-GUI spectra plotting
- [ ] 3D viewport for dipole lattice visualisation
- [ ] Filtered Coupled Dipole (FCD) for metallic convergence
- [ ] Surface averaging of polarisabilities
- [ ] GMRES iterative solver for large N

## v0.3 — Periodic Structures & HPC

- [ ] Ewald summation for 1D/2D/3D periodicity
- [ ] Bloch boundary conditions for oblique incidence
- [ ] MPI distributed computing backend
- [ ] Hybrid MPI + Rayon parallelism
- [ ] SLURM job submission helpers
- [ ] Palik dielectric data library
- [ ] DFT polarisability import (VASP, Gaussian)

## v1.0 — Multi-Method & Cross-Platform

- [ ] Boundary Element Method (BEM) solver
- [ ] T-matrix solver
- [ ] FDTD time-domain solver
- [ ] CUDA-specific compute backend (cuSOLVER/cuBLAS)
- [ ] Metal Performance Shaders backend
- [ ] Cross-platform binary releases
- [ ] Comprehensive test suite with property-based testing
- [ ] Full documentation with tutorials and examples
