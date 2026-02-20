# Version Roadmap

## v0.1.0 — Linear Solver Foundation

- [x] Workspace structure with six crates
- [x] Dyadic Green's function implementation
- [x] Interaction matrix assembly with Rayon parallelism
- [x] Direct solver (LU via `faer`) for N <= 1000
- [x] Clausius-Mossotti polarisability with radiative corrections (RRCM + LDR)
- [x] Extinction, absorption, scattering cross-section computation
- [x] Johnson & Christy Au data (43 points, 188–892 nm) with cubic spline interpolation
- [x] Geometry primitives (sphere, cuboid, ellipsoid)
- [x] Discretisation: centred cubic lattice from primitives
- [x] GUI dashboard with geometry, materials, simulation, and results panels
- [x] CLI with TOML configuration
- [x] Mie theory validation (10 nm Au sphere, dielectric spheres)
- [x] mdBook documentation with validation results

## v0.1.1 — Accurate & Scalable Linear Solver

- [x] GMRES(m) iterative solver with restart for N > 1000
- [x] Filtered Coupled Dipole (FCD) Green's function (IGT method)
- [x] FCD wired as default in CdaSolver
- [x] Cylinder containment and bounding box
- [x] Helix containment and bounding box
- [x] GMRES vs direct LU validation (machine-precision agreement)
- [x] FCD gold sphere full-spectrum validation (420–800 nm)
- [x] Cylinder and helix discretisation tests

### v0.1.1 Validation Summary

- **GMRES**: Agrees with direct LU to ~10⁻¹³ relative error
- **FCD interband (420–510 nm)**: 13–24% error (improved from 22–28% with point-dipole)
- **FCD Drude (>550 nm)**: Still >80% error — known limitation at d=2nm
- **Cylinder/Helix**: Correct discretisation with expected dipole counts

## v0.1.2 — Plotting & Visualisation

- [ ] egui_plot integration for interactive spectra plots (C_ext, C_abs, C_sca vs λ)
- [ ] ε(λ) dielectric function plots in the materials panel
- [ ] Near-field |E|² heatmap in the results panel
- [ ] 2D dipole lattice scatter preview (xy/xz/yz projection)
- [ ] Shape parameter UI for all primitives (cylinder, cuboid, ellipsoid)

## v0.1.3 — I/O & Export

- [ ] Far-field radiation pattern computation and polar plots
- [ ] File dialog for CSV export path
- [ ] JSON export option for simulation results
- [ ] Metadata header in CSV (simulation params, date, version)
- [ ] OBJ mesh parser (vertex + face → volume-filling dipole lattice)
- [ ] .xyz file parser

## v0.2.0 — GPU Compute Engine

- [ ] wgpu compute shaders for matrix assembly
- [ ] FFT-accelerated matvec for regular lattices (block-Toeplitz)
- [ ] GPU GMRES with CPU Arnoldi + GPU matvec
- [ ] CPU fallback parity via ComputeBackend trait
- [ ] Performance benchmarks (CPU vs GPU, N = 1k–50k)
- [ ] Near-field |E|² maps on observation planes

## v0.3.0 — Nonlinear Optics & Periodicity

- [ ] SHG source term computation from local fields
- [ ] THG source term computation
- [ ] χ(2) symmetry group handling
- [ ] Time-domain field reconstruction via inverse FFT
- [ ] Ewald summation for 1D/2D/3D periodicity
- [ ] Bloch boundary conditions for oblique incidence
- [ ] MPI distributed computing backend
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
