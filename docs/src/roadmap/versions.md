# Version Roadmap

## v0.1 — Linear Solver Foundation

- [x] Workspace structure with six crates
- [ ] Dyadic Green's function implementation
- [ ] Interaction matrix assembly with Rayon parallelism
- [ ] Direct solver (LU via `faer`) for N <= 1000
- [ ] GMRES iterative solver for N > 1000
- [ ] Clausius-Mossotti polarisability with radiative corrections
- [ ] Extinction, absorption, scattering cross-section computation
- [ ] Near-field |E|^2 maps on observation planes
- [ ] Johnson & Christy Au/Ag/Cu with cubic spline interpolation
- [ ] Geometry primitives (sphere, cylinder, cuboid, ellipsoid, helix)
- [ ] Discretisation: primitive to dipole lattice
- [ ] .xyz file parser
- [ ] GUI dashboard with geometry, materials, simulation, and results panels
- [ ] CLI with TOML configuration
- [ ] Mie theory validation (20 nm Au sphere)
- [ ] GPU matrix assembly via wgpu (CPU fallback)
- [ ] mdBook documentation

## v0.2 — Nonlinear & Time Domain

- [ ] SHG source term computation from local fields
- [ ] THG source term computation
- [ ] Symmetry group handling for chi(2) tensors
- [ ] Time-domain field reconstruction via inverse FFT
- [ ] Far-field radiation pattern computation
- [ ] Result export: CSV, JSON, HDF5
- [ ] egui_plot integration for in-GUI spectra plotting
- [ ] 3D viewport for dipole lattice visualisation

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
