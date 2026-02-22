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

- [x] egui_plot integration for interactive spectra plots (C_ext, C_abs, C_sca vs λ)
- [x] ε(λ) dielectric function plots in the materials panel (Au, Ag, Cu with correct data)
- [x] Near-field |E|² heatmap in the results panel
- [x] 2D dipole lattice scatter preview (xy/xz/yz projection)
- [x] Shape parameter UI for all primitives (sphere, cylinder, cuboid, ellipsoid, helix)
- [x] Helix added to GUI geometry panel (v0.1.1 primitive now exposed in GUI)
- [x] Simulation wired to actual material and geometry panel selections (not hardcoded)

## v0.1.3 — I/O & Export

- [x] Far-field radiation pattern computation (`compute_far_field` in `fields.rs`)
- [x] Far-field polar plot (E-plane / H-plane cuts) in the GUI results panel
- [x] Circular dichroism ΔC_ext computation (`compute_circular_dichroism` in `fields.rs`)
- [x] Incident polarisation selector in GUI (x-pol / y-pol / circular)
- [x] CD checkbox in GUI simulation panel; CD line on spectra plot
- [x] File dialog for CSV export (via `rfd`)
- [x] JSON export option for simulation results (GUI + CLI)
- [x] Metadata header in CSV (version, object list, material, spacing)
- [x] Palik TiO₂ and SiO₂ dielectric data (300–1000 nm) in `lumina-materials`
- [x] TiO₂ / SiO₂ available in GUI materials panel and CLI (`TiO2_Palik`, `SiO2_Palik`)
- [x] Ag_JC and Cu_JC fully wired in CLI material resolver
- [x] CLI geometry expanded: cylinder, cuboid, ellipsoid, helix in TOML configs
- [x] Near-field export to CSV (`save_near_field = true` in TOML output config)
- [x] .xyz file import in GUI geometry panel (Ångström → nm conversion)
- [x] OBJ mesh parser (vertex + face → volume-filling dipole lattice via ray-casting)

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
- [x] Palik dielectric data (TiO₂, SiO₂) — shipped in v0.1.3
- [ ] Extended Palik library (Al₂O₃, Si, GaAs, etc.)
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
