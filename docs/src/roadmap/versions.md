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

- [x] Matrix-free GMRES interface (matvec closure instead of matrix reference)
- [x] ComputeBackend trait wired into CdaSolver (CPU default, GPU optional)
- [x] Rayon-parallel interaction matrix assembly (off-diagonal Green's tensor blocks)
- [x] wgpu WGSL compute shader for complex matrix-vector product (f32)
- [x] GpuBackend with matrix caching (upload once, reuse for ~300 GMRES iterations)
- [x] GPU feature forwarding through crate chain (lumina-compute → core → gui/cli)
- [x] GPU toggle in GUI simulation panel
- [x] Backend selection in CLI TOML config (`backend = "auto" | "cpu" | "gpu"`)
- [x] Automatic CPU fallback when GPU unavailable
- [x] GPU vs CPU benchmark harness (matvec + GMRES, N = 300–3000)
- [x] All 85 tests passing (81 existing + 4 GPU)

### v0.2.0 Performance Notes

- GPU matvec is correct (agrees with CPU to ~1e-6, the f32 precision limit)
- At N ≤ 3000, GPU is slower than CPU due to per-call buffer overhead
- Speedup expected at N > 5000 with persistent buffer reuse (v0.2.1)

## v0.2.1 — Performance & XYZ Workflow

### Performance

- [x] Parallel wavelength sweep in GUI (rayon `par_iter` — matches CLI performance)
- [x] Parallel wavelength sweep in CLI (rayon `par_iter`, `AtomicUsize` progress)
- [x] Stack-allocated Green's tensors (`[[Complex64; 3]; 3]` — eliminates ~N² heap allocations per wavelength)
- [x] Direct-write matrix assembly (raw pointer writes to disjoint row groups — removes intermediate `Vec` allocation)

### XYZ Import Workflow

- [x] Auto-detected nearest-neighbour distance for XYZ imports (correct FCD cell size and CM volume)
- [x] Directory-scoped file browser in GUI geometry panel (pick folder, then select from dropdown)
- [x] Example Au Mackay icosahedron (k=6, 923 atoms) supplied in `examples/au_icosahedron_k6.xyz`
- [x] Mackay icosahedron test suite in `lumina-geometry` (atom counts, symmetry, distances)

### GUI Enhancements

- [x] Debug output toggle in simulation panel (per-wavelength timing, matrix size, ε values, solver method)
- [x] Scrollable monospace debug log area with auto-scroll
- [x] Folder browser + file dropdown for `.xyz`/`.obj` import (with manual refresh)

### v0.2.1 Performance Notes

- GUI wavelength sweep now runs in parallel on all CPU cores (previously sequential)
- For a 923-dipole icosahedron with 51 wavelengths, runtime drops from ~90 minutes to ~30 seconds (release mode)
- Stack-allocated Green's tensors eliminate ~850K heap allocations per wavelength for N=923
- Direct-write assembly removes ~136 MB intermediate buffer for N=923

### Deferred to v0.2.2+

- [x] Persistent GPU buffers (eliminate per-matvec allocation overhead)
- [x] Complex assembly builder (SceneSpec, multi-object scenes)
- [ ] FFT-accelerated matvec for regular lattices (block-Toeplitz)
- [ ] GPU matrix assembly shader
- [ ] Surface averaging for metallic convergence

## v0.2.2 — Persistent Buffers & Scene Assembly

### Persistent GPU Session Buffers

- [x] `MatvecSession` trait: `upload_matrix` + `matvec` — matrix uploaded once per wavelength
- [x] `GpuSession`: 5 pre-allocated WGSL buffers + pre-built bind group, reused across all GMRES iterations
- [x] `GpuBackend` holds `Arc<Mutex<GpuState>>` — thread-safe across Rayon wavelength workers
- [x] Per-matvec overhead eliminated: 7 allocations → 0, 300 matrix writes → 1 per wavelength

### SceneSpec — Complex Assembly Builder

- [x] `SceneSpec` in `lumina-geometry`: shared GUI+CLI scene description, serialisable to/from TOML
- [x] `[[geometry.object]]` arrays: each entry is an `ObjectSpec` with name, shape, material, spacing, transform
- [x] `ObjectTransform`: position (nm), rotation\_deg (Euler XYZ), scale factor
- [x] CoreShell geometry: concentric shell layers, each with its own material
- [x] `SceneSpec::build_geometry()` → flat `SceneDipoles { positions, materials, spacings, hints }`
- [x] GUI: `ScenePanel` replaces old `GeometryPanel` + `MaterialsPanel`
- [x] CLI: `JobConfig.geometry: SceneSpec`; per-dipole material resolution in runner

### v0.2.2 Performance Notes

- GPU session buffers eliminate all per-iteration allocation overhead on the GPU path
- At N = 1000–3000, GPU path now competitive with CPU; speedup expected at N > 5000

## v0.3.0 — Periodic Structures & Extended Geometry

### Periodic Structures

- [x] `LatticeSpec` in `lumina-geometry`: 1D chain or 2D planar lattice specification
- [x] Direct real-space lattice sum in `EwaldGreens` (finite shell cutoff, exact for small unit cells)
- [x] `SimulationParams.k_bloch: [f64; 3]` — Bloch wavevector for oblique incidence
- [x] New result types: `BlochCrossSections`, `DispersionMap` for band structure sweeps

### Extended Geometry

- [x] CoreShell (concentric ellipsoidal shells, arbitrary number of layers)
- [x] Anisotropic ellipsoid dipoles: per-dipole 3×3 polarisability tensor from depolarisation factors
- [x] Depolarisation factors via 16-point Gauss–Legendre quadrature (`ellipsoid_depol_factors`)
- [x] `SubstrateSpec`: image-dipole substrate (quasi-static Fresnel reflection, flat interface below the structure)

### Other

- [x] Palik dielectric data (TiO₂, SiO₂) — shipped in v0.1.3
- [ ] Extended Palik library (Al₂O₃, Si, GaAs, etc.)

## v0.4.0 — Nonlinear Optics & Accelerated Periodic Solver

### Ewald Acceleration

- [x] `EwaldGreens` upgraded from direct sum to full Ewald split
- [x] Real-space sum: erfc-damped (`erfc_real(η·|R|)` multiplied per image) — exponential convergence, ≤ 5 shells sufficient
- [x] Reciprocal-space sum (2D lattices): spectral G-sum `(i/2A)·Σ_G exp(iQ·ρ)·exp(-q|z|)/q · M(Q,q)`
- [x] `csqrt_positive_imag`: Im≥0 branch selection for evanescent vs. propagating wave convention
- [x] Combined sum: `evaluate()` = `real_space() + recip_space()` (API unchanged)

### Retarded Substrate

- [x] `fresnel_reflection_amplitude(eps_sub, eps_env)` — retarded normal-incidence reflection coefficient r = (n_sub − n_env)/(n_sub + n_env)
- [x] `substrate_reflection_factor(..., use_retarded)` dispatcher
- [x] `SubstrateSpec.use_retarded: bool` (default `true`, backward-compatible TOML)
- [x] GUI and CLI updated to use retarded coefficient by default

### Second-Harmonic Generation (SHG)

- [x] `Chi2Tensor` in `lumina-core::types`: rank-3 [Complex64; 27], `zero()`, `isotropic_surface(chi_zzz, chi_zxx)`, `contract()`
- [x] `compute_shg_sources(local_fields, chi2_tensors)` — χ²:EE source term per dipole
- [x] `compute_shg_response(solver, omega_response, dipoles_2omega, chi2_tensors, params_2omega, calc_far_field)` — full self-consistent SHG solve
- [x] `ShgResult`: source_moments, driven_moments, shg_intensity (nm⁶, ∝|E|⁴), optional far-field pattern
- [x] `CdaSolver::cross_sections_from_response` — avoids double-solve when SHG is also computed
- [x] CLI: `[nonlinear]` TOML section; `shg_spectrum.csv` output
- [x] GUI: SHG panel with symmetry selector, χ DragValue inputs, SHG spectra tab in results
- [x] 5 unit tests covering zero χ², χ_xxx contraction, C∞v tensor structure, end-to-end solve, E⁴ scaling

### Third-Harmonic Generation (THG)

- [x] `Chi3Tensor` in `lumina-core::types`: rank-4 [Complex64; 81], index a*27+b*9+c*3+d, `zero()`, `isotropic_bulk(chi_xxxx, chi_xxyy)`, `contract()`
- [x] `compute_thg_sources(local_fields, chi3_tensors)` — χ^(3):EEE triple contraction per dipole
- [x] `compute_thg_response(solver, omega_response, dipoles_3omega, chi3_tensors, params_3omega, calc_far_field)` — self-consistent THG solve at λ/3
- [x] `ThgResult`: source_moments, driven_moments, thg_intensity (nm⁹, ∝|E|⁶), optional far-field
- [x] CLI: `enable_thg`, `chi3_symmetry`, `chi3_xxxx`, `chi3_xxyy` in `[nonlinear]`; `thg_spectrum.csv` output
- [x] GUI: THG panel alongside SHG; ThgSpectra tab in results
- [x] 5 unit tests: zero χ^(3), χ_xxxx contraction, isotropic_bulk tensor structure, end-to-end solve, E⁶ scaling

### Deferred to v0.5

- [ ] Full Sommerfeld integral substrate (k∥-dependent Fresnel coefficients)
- [ ] 1D reciprocal Ewald sum
- [ ] DFT polarisability import (VASP WAVEDER, Gaussian)
- [ ] MPI distributed wavelength sweep

## v1.0 — Multi-Method & Cross-Platform

- [ ] Boundary Element Method (BEM) solver
- [ ] T-matrix solver
- [ ] FDTD time-domain solver
- [ ] CUDA-specific compute backend (cuSOLVER/cuBLAS)
- [ ] Metal Performance Shaders backend
- [ ] Cross-platform binary releases
- [ ] Comprehensive test suite with property-based testing
- [ ] Full documentation with tutorials and examples
