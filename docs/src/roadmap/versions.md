# Version Roadmap

## v0.1.0 -- Linear Solver Foundation

- [x] Workspace structure with six crates
- [x] Dyadic Green's function implementation
- [x] Interaction matrix assembly with Rayon parallelism
- [x] Direct solver (LU via `faer`) for N <= 1000
- [x] Clausius-Mossotti polarisability with radiative corrections (RRCM + LDR)
- [x] Extinction, absorption, scattering cross-section computation
- [x] Johnson & Christy Au data (43 points, 188-892 nm) with cubic spline interpolation
- [x] Geometry primitives (sphere, cuboid, ellipsoid)
- [x] Discretisation: centred cubic lattice from primitives
- [x] GUI dashboard with geometry, materials, simulation, and results panels
- [x] CLI with TOML configuration
- [x] Mie theory validation (10 nm Au sphere, dielectric spheres)
- [x] mdBook documentation with validation results

## v0.1.1 -- Accurate & Scalable Linear Solver

- [x] GMRES(m) iterative solver with restart for N > 1000
- [x] Filtered Coupled Dipole (FCD) Green's function (IGT method)
- [x] FCD wired as default in CdaSolver
- [x] Cylinder containment and bounding box
- [x] Helix containment and bounding box
- [x] GMRES vs direct LU validation (machine-precision agreement)
- [x] FCD gold sphere full-spectrum validation (420-800 nm)
- [x] Cylinder and helix discretisation tests

### v0.1.1 Validation Summary

- **GMRES**: Agrees with direct LU to ~10^-13 relative error
- **FCD interband (420-510 nm)**: 13-24% error (improved from 22-28% with point-dipole)
- **FCD Drude (>550 nm)**: Still >80% error -- known limitation at d=2nm
- **Cylinder/Helix**: Correct discretisation with expected dipole counts

## v0.1.2 -- Plotting & Visualisation

- [x] egui_plot integration for interactive spectra plots (C_ext, C_abs, C_sca vs lambda)
- [x] eps(lambda) dielectric function plots in the materials panel (Au, Ag, Cu with correct data)
- [x] Near-field |E|^2 heatmap in the results panel
- [x] 2D dipole lattice scatter preview (xy/xz/yz projection)
- [x] Shape parameter UI for all primitives (sphere, cylinder, cuboid, ellipsoid, helix)
- [x] Helix added to GUI geometry panel (v0.1.1 primitive now exposed in GUI)
- [x] Simulation wired to actual material and geometry panel selections (not hardcoded)

## v0.1.3 -- I/O & Export

- [x] Far-field radiation pattern computation (`compute_far_field` in `fields.rs`)
- [x] Far-field polar plot (E-plane / H-plane cuts) in the GUI results panel
- [x] Circular dichroism DeltaC_ext computation (`compute_circular_dichroism` in `fields.rs`)
- [x] Incident polarisation selector in GUI (x-pol / y-pol / circular)
- [x] CD checkbox in GUI simulation panel; CD line on spectra plot
- [x] File dialog for CSV export (via `rfd`)
- [x] JSON export option for simulation results (GUI + CLI)
- [x] Metadata header in CSV (version, object list, material, spacing)
- [x] Palik TiO2 and SiO2 dielectric data (300-1000 nm) in `lumina-materials`
- [x] TiO2 / SiO2 available in GUI materials panel and CLI (`TiO2_Palik`, `SiO2_Palik`)
- [x] Ag_JC and Cu_JC fully wired in CLI material resolver
- [x] CLI geometry expanded: cylinder, cuboid, ellipsoid, helix in TOML configs
- [x] Near-field export to CSV (`save_near_field = true` in TOML output config)
- [x] .xyz file import in GUI geometry panel (Angstromngstroem -> nm conversion)
- [x] OBJ mesh parser (vertex + face -> volume-filling dipole lattice via ray-casting)

## v0.2.0 -- GPU Compute Engine

- [x] Matrix-free GMRES interface (matvec closure instead of matrix reference)
- [x] ComputeBackend trait wired into CdaSolver (CPU default, GPU optional)
- [x] Rayon-parallel interaction matrix assembly (off-diagonal Green's tensor blocks)
- [x] wgpu WGSL compute shader for complex matrix-vector product (f32)
- [x] GpuBackend with matrix caching (upload once, reuse for ~300 GMRES iterations)
- [x] GPU feature forwarding through crate chain (lumina-compute -> core -> gui/cli)
- [x] GPU toggle in GUI simulation panel
- [x] Backend selection in CLI TOML config (`backend = "auto" | "cpu" | "gpu"`)
- [x] Automatic CPU fallback when GPU unavailable
- [x] GPU vs CPU benchmark harness (matvec + GMRES, N = 300-3000)
- [x] All 85 tests passing (81 existing + 4 GPU)

### v0.2.0 Performance Notes

- GPU matvec is correct (agrees with CPU to ~1e-6, the f32 precision limit)
- At N <= 3000, GPU is slower than CPU due to per-call buffer overhead
- Speedup expected at N > 5000 with persistent buffer reuse (v0.2.1)

## v0.2.1 -- Performance & XYZ Workflow

### Performance

- [x] Parallel wavelength sweep in GUI (rayon `par_iter` -- matches CLI performance)
- [x] Parallel wavelength sweep in CLI (rayon `par_iter`, `AtomicUsize` progress)
- [x] Stack-allocated Green's tensors (`[[Complex64; 3]; 3]` -- eliminates ~N^2 heap allocations per wavelength)
- [x] Direct-write matrix assembly (raw pointer writes to disjoint row groups -- removes intermediate `Vec` allocation)

### XYZ Import Workflow

- [x] Auto-detected nearest-neighbour distance for XYZ imports (correct FCD cell size and CM volume)
- [x] Directory-scoped file browser in GUI geometry panel (pick folder, then select from dropdown)
- [x] Example Au Mackay icosahedron (k=6, 923 atoms) supplied in `examples/au_icosahedron_k6.xyz`
- [x] Mackay icosahedron test suite in `lumina-geometry` (atom counts, symmetry, distances)

### GUI Enhancements

- [x] Debug output toggle in simulation panel (per-wavelength timing, matrix size, eps values, solver method)
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

## v0.2.2 -- Persistent Buffers & Scene Assembly

### Persistent GPU Session Buffers

- [x] `MatvecSession` trait: `upload_matrix` + `matvec` -- matrix uploaded once per wavelength
- [x] `GpuSession`: 5 pre-allocated WGSL buffers + pre-built bind group, reused across all GMRES iterations
- [x] `GpuBackend` holds `Arc<Mutex<GpuState>>` -- thread-safe across Rayon wavelength workers
- [x] Per-matvec overhead eliminated: 7 allocations -> 0, 300 matrix writes -> 1 per wavelength

### SceneSpec -- Complex Assembly Builder

- [x] `SceneSpec` in `lumina-geometry`: shared GUI+CLI scene description, serialisable to/from TOML
- [x] `[[geometry.object]]` arrays: each entry is an `ObjectSpec` with name, shape, material, spacing, transform
- [x] `ObjectTransform`: position (nm), rotation\_deg (Euler XYZ), scale factor
- [x] CoreShell geometry: concentric shell layers, each with its own material
- [x] `SceneSpec::build_geometry()` -> flat `SceneDipoles { positions, materials, spacings, hints }`
- [x] GUI: `ScenePanel` replaces old `GeometryPanel` + `MaterialsPanel`
- [x] CLI: `JobConfig.geometry: SceneSpec`; per-dipole material resolution in runner

### v0.2.2 Performance Notes

- GPU session buffers eliminate all per-iteration allocation overhead on the GPU path
- At N = 1000-3000, GPU path now competitive with CPU; speedup expected at N > 5000

## v0.3.0 -- Periodic Structures & Extended Geometry

### Periodic Structures

- [x] `LatticeSpec` in `lumina-geometry`: 1D chain or 2D planar lattice specification
- [x] Direct real-space lattice sum in `EwaldGreens` (finite shell cutoff, exact for small unit cells)
- [x] `SimulationParams.k_bloch: [f64; 3]` -- Bloch wavevector for oblique incidence
- [x] New result types: `BlochCrossSections`, `DispersionMap` for band structure sweeps

### Extended Geometry

- [x] CoreShell (concentric ellipsoidal shells, arbitrary number of layers)
- [x] Anisotropic ellipsoid dipoles: per-dipole 3x3 polarisability tensor from depolarisation factors
- [x] Depolarisation factors via 16-point Gauss-Legendre quadrature (`ellipsoid_depol_factors`)
- [x] `SubstrateSpec`: image-dipole substrate (quasi-static Fresnel reflection, flat interface below the structure)

### Other

- [x] Palik dielectric data (TiO2, SiO2) -- shipped in v0.1.3
- [ ] Extended Palik library (Al2O3, Si, GaAs, etc.)

## v0.4.0 -- Nonlinear Optics & Accelerated Periodic Solver

### Ewald Acceleration

- [x] `EwaldGreens` upgraded from direct sum to full Ewald split
- [x] Real-space sum: erfc-damped (`erfc_real(eta*|R|)` multiplied per image) -- exponential convergence, <= 5 shells sufficient
- [x] Reciprocal-space sum (2D lattices): spectral G-sum `(i/2A)*sum_G exp(iQ*rho)*exp(-q|z|)/q * M(Q,q)`
- [x] `csqrt_positive_imag`: Im>=0 branch selection for evanescent vs. propagating wave convention
- [x] Combined sum: `evaluate()` = `real_space() + recip_space()` (API unchanged)

### Retarded Substrate

- [x] `fresnel_reflection_amplitude(eps_sub, eps_env)` -- retarded normal-incidence reflection coefficient r = (n_sub - n_env)/(n_sub + n_env)
- [x] `substrate_reflection_factor(..., use_retarded)` dispatcher
- [x] `SubstrateSpec.use_retarded: bool` (default `true`, backward-compatible TOML)
- [x] GUI and CLI updated to use retarded coefficient by default

### Second-Harmonic Generation (SHG)

- [x] `Chi2Tensor` in `lumina-core::types`: rank-3 [Complex64; 27], `zero()`, `isotropic_surface(chi_zzz, chi_zxx)`, `contract()`
- [x] `compute_shg_sources(local_fields, chi2_tensors)` -- chi^2:EE source term per dipole
- [x] `compute_shg_response(solver, omega_response, dipoles_2omega, chi2_tensors, params_2omega, calc_far_field)` -- full self-consistent SHG solve
- [x] `ShgResult`: source_moments, driven_moments, shg_intensity (nm^6, prop to|E|^4), optional far-field pattern
- [x] `CdaSolver::cross_sections_from_response` -- avoids double-solve when SHG is also computed
- [x] CLI: `[nonlinear]` TOML section; `shg_spectrum.csv` output
- [x] GUI: SHG panel with symmetry selector, chi DragValue inputs, SHG spectra tab in results
- [x] 5 unit tests covering zero chi^2, chi_xxx contraction, Cinfv tensor structure, end-to-end solve, E^4 scaling

### Third-Harmonic Generation (THG)

- [x] `Chi3Tensor` in `lumina-core::types`: rank-4 [Complex64; 81], index a*27+b*9+c*3+d, `zero()`, `isotropic_bulk(chi_xxxx, chi_xxyy)`, `contract()`
- [x] `compute_thg_sources(local_fields, chi3_tensors)` -- chi^(3):EEE triple contraction per dipole
- [x] `compute_thg_response(solver, omega_response, dipoles_3omega, chi3_tensors, params_3omega, calc_far_field)` -- self-consistent THG solve at lambda/3
- [x] `ThgResult`: source_moments, driven_moments, thg_intensity (nm^9, prop to|E|^6), optional far-field
- [x] CLI: `enable_thg`, `chi3_symmetry`, `chi3_xxxx`, `chi3_xxyy` in `[nonlinear]`; `thg_spectrum.csv` output
- [x] GUI: THG panel alongside SHG; ThgSpectra tab in results
- [x] 5 unit tests: zero chi^(3), chi_xxxx contraction, isotropic_bulk tensor structure, end-to-end solve, E^6 scaling

### Deferred to v0.5

- [ ] Full Sommerfeld integral substrate (k||-dependent Fresnel coefficients)
- [x] 1D Ewald direct sum fix (v0.4.1) -- Bessel spectral correction still deferred
- [ ] DFT polarisability import (VASP WAVEDER, Gaussian)
- [ ] MPI distributed wavelength sweep

## v0.4.1 -- Memory & Periodic Fixes

### On-the-fly Matvec (CPU GMRES path)

- [x] `assembly::matvec_on_the_fly` -- computes A*x by evaluating G(r_i,r_j) on-the-fly per GMRES iteration; no matrix assembled
- [x] Peak memory drops from O(N^2) to O(N) (Krylov vectors only, ~= 4 MB for N = 3 000)
- [x] `CdaSolver.matrix_memory_budget: u64` (default 2 GiB) -- adaptive strategy: cached matrix if it fits, on-the-fly otherwise
- [x] `matrix_memory_gib` exposed in CLI TOML `[simulation]`

### GPU Memory Budget Integration

- [x] Budget gate applied **before** GPU path -- large systems skip GPU assembly entirely (previously: silent O(N^2) CPU alloc before GPU attempt)
- [x] Restructured `solve_dipoles` branch: direct LU -> (matrix fits budget? -> GPU or CPU cached) -> on-the-fly
- [x] Memory-aware wavelength thread pool cap: `min(cpus, 8 GiB / matrix_bytes)` -- prevents OOM on parallel wavelength sweeps
- [x] Applied in GUI (`app.rs`) and CLI (`runner.rs`) including SHG and THG sweeps

### 1D Ewald Direct Sum Fix

- [x] **Bug fixed**: eta = sqrt(pi/|**a**1|) zeroed all real-space terms (erfc ~= 10^-109 at nearest image) -- 1D periodic G was silently ~= 0
- [x] 1D chains now use a direct lattice sum (no erfc damping) with >= 20 shells (~2 % accuracy at optical wavelengths)
- [x] 2D lattices unaffected -- erfc Ewald still used for 2D
- [x] New tests: `test_1d_chain_finite` strengthened (asserts max_mag > 1e-12), `test_1d_chain_reciprocity` added
- [x] 1D Bessel spectral sum (K0/K1 for complex argument) documented as deferred

### Extended Palik Materials

- [x] `PalikMaterial::al2o3()` -- sapphire (Al2O3) ordinary ray, 300-1000 nm, n ~= 1.76-1.81
- [x] `PalikMaterial::silicon()` -- crystalline Si at 300 K, 400-1100 nm (starts at 400 nm to avoid UV gap structure)
- [x] `PalikMaterial::gaas()` -- GaAs at 300 K, 300-1000 nm, band edge ~870 nm
- [x] Identifiers: `Al2O3_Palik`, `Si_Palik`, `GaAs_Palik` in GUI and CLI

### Documentation

- [x] `docs/src/theory/periodic.md` -- 1D convergence table, Bessel deferral note, updated status/limitations table
- [x] `docs/src/developer/compute-backends.md` -- memory-aware dispatch section, solver path flowchart, memory scaling diagram
- [x] `docs/images/solver_dispatch.svg` -- flowchart of `solve_dipoles` path selection
- [x] `docs/images/memory_scaling.svg` -- log-log memory vs N for cached vs on-the-fly

### Memory Scaling After v0.4.1

| N dipoles | Matrix size | CPU on-the-fly peak | GPU budget cap (8 GiB pool) |
|-----------|------------|---------------------|------------------------------|
| 1 000 | 144 MB | -- (direct LU) | -- (direct LU) |
| 3 000 | 1.3 GB | **~4 MB** | 6 parallel wavelength threads |
| 5 000 | 3.6 GB | **~11 MB** | 2 parallel wavelength threads |
| 10 000 | 14.4 GB | **~43 MB** | 1 thread at a time |

## v0.5.0 -- Spectral Ewald, FFT Matvec & Sommerfeld Substrate

### 1D Bessel Spectral Ewald Sum

- [x] `recip_space_1d` in `EwaldGreens`: spectral reciprocal sum for 1D chains using K0/K1 Bessel functions at complex argument
- [x] Branch-cut handling: K0(i*alpha*rho) via (pi/2)(-Y0+iJ0), K1(i*alpha*rho) via -(pi/2)(J1+iY1)
- [x] Polynomial Bessel approximations (A&S) for Y0, Y1, J0, J1, J2
- [x] Reciprocity test for 1D chain reciprocal sum

### FFT Block-Toeplitz Matvec

- [x] `FftMatvecPlan` in `lumina-core`: circulant embedding in (2Nxx2Nyx2Nz), O(N log N) matvec
- [x] `detect_regular_grid`: detects full-occupancy cubic dipole grids (N = Nx*Ny*Nz)
- [x] Activated automatically when grid is regular and no lattice/substrate/FCD is in use
- [x] Pre-built FFT plans via `rustfft`; circulant displacement uses `rem_euclid`

### GPU Matrix Assembly Shader

- [x] `matrix_assemble.wgsl`: `@workgroup_size(16,16)` WGSL shader fills off-diagonal Green's tensor blocks on GPU
- [x] `assemble_and_create_session` in `GpuBackend`: full GPU assembly + CPU diagonal patch via `write_buffer`
- [x] `compute_diagonal_blocks` in `assembly.rs`: per-dipole 3x3 inversion for radiative correction
- [x] Dispatch priority: FFT -> GPU assembly -> CPU assembly (all within the no-FCD, no-lattice, no-substrate branch)
- [x] GPU state mutex held in single scope to prevent interleave under parallel wavelength sweep

### Sommerfeld Substrate

- [x] `SommerfeldGreens` in `lumina-core`: G7K15 adaptive quadrature along Sommerfeld path, 3 segments per light line
- [x] Bessel J0/J1/J2 for the lateral spectral integrals
- [x] `SubstrateRuntime` enum: `Fresnel { z_interface, delta_eps }` | `Sommerfeld { z_interface, greens }`
- [x] `SubstrateSpec.use_sommerfeld: bool` (default `false`) -- opt-in at CLI and GUI
- [x] Sommerfeld path is free-space only, with no lattice sum; `SubstrateRuntime::from_spec` constructs
      the correct variant per wavelength. As of v0.7.0, combining it with a periodic lattice is
      refused rather than silently approximated.

## v0.5.1 -- AI Chat Panel

### New crate: `lumina-chat`

- [x] `ChatConfig`: provider name, base URL (no trailing slash), model, optional API key; serde round-trips to TOML
- [x] Built-in presets: Ollama (`http://localhost:11434/v1`), Groq, OpenRouter
- [x] `load()` / `save()` to `<config_dir>/lumina/chat.toml`; silent fallback to defaults on any error
- [x] `Role` / `Message` types with `serde(rename_all = "lowercase")` for OpenAI wire format
- [x] `ChatEvent`: `Token(String)`, `Done`, `Error(String)`
- [x] `ChatClient`: background OS thread + `current_thread` tokio runtime; `std::sync::mpsc` task/event channels
- [x] SSE streaming via `reqwest 0.12`; `parse_data_line` extracts `choices[0].delta.content`
- [x] `ContextSnapshot`: clones scene + spectra into a markdown message for the "Attach context" action
- [x] 11 unit tests (config round-trip, presets, SSE parsing, context message format)

### GUI changes (`lumina-gui`)

- [x] `Panel::Chat` added; **Chat** entry in left nav
- [x] `ChatPanel`: Settings tab (provider ComboBox, URL/model/key fields, Save) + Chat tab (streaming message history)
- [x] "Attach context" button injects `ContextSnapshot` as a user message and triggers `client.send`
- [x] Auto-scroll to bottom on new tokens; spinner on last message while streaming
- [x] `LuminaApp::attach_context()` accesses `scene_panel` and `chat_panel` as disjoint field projections (no borrow conflict)
- [x] `egui::Context` injected lazily on first `update()` frame (no eframe dependency at `ChatClient::new`)

## v0.6.0 -- Surface Averaging

Recorded retrospectively: this work shipped but was never given a version entry.

### Geometry (`lumina-geometry`)

- [x] `fill_fraction_mc`: 5x5x5 sub-sampling of each cubic cell, with an 8-corner fast path
- [x] `LatticePoint.fill_fraction`, populated by `discretise_primitive` and `discretise_mesh`
- [x] `SceneDipoles.fill_fractions`, parallel to `positions`

### Physics

- [x] Surface-averaged polarisability alpha_SA = f * alpha_RRCM, applied per dipole
- [x] Polarisability scaling rather than eps-mixing: eps-mixing produces a spurious
      Clausius-Mossotti resonance for Drude metals at fill fractions f ~ 0.2-0.3

### Notes

- Improvement is asymptotic in the dipole spacing. For metals at a single coarse
  spacing it is not guaranteed to be monotonic.

---

## v0.7.0 -- Correctness Pass: FCD Cell Size, 1D Ewald, Solver Performance

A correctness and performance release. Three changes alter numerical results; see
CHANGELOG.md for the migration notes.

### Corrected physics

- [x] FCD cell size is now a per-dipole property (`Dipole.cell_size`) taken from the
      geometry. The CLI previously used a fixed 3.0 nm cell regardless of the job's
      dipole spacing, shifting C_ext by 4-10 per cent for other spacings.
- [x] 1D Ewald spectral sum: removed a spurious 1/rho from the cross coefficient. The
      term derives from d_x d_alpha K0 = (i g_m)(-kappa K1(kappa rho) rho_hat_alpha)
      and carries a bare kappa*K1. Off-axis pairs were wrong by a factor of rho in nm.
- [x] 1D Ewald mode count now scales with a1/rho. A rho-independent count
      under-resolved small rho (about 1e-3 relative error at rho/a1 = 0.14).
- [x] Collinear pairs in a 1D chain are summed directly in real space. The Bessel
      representation is not valid at rho = 0, and nothing covered that case, so
      on-axis chains previously felt no lattice at all.
- [x] A Sommerfeld substrate combined with a periodic lattice is now refused with
      `SolverError::UnsupportedConfiguration` instead of mixing a periodic direct
      term with a non-periodic substrate term.

### Performance

- [x] Direct LU is used up to N = 3000 (previously 1000) and is gated on the matrix
      memory budget. Measured 8.2x to 13.5x faster than GMRES between N = 1021 and
      N = 1791, where the iterative path had been slower.
- [x] Block-Jacobi preconditioning on every GMRES path, including SHG and THG.
      A_ii = alpha^-1, so M^-1 = diag(alpha). This helps only when polarisability
      varies between dipoles; for a single material it is a scalar and has no effect.
- [x] The O(N log N) FFT matvec is reachable from a job file via `use_fcd = false`,
      measured 9.3x faster on a 1331-dipole cube.

### Structure

- [x] `lumina_core::build_dipoles`: one dipole-construction path, replacing the same
      loop duplicated at 11 call sites across the CLI and GUI
- [x] `SceneDipoles::uniform_spacing()`
- [x] `use_fcd` and `iterative_threshold` exposed in the job schema
- [x] `SolverError` is `#[non_exhaustive]`; `solver::cda::assembly` is now crate-private

### Housekeeping

- [x] Workspace version corrected to match the shipped release. It had read 0.3.0,
      and was stamped into every exported data file.
- [x] Bundled examples and the quick-start and geometry guides updated to the current
      scene schema. None of the three examples parsed before this release.
- [x] Source, documentation and diagrams are ASCII-only.


## v1.0 -- Multi-Method & Cross-Platform

- [ ] Boundary Element Method (BEM) solver
- [ ] T-matrix solver
- [ ] FDTD time-domain solver
- [ ] CUDA-specific compute backend (cuSOLVER/cuBLAS)
- [ ] Metal Performance Shaders backend
- [ ] Cross-platform binary releases
- [ ] Comprehensive test suite with property-based testing
- [ ] Full documentation with tutorials and examples
