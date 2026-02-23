# Architecture Overview

Lumina is organised as a Rust workspace with six crates.

## Crate Dependency Graph

```
lumina-gui ──┐
lumina-cli ──┤
             ├── lumina-core ──┬── lumina-compute
             │                 ├── lumina-geometry
             │                 └── lumina-materials
             └─────────────────┘
```

## Crate Responsibilities

| Crate | Type | Purpose |
|-------|------|---------|
| `lumina-core` | Library | Physics engine: solvers, fields, cross-sections |
| `lumina-compute` | Library | Device abstraction: CPU, GPU, MPI |
| `lumina-geometry` | Library | Shapes, file parsers (OBJ, XYZ), discretisation, transforms |
| `lumina-materials` | Library | Material data providers (J&C, Palik), spline interpolation |
| `lumina-gui` | Binary | Interactive egui dashboard |
| `lumina-cli` | Binary | Headless CLI for HPC/batch runs |

## Key Abstractions

### `OpticalSolver` (lumina-core)

All simulation methods implement this trait. The CDA is the first; BEM, FDTD, and T-matrix will follow.

```rust
pub trait OpticalSolver {
    fn compute_cross_sections(...) -> Result<CrossSections, SolverError>;
    fn solve_dipoles(...)          -> Result<DipoleResponse, SolverError>;
    fn compute_near_field(...)     -> Result<NearFieldMap, SolverError>;
    fn compute_far_field(...)      -> FarFieldMap;
}
```

### `ComputeBackend` (lumina-compute)

Abstracts over execution environments. Physics code calls `parallel_matrix_fill()` and `matvec()` without knowing the device. See [Compute Backends](./compute-backends.md) for full details.

The `CdaSolver` holds an `Arc<dyn ComputeBackend>` and delegates matvec calls to whatever backend is configured (CPU or GPU).

### `MaterialProvider` (lumina-materials)

Returns frequency-dependent \\(\epsilon(\omega)\\) for any material source. Built-in providers include Johnson & Christy metals (Au, Ag, Cu) and Palik dielectrics (TiO\\(_2\\), SiO\\(_2\\)).

## Solver Pipeline

For a single wavelength, the CDA solver pipeline is:

1. **Discretise** — fill the shape with a centred cubic lattice of dipole positions (`lumina-geometry`)
2. **Look up \\(\epsilon(\lambda)\\)** — query the material provider (`lumina-materials`)
3. **Assemble** — build the 3N\\(\times\\)3N interaction matrix using the dyadic Green's function, parallelised with Rayon (`lumina-core`)
4. **Solve** — compute self-consistent dipole moments:
   - N \\(\leq\\) 1000: direct LU via `faer`
   - N \\(>\\) 1000: matrix-free GMRES(m=30) with matvec delegated to the `ComputeBackend`
5. **Post-process** — compute cross-sections, near-field maps, far-field patterns (`lumina-core`)

## Feature Flags

The `gpu` feature propagates through the crate chain:

```
lumina-compute/gpu  →  lumina-core/gpu  →  lumina-gui/gpu
                                        →  lumina-cli/gpu
```

When `gpu` is enabled, `GpuBackend` becomes available and the GUI/CLI can offer GPU acceleration. Without the flag, only `CpuBackend` is compiled.
