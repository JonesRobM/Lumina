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
| `lumina-geometry` | Library | Shapes, file parsers, discretisation, transforms |
| `lumina-materials` | Library | Material data providers, spline interpolation |
| `lumina-gui` | Binary | Interactive egui dashboard |
| `lumina-cli` | Binary | Headless CLI for HPC/batch runs |

## Key Abstractions

### `OpticalSolver` (lumina-core)
All simulation methods implement this trait. The CDA is the first; BEM, FDTD, and T-matrix will follow.

### `ComputeBackend` (lumina-compute)
Abstracts over execution environments. Physics code calls `parallel_matrix_fill()` and `matvec()` without knowing the device.

### `MaterialProvider` (lumina-materials)
Returns frequency-dependent \\(\epsilon(\omega)\\) for any material source.
