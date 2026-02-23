# Compute Backends

The `lumina-compute` crate provides a `ComputeBackend` trait that abstracts over different execution environments. Physics code in `lumina-core` calls `parallel_matrix_fill()` and `matvec()` without knowing whether the work runs on CPU or GPU.

## Available Backends

| Backend | Feature flag | Status | Description |
|---------|-------------|--------|-------------|
| CPU (Rayon) | `cpu` (default) | Implemented | Shared-memory parallelism via Rayon |
| GPU (wgpu) | `gpu` | Implemented (v0.2.0) | Compute shaders for matrix-vector products |
| Distributed (MPI) | `distributed` | Planned | Multi-node parallelism for HPC clusters |

### CPU (default)

Uses Rayon for shared-memory parallelism. Enabled by default via the `cpu` feature.

- **Matrix assembly:** Off-diagonal Green's tensor blocks computed in parallel via `rayon::par_iter`, then placed into the interaction matrix sequentially.
- **Matvec:** `ndarray` dense matrix-vector product (`matrix.dot(&vector)`).
- **Direct solve:** LU decomposition via `faer` for systems with N \\(\leq\\) 1000 dipoles.

### GPU (wgpu)

Uses wgpu v24 compute shaders for the GMRES matrix-vector product — the dominant cost for large systems (\\(N > 1000\\)).

Enabled via the `gpu` feature flag:

```bash
cargo build --release --features gpu
```

**Architecture:**

1. The interaction matrix is assembled on the CPU (Rayon-parallel).
2. The matrix is uploaded to the GPU once and cached in VRAM.
3. Each GMRES iteration dispatches a WGSL compute shader for the complex matvec.
4. The Arnoldi orthogonalisation and Givens rotations remain on the CPU in f64.

**Key design decisions:**

- **f32 precision on GPU:** wgpu/WGSL does not support f64. The matrix and vectors are converted to `Complex<f32>` (stored as `vec2<f32>`) for the GPU matvec. The conversion overhead is O(N), negligible compared to the O(N\\(^2\\)) matvec.
- **Matrix caching:** The `GpuBackend` caches the uploaded matrix buffer. For a typical GMRES solve (~30–300 iterations), the matrix is uploaded once and reused for every iteration.
- **Automatic fallback:** If GPU initialisation fails (no compatible GPU, missing drivers), the solver falls back to `CpuBackend` with a log warning.

**Precision implications:**

The f32 matvec limits the achievable GMRES residual to approximately \\(10^{-7}\\). For most CDA applications this is more than adequate — the discretisation error (typically 5–25%) dominates. The solver tolerance should be set to \\(\geq 10^{-6}\\) when using GPU acceleration.

**WGSL shader:**

The compute shader (`crates/lumina-compute/src/shaders/matvec.wgsl`) performs one complex dot product per thread (one row of the matrix). Complex multiplication uses the standard `(a+bi)(c+di) = (ac-bd) + (ad+bc)i` formula with `vec2<f32>` storage.

```wgsl
@compute @workgroup_size(256)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let row = gid.x;
    // ... complex dot product for row ...
    output_vec[row] = sum;
}
```

### Distributed (planned)

MPI-based multi-node parallelism for HPC clusters. Block-row distribution of the interaction matrix with hybrid MPI+Rayon within each node.

Enabled via the `distributed` feature flag.

## Backend Selection

### In TOML configuration (CLI)

```toml
[simulation]
backend = "auto"   # "cpu", "gpu", or "auto" (default)
```

- `"auto"` — tries GPU first, falls back to CPU if unavailable
- `"cpu"` — always use CPU (useful for reproducibility or debugging)
- `"gpu"` — require GPU (errors if unavailable)

### In the GUI

The simulation panel includes a **GPU acceleration** checkbox. This is only visible when the binary is built with `--features gpu`.

### Programmatic (Rust API)

```rust
use std::sync::Arc;
use lumina_compute::{ComputeBackend, CpuBackend};
use lumina_core::solver::cda::CdaSolver;

// Default: CpuBackend
let solver = CdaSolver::default();

// Explicit GPU backend
#[cfg(feature = "gpu")]
{
    use lumina_compute::GpuBackend;
    let gpu = Arc::new(GpuBackend::new_blocking().expect("GPU init failed"));
    let solver = CdaSolver::with_backend(gpu);
}
```

## Performance

Benchmark results (NVIDIA RTX 4050 Laptop GPU, release mode):

### Matvec (single A\\(\cdot\\)x)

| N (3N\\(\times\\)3N) | CPU (ms) | GPU (ms) | Speedup |
|----------------------|----------|----------|---------|
| 300 | 0.06 | 1.1 | 0.05\\(\times\\) |
| 1 000 | 1.1 | 6.1 | 0.18\\(\times\\) |
| 3 000 | 16.9 | 55.5 | 0.30\\(\times\\) |

### GMRES solve

| N (3N\\(\times\\)3N) | CPU (ms) | GPU (ms) | Speedup |
|----------------------|----------|----------|---------|
| 300 | 0.5 | 13.1 | 0.04\\(\times\\) |
| 1 000 | 9.1 | 35.7 | 0.26\\(\times\\) |
| 3 000 | 83.8 | 406.5 | 0.21\\(\times\\) |

At N \\(\leq\\) 3 000, the GPU is slower than CPU due to per-call buffer overhead (input/output/staging buffer creation and readback latency). The crossover point is expected at N \\(>\\) 5 000, where the O(N\\(^2\\)) computation dominates over the constant overhead. Persistent buffer reuse is planned for v0.2.1.

## Implementing a New Backend

Implement the `ComputeBackend` trait:

```rust
pub trait ComputeBackend: Send + Sync {
    fn device_info(&self) -> DeviceInfo;
    fn parallel_matrix_fill(
        &self, rows: usize, cols: usize,
        fill_fn: &(dyn Fn(usize, usize) -> Complex64 + Send + Sync),
    ) -> Result<Array2<Complex64>, ComputeError>;
    fn matvec(
        &self, matrix: &Array2<Complex64>, vector: &Array1<Complex64>,
    ) -> Result<Array1<Complex64>, ComputeError>;
    fn dense_solve(
        &self, matrix: &Array2<Complex64>, rhs: &Array1<Complex64>,
    ) -> Result<Array1<Complex64>, ComputeError>;
}
```

The `dense_solve` method has a default implementation that returns `Err(ComputeError::Unavailable)`. Override it if your backend supports dense linear solves.
