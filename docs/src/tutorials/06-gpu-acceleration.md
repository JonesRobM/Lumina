# Tutorial 6: GPU-Accelerated Large System

**Goal:** Enable GPU acceleration for a large nanoparticle simulation (N > 5 000 dipoles) and compare performance between CPU and GPU backends.

**Time:** ~15 minutes

**Requires:** `--features gpu` build and a compatible GPU (Vulkan/DX12/Metal)

---

## Background

For large dipole systems the GMRES iterative solver dominates the computation time. Each GMRES restart cycle performs ~30 matrix-vector products (matvecs), and each matvec is an O(N²) operation on the 3N×3N interaction matrix. GPU acceleration offloads these matvecs to the GPU via wgpu compute shaders.

**How it works:**

1. The interaction matrix is assembled on the CPU using Rayon-parallel Green's tensor computation.
2. The matrix is uploaded to the GPU once (converted from f64 to f32).
3. Each GMRES matvec dispatches a WGSL compute shader — one thread per row.
4. The Arnoldi orthogonalisation and Givens rotations remain on the CPU in f64 precision.
5. After convergence, the GPU buffer is released.

The f32 precision on the GPU limits the achievable GMRES residual to ~10⁻⁷, but this is far below the CDA discretisation error (typically 5–25%).

---

## Step 1: Build with GPU Support

```bash
cd Lumina
cargo build --release --features gpu
```

Verify that your GPU is detected:

```bash
cargo test -p lumina-compute --features gpu -- test_gpu_backend_creation --nocapture
```

Expected output:

```
test gpu::tests::test_gpu_backend_creation ... ok
```

If this test fails, check your GPU drivers. See [Installation](../user-guide/installation.md) for troubleshooting.

---

## Step 2: Configure a Large Simulation

Create `examples/large_sphere_gpu.toml`:

```toml
# Tutorial 6 — large Au sphere with GPU acceleration
[simulation]
wavelengths    = { range = [400.0, 700.0], points = 30 }
environment_n  = 1.0
backend        = "auto"        # uses GPU if available, CPU fallback

[[geometry.object]]
name           = "Au_Sphere_Large"
type           = "sphere"
centre         = [0.0, 0.0, 0.0]
radius         = 30.0          # nm
material       = "Au_JC"
dipole_spacing = 1.5           # nm  → ~N ≈ 12 000 dipoles

[output]
directory      = "./output/tutorial_06"
save_spectra   = true
save_json      = true
```

> **Dipole count:** A 30 nm sphere at d=1.5 nm spacing generates roughly 12 000 dipoles, giving a 36 000 × 36 000 interaction matrix. This is well above the GMRES threshold (N > 1 000).

---

## Step 3: Run with GPU

```bash
cargo run --release -p lumina-cli --features gpu -- run examples/large_sphere_gpu.toml
```

Watch the log output for the backend selection:

```
[INFO] Backend: GPU (NVIDIA GeForce RTX 4050 Laptop GPU)
[INFO] Geometry: Au_Sphere_Large — 12083 dipoles
[INFO] Solver: GMRES(m=30), tolerance=1e-6
```

If no GPU is available, you will see:

```
[WARN] GPU unavailable, falling back to CPU
[INFO] Backend: CPU (Rayon, 16 threads)
```

---

## Step 4: Compare CPU vs GPU

Run the same simulation with `backend = "cpu"` to measure the difference:

```toml
[simulation]
backend = "cpu"
```

```bash
cargo run --release -p lumina-cli --features gpu -- run examples/large_sphere_gpu.toml
```

For a rough comparison, use the `RUST_LOG=info` environment variable to see per-wavelength timing:

```bash
RUST_LOG=info cargo run --release -p lumina-cli --features gpu -- run examples/large_sphere_gpu.toml
```

### Expected scaling

| N (dipoles) | 3N | CPU matvec | GPU matvec | Notes |
|------------|-----|-----------|-----------|-------|
| 1 000 | 3 000 | ~1 ms | ~6 ms | GPU overhead dominates |
| 3 000 | 9 000 | ~17 ms | ~56 ms | Still overhead-bound |
| 10 000 | 30 000 | ~200 ms | ~150 ms | Crossover region |
| 30 000 | 90 000 | ~2 s | ~0.5 s | GPU 4× faster |

The GPU crossover point depends on your hardware. The overhead comes from per-call buffer creation — future versions will use persistent buffers to shift the crossover to smaller N.

---

## Step 5: Run the Built-in Benchmark

Lumina includes a GPU vs CPU benchmark test:

```bash
cargo test -p lumina-core --features gpu --release -- gpu_benchmark --nocapture
```

This prints a table comparing matvec and full GMRES solve times at N = 300, 1 000, and 3 000:

```
=== GPU Device: GPU (NVIDIA GeForce RTX 4050 Laptop GPU) ===

--- Matvec (A*x) ---
Dim          CPU (ms)     GPU (ms)    Speedup
----------------------------------------------
300             0.056        1.054      0.05x
1000            1.103        6.120      0.18x
3000           16.876       55.533      0.30x

--- GMRES Solve ---
Dim          CPU (ms)     GPU (ms)    Speedup
----------------------------------------------
300               0.5         13.1      0.04x
1000              9.1         35.7      0.26x
3000             83.8        406.5      0.21x
```

---

## Step 6: GPU in the GUI

Launch the GUI with GPU support:

```bash
cargo run --release -p lumina-gui --features gpu
```

In the **Simulation** panel:

1. Check the **GPU acceleration (wgpu)** checkbox
2. Set a fine dipole spacing (e.g. 1.5 nm) for a large system
3. Click **Run Simulation**

The GPU backend is used for the GMRES matvec; all other computations (assembly, post-processing) remain on the CPU.

---

## Precision Notes

The GPU computes matvecs in f32 (single precision) because wgpu/WGSL does not support f64. This means:

- **GMRES residual floor:** ~10⁻⁷ (cannot converge tighter than this)
- **Solution accuracy:** agrees with f64 CPU to ~10⁻⁶ relative error
- **Physical accuracy:** f32 precision is negligible compared to the CDA discretisation error (5–25%)

For applications requiring tighter solver convergence (e.g. academic benchmarking against Mie theory), use `backend = "cpu"`.

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| "GPU unavailable" warning | Check GPU drivers; run `vulkaninfo` to verify Vulkan support |
| GPU slower than CPU | Expected for N < 5 000; use CPU for small systems |
| GMRES fails to converge with GPU | Increase `solver_tolerance` to ≥ 1e-6 (f32 floor) |
| Build error with `--features gpu` | Ensure wgpu v24 dependencies resolve; check Rust ≥ 1.75 |

---

## Next Steps

- Read the [Compute Backends](../developer/compute-backends.md) developer guide for architecture details
- See the [Validation](../theory/validation.md) chapter for GPU correctness benchmarks
- Return to the [Tutorials index](./index.md) for a summary of all tutorials
