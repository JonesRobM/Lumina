# v0.2.2 Design: Persistent GPU Buffers

**Date:** 2026-03-01
**Status:** Approved
**Scope:** `lumina-compute`, `lumina-core`

---

## Problem

The current `GpuBackend::matvec()` implementation allocates 7 GPU buffers and creates 1 bind group per call. For a typical GMRES run (30 iterations × 10 restarts = 300 matvecs per wavelength), this produces ~1800 allocations and ~300 matrix uploads per wavelength — the dominant GPU overhead for our N range.

Additionally, the GUI parallel wavelength sweep shares a `GpuBackend` across Rayon threads with no command queue synchronisation, creating a potential race condition on command submission.

---

## Goals

1. Eliminate per-matvec buffer allocations in GMRES via pre-allocated persistent buffers.
2. Upload the interaction matrix to the GPU once per solve, not once per iteration.
3. Fix thread safety for parallel wavelength sweeps in the GUI and CLI.

---

## Non-Goals

- GPU matrix assembly (Green's function branching is GPU-unfriendly).
- FFT-accelerated matvec (deferred to separate v0.2.2 task).
- Surface averaging (separate v0.2.2 task).
- Async/pipelined GPU submission (complexity not justified for N < 10K).

---

## Design

### 1. New Data Structures

**`GpuState`** (internal to `gpu.rs`, wrapped in `Arc<Mutex<_>>`):
```rust
struct GpuState {
    device: wgpu::Device,
    queue: wgpu::Queue,
    pipeline: wgpu::ComputePipeline,
    bind_group_layout: wgpu::BindGroupLayout,
}
```
Serialises all device/queue access across Rayon threads.

**`GpuSession`** (one per solve, returned to caller):
```rust
pub struct GpuSession {
    matrix_buf:  wgpu::Buffer,   // 9N² × 8 bytes (f32 complex pairs)
    input_buf:   wgpu::Buffer,   // 3N × 8 bytes
    output_buf:  wgpu::Buffer,   // 3N × 8 bytes
    staging_buf: wgpu::Buffer,   // 3N × 8 bytes (MAP_READ)
    params_buf:  wgpu::Buffer,   // 16 bytes (uniform)
    bind_group:  wgpu::BindGroup,
    dim: usize,
    state: Arc<Mutex<GpuState>>,
}
```
Buffers are owned by the session; bind group is pre-built at creation time.

**`GpuBackend`** (simplified):
```rust
pub struct GpuBackend {
    state: Arc<Mutex<GpuState>>,
}
```

### 2. API Changes

**`ComputeBackend` trait** gains:
```rust
fn create_gpu_session(&self, n: usize) -> Option<GpuSession>;
```
- `CpuBackend` → returns `None`
- `GpuBackend` → returns `Some(GpuSession)` with pre-allocated buffers

New companion method on `GpuSession`:
```rust
impl GpuSession {
    pub fn upload_matrix(&mut self, matrix: &Array2<Complex64>);
    pub fn matvec(&mut self, x: ArrayView1<Complex64>) -> Array1<Complex64>;
}
```

The existing `ComputeBackend::matvec(&matrix, x)` is retained for the direct-solve (small N) path where no session overhead is desired.

### 3. Execution Flow

**Session creation** (`GpuBackend::create_gpu_session(n)`):
1. Lock `Mutex<GpuState>`
2. Allocate all 5 buffers at fixed size 3N
3. Build bind group binding all 5 buffers
4. Release lock; return `GpuSession` (buffers owned by session)

**Matrix upload** (`session.upload_matrix(&matrix)`) — called once before GMRES:
1. Lock `Mutex<GpuState>` (queue access)
2. Convert matrix `Complex64` → `[f32; 2]` pairs once
3. `queue.write_buffer(&matrix_buf, 0, data)`
4. Release lock

**Per-iteration matvec** (`session.matvec(x)`):
1. Lock `Mutex<GpuState>`
2. Write `x` into `input_buf` (`queue.write_buffer`)
3. Submit compute dispatch (bind group already set, matrix already on GPU)
4. Copy `output_buf` → `staging_buf`
5. Map staging, poll, read result
6. Release lock; return `Array1<Complex64>`

**`CdaSolver` GMRES path** changes:
```
Before loop:  session = backend.create_gpu_session(3*N)  // Option<GpuSession>
              session.upload_matrix(&matrix)
Per iteration: y = session.matvec(x)                     // replaces backend.matvec(&matrix, x)
After loop:   session dropped → GPU buffers freed
```

### 4. Per-Call Savings

| Operation         | Current | v0.2.2 |
|-------------------|---------|--------|
| Buffer allocations | 7      | 0      |
| Bind group creates | 1      | 0      |
| Matrix writes      | 1      | 0      |
| Vector writes      | 1      | 1      |
| Staging maps       | 1      | 1      |

For 300 matvecs/wavelength: eliminates ~1800 allocations and ~300 matrix uploads.

### 5. Thread Safety

- `GpuSession` is **not** shared between threads — each wavelength's `solve_dipoles()` creates and owns its own session.
- `Arc<Mutex<GpuState>>` is cloned into each session at creation; Mutex serialises device/queue access.
- Rayon threads for different wavelengths queue up on the Mutex — correct serialised GPU access.
- No unsafe code added.

---

## Testing

| Test | Description | Pass Criterion |
|------|-------------|----------------|
| `test_persistent_session_matches_direct` | 6×6 system: session matvec vs CPU matvec | Agree to 1e-5 (f32 floor) |
| `test_mie_sphere_gpu_session` | Gold sphere Mie validation with GPU backend | C_ext within 20% of Mie theory |
| `test_parallel_wavelengths_gpu` | 4 Rayon threads, shared `Arc<GpuBackend>`, small system | No panics; results match CPU |

---

## Files Affected

| File | Change |
|------|--------|
| `crates/lumina-compute/src/gpu.rs` | Refactor `GpuBackend`; add `GpuState`, `GpuSession` |
| `crates/lumina-compute/src/backend.rs` | Add `create_gpu_session()` to trait; update `CpuBackend` stub |
| `crates/lumina-core/src/solver/cda/mod.rs` | Use session in GMRES path |
| `crates/lumina-compute/src/tests/` | Add 3 new tests |
