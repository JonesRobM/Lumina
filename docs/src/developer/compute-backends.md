# Compute Backends

The `lumina-compute` crate provides a `ComputeBackend` trait that abstracts over different execution environments.

## Available Backends

### CPU (default)
Uses Rayon for shared-memory parallelism. Enabled by default via the `cpu` feature.

### GPU (planned)
Uses wgpu compute shaders for matrix assembly and matrix-vector products. Dense linear solves fall back to CPU via `faer` in v0.1; future versions will use platform-specific solvers (cuSOLVER, Metal Performance Shaders).

Enabled via the `gpu` feature flag.

### Distributed (planned)
MPI-based multi-node parallelism for HPC clusters. Block-row distribution of the interaction matrix with hybrid MPI+Rayon within each node.

Enabled via the `distributed` feature flag.

## Implementing a New Backend

Implement the `ComputeBackend` trait:

```rust
pub trait ComputeBackend: Send + Sync {
    fn device_info(&self) -> DeviceInfo;
    fn parallel_matrix_fill(...) -> Result<Array2<Complex64>, ComputeError>;
    fn matvec(...) -> Result<Array1<Complex64>, ComputeError>;
    fn dense_solve(...) -> Result<Array1<Complex64>, ComputeError>;
}
```
