//! Compute backend trait and device abstraction.
//!
//! The [`ComputeBackend`] trait abstracts over different execution environments
//! (CPU, GPU, distributed) so that the physics code in `lumina-core` remains
//! device-agnostic.

use ndarray::{Array1, Array2, ArrayView1};
use num_complex::Complex64;
use thiserror::Error;

/// Errors originating from compute backends.
#[derive(Debug, Error)]
pub enum ComputeError {
    #[error("Backend not available: {0}")]
    Unavailable(String),

    #[error("Device error: {0}")]
    DeviceError(String),

    #[error("Out of memory: requested {requested} bytes, available {available}")]
    OutOfMemory { requested: usize, available: usize },
}

/// Describes the capabilities of a compute backend.
#[derive(Debug, Clone)]
pub struct DeviceInfo {
    pub name: String,
    pub backend_type: BackendType,
    pub memory_bytes: Option<usize>,
    pub compute_units: Option<usize>,
}

/// The type of compute backend.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BackendType {
    Cpu,
    Gpu,
    Distributed,
}

/// A persistent GPU solve session with pre-allocated device buffers.
///
/// Obtained via [`ComputeBackend::create_session`]. Upload the interaction
/// matrix once with [`upload_matrix`](MatvecSession::upload_matrix), then call
/// [`matvec`](MatvecSession::matvec) for each GMRES iteration without any
/// further device allocations.
pub trait MatvecSession: Send + Sync {
    /// Upload the interaction matrix to device memory.
    ///
    /// Must be called once before the first [`matvec`](MatvecSession::matvec).
    fn upload_matrix(&self, matrix: &Array2<Complex64>) -> Result<(), ComputeError>;

    /// Perform $\mathbf{y} = \mathbf{A}\mathbf{x}$ using the pre-uploaded matrix.
    fn matvec(&self, x: ArrayView1<Complex64>) -> Result<Array1<Complex64>, ComputeError>;
}

/// Abstraction over compute backends.
///
/// Physics code in `lumina-core` operates against this trait. Implementations
/// provide device-specific execution for the hot-path operations (matrix
/// assembly, matrix-vector products, linear solves).
pub trait ComputeBackend: Send + Sync {
    /// Return information about the device.
    fn device_info(&self) -> DeviceInfo;

    /// Create a persistent solve session with pre-allocated buffers for `n` unknowns.
    ///
    /// GPU backends return a [`MatvecSession`] that eliminates per-matvec buffer
    /// allocations. CPU backends return `None` (no allocation overhead to eliminate).
    fn create_session(&self, _n: usize) -> Option<Box<dyn MatvecSession>> {
        None
    }

    /// Perform a parallel element-wise operation over a matrix.
    ///
    /// This is the primary entry point for parallelising the Green's function
    /// matrix assembly: each (i, j) block can be computed independently.
    fn parallel_matrix_fill(
        &self,
        rows: usize,
        cols: usize,
        fill_fn: &(dyn Fn(usize, usize) -> Complex64 + Send + Sync),
    ) -> Result<Array2<Complex64>, ComputeError>;

    /// Perform a complex matrix-vector product $\mathbf{y} = \mathbf{A}\mathbf{x}$.
    ///
    /// Used by iterative solvers (GMRES). GPU backends can keep the matrix
    /// resident on-device and perform the matvec without host transfer.
    fn matvec(
        &self,
        matrix: &Array2<Complex64>,
        vector: &ndarray::Array1<Complex64>,
    ) -> Result<ndarray::Array1<Complex64>, ComputeError>;

    /// Assemble the interaction matrix on-device and return a ready session.
    ///
    /// GPU backends implement this to eliminate CPU assembly + PCIe upload.
    /// The default returns `None` (CPU backends do not override this).
    ///
    /// Activation guards (enforced by the caller in `cda/mod.rs`):
    /// - `use_fcd == false`
    /// - `lattice.is_none()` (no periodic BC)
    /// - `substrate_runtime.is_none()` (no substrate)
    ///
    /// # Arguments
    /// * `positions`        — Dipole positions in nm, one `[f64; 3]` per dipole.
    /// * `k`                — Wavenumber in the medium (nm⁻¹, always real).
    /// * `diagonal_blocks`  — Pre-computed α⁻¹ per dipole (output of `compute_diagonal_blocks`).
    fn assemble_and_create_session(
        &self,
        positions: &[[f64; 3]],
        k: f64,
        diagonal_blocks: &[[Complex64; 9]],
    ) -> Option<Result<Box<dyn MatvecSession>, ComputeError>> {
        let _ = (positions, k, diagonal_blocks);
        None
    }

    /// Solve a dense linear system $\mathbf{A}\mathbf{x} = \mathbf{b}$.
    ///
    /// For CPU, delegates to `faer` LU. GPU backends may use cuSOLVER or
    /// equivalent. This is reserved for future GPU-resident solves.
    fn dense_solve(
        &self,
        matrix: &Array2<Complex64>,
        rhs: &ndarray::Array1<Complex64>,
    ) -> Result<ndarray::Array1<Complex64>, ComputeError> {
        // Default: not supported. Implementations override this.
        let _ = (matrix, rhs);
        Err(ComputeError::Unavailable(
            "Dense solve not implemented for this backend".into(),
        ))
    }
}
