//! Compute backend trait and device abstraction.
//!
//! The [`ComputeBackend`] trait abstracts over different execution environments
//! (CPU, GPU, distributed) so that the physics code in `lumina-core` remains
//! device-agnostic.

use ndarray::Array2;
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

/// Abstraction over compute backends.
///
/// Physics code in `lumina-core` operates against this trait. Implementations
/// provide device-specific execution for the hot-path operations (matrix
/// assembly, matrix-vector products, linear solves).
pub trait ComputeBackend: Send + Sync {
    /// Return information about the device.
    fn device_info(&self) -> DeviceInfo;

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
