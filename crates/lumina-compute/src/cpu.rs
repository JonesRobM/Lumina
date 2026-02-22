//! CPU compute backend using Rayon for shared-memory parallelism.

use ndarray::{Array1, Array2};
use num_complex::Complex64;

use crate::backend::{BackendType, ComputeBackend, ComputeError, DeviceInfo};

/// CPU backend that parallelises work across threads via Rayon.
pub struct CpuBackend {
    num_threads: usize,
}

impl CpuBackend {
    /// Create a new CPU backend using all available threads.
    pub fn new() -> Self {
        Self {
            num_threads: rayon::current_num_threads(),
        }
    }

    /// Create a CPU backend with a specified thread count.
    pub fn with_threads(num_threads: usize) -> Self {
        Self { num_threads }
    }
}

impl Default for CpuBackend {
    fn default() -> Self {
        Self::new()
    }
}

impl ComputeBackend for CpuBackend {
    fn device_info(&self) -> DeviceInfo {
        DeviceInfo {
            name: format!("CPU ({} threads)", self.num_threads),
            backend_type: BackendType::Cpu,
            memory_bytes: None,
            compute_units: Some(self.num_threads),
        }
    }

    fn parallel_matrix_fill(
        &self,
        rows: usize,
        cols: usize,
        fill_fn: &(dyn Fn(usize, usize) -> Complex64 + Send + Sync),
    ) -> Result<Array2<Complex64>, ComputeError> {
        use rayon::prelude::*;

        let data: Vec<Complex64> = (0..rows * cols)
            .into_par_iter()
            .map(|idx| {
                let i = idx / cols;
                let j = idx % cols;
                fill_fn(i, j)
            })
            .collect();

        Array2::from_shape_vec((rows, cols), data)
            .map_err(|e| ComputeError::DeviceError(e.to_string()))
    }

    fn matvec(
        &self,
        matrix: &Array2<Complex64>,
        vector: &Array1<Complex64>,
    ) -> Result<Array1<Complex64>, ComputeError> {
        Ok(matrix.dot(vector))
    }
}
