//! Iterative GMRES solver for large CDA systems.
//!
//! For systems with $N > 1000$ dipoles, direct LU decomposition becomes
//! prohibitively expensive ($O(N^3)$). GMRES (Generalised Minimal Residual
//! method) solves the system iteratively using only matrix-vector products,
//! which can be accelerated via GPU or parallelised across nodes.

use ndarray::Array1;
use num_complex::Complex64;

use super::super::SolverError;

/// Solve the interaction system using GMRES.
///
/// # Arguments
/// * `matrix` - The $3N \times 3N$ interaction matrix (or a closure for matrix-free operation).
/// * `rhs` - The incident field vector (length $3N$).
/// * `tolerance` - Convergence tolerance on the relative residual.
/// * `max_iterations` - Maximum number of GMRES iterations.
///
/// # Returns
/// The solved dipole moment vector $\mathbf{p}$ (length $3N$).
pub fn solve_gmres(
    _matrix: &ndarray::Array2<Complex64>,
    _rhs: &Array1<Complex64>,
    _tolerance: f64,
    _max_iterations: usize,
) -> Result<Array1<Complex64>, SolverError> {
    // TODO: Implement GMRES with Arnoldi iteration.
    // Future: accept a matrix-free operator trait for GPU-accelerated matvec.
    todo!("GMRES iterative solve")
}
