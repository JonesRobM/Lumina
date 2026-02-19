//! Direct linear solver for small CDA systems.
//!
//! Uses LU decomposition via `faer` to solve the $3N \times 3N$ system
//! $\mathbf{A}\mathbf{p} = \mathbf{E}_{\text{inc}}$ exactly.
//!
//! Appropriate for systems with $N \leq 1000$ dipoles (i.e. matrices up to
//! $3000 \times 3000$). For larger systems, use the iterative GMRES solver.

use ndarray::Array1;
use num_complex::Complex64;

use super::super::SolverError;

/// Solve the interaction system using direct LU decomposition.
///
/// # Arguments
/// * `matrix` - The $3N \times 3N$ interaction matrix $\mathbf{A}$.
/// * `rhs` - The incident field vector $\mathbf{E}_{\text{inc}}$ (length $3N$).
///
/// # Returns
/// The solved dipole moment vector $\mathbf{p}$ (length $3N$).
pub fn solve_direct(
    _matrix: &ndarray::Array2<Complex64>,
    _rhs: &Array1<Complex64>,
) -> Result<Array1<Complex64>, SolverError> {
    // TODO: Convert ndarray → faer matrix, perform LU, solve, convert back.
    todo!("Direct LU solve via faer")
}
