//! Direct linear solver for small CDA systems.
//!
//! Uses LU decomposition via `faer` to solve the $3N \times 3N$ system
//! $\mathbf{A}\mathbf{p} = \mathbf{E}_{\text{inc}}$ exactly.
//!
//! Appropriate for systems with $N \leq 1000$ dipoles (i.e. matrices up to
//! $3000 \times 3000$). For larger systems, use the iterative GMRES solver.

use faer::linalg::solvers::SpSolver;
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
    matrix: &ndarray::Array2<Complex64>,
    rhs: &Array1<Complex64>,
) -> Result<Array1<Complex64>, SolverError> {
    let dim = matrix.nrows();
    assert_eq!(dim, matrix.ncols(), "Matrix must be square");
    assert_eq!(dim, rhs.len(), "RHS length must match matrix dimension");

    // Convert ndarray to faer Mat<c64>
    let faer_mat = faer::Mat::<faer::complex_native::c64>::from_fn(dim, dim, |i, j| {
        let c = matrix[[i, j]];
        faer::complex_native::c64::new(c.re, c.im)
    });

    // Convert RHS to faer Col<c64>
    let faer_rhs = faer::Col::<faer::complex_native::c64>::from_fn(dim, |i| {
        let c = rhs[i];
        faer::complex_native::c64::new(c.re, c.im)
    });

    // LU decomposition with partial pivoting
    let lu = faer_mat.partial_piv_lu();

    // Solve Ax = b
    let faer_sol = lu.solve(&faer_rhs);

    // Convert back to ndarray
    let solution = Array1::from_vec(
        (0..dim)
            .map(|i| {
                let c = faer_sol[i];
                Complex64::new(c.re, c.im)
            })
            .collect(),
    );

    Ok(solution)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    #[test]
    fn test_solve_identity_system() {
        // Ax = b where A = I, so x = b
        let dim = 6;
        let mut matrix = ndarray::Array2::<Complex64>::zeros((dim, dim));
        for i in 0..dim {
            matrix[[i, i]] = Complex64::from(1.0);
        }
        let rhs = Array1::from_vec(
            (0..dim)
                .map(|i| Complex64::new(i as f64, 0.0))
                .collect(),
        );

        let sol = solve_direct(&matrix, &rhs).unwrap();
        for i in 0..dim {
            assert!((sol[i] - rhs[i]).norm() < 1e-12);
        }
    }

    #[test]
    fn test_solve_complex_system() {
        // 2x2 complex system
        let matrix = ndarray::Array2::from_shape_vec(
            (2, 2),
            vec![
                Complex64::new(1.0, 1.0),
                Complex64::new(2.0, 0.0),
                Complex64::new(0.0, 1.0),
                Complex64::new(3.0, -1.0),
            ],
        )
        .unwrap();

        let rhs = array![Complex64::new(5.0, 1.0), Complex64::new(4.0, 2.0)];

        let sol = solve_direct(&matrix, &rhs).unwrap();

        // Verify: A * sol should equal rhs
        let check = matrix.dot(&sol);
        for i in 0..2 {
            assert!(
                (check[i] - rhs[i]).norm() < 1e-10,
                "Mismatch at {}: got {:?}, expected {:?}",
                i,
                check[i],
                rhs[i]
            );
        }
    }
}
