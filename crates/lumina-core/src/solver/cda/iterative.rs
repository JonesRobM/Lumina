//! Iterative GMRES solver for large CDA systems.
//!
//! For systems with $N > 1000$ dipoles, direct LU decomposition becomes
//! prohibitively expensive ($O(N^3)$). GMRES (Generalised Minimal Residual
//! method) solves the system iteratively using only matrix-vector products,
//! which can be accelerated via GPU or parallelised across nodes.
//!
//! This implements GMRES(m) with restart, using Arnoldi iteration and
//! Givens rotations for the least-squares solve.

use ndarray::{Array1, Array2};
use num_complex::Complex64;

use super::super::SolverError;

/// Solve the interaction system using restarted GMRES(m).
///
/// # Arguments
/// * `matrix` - The $3N \times 3N$ interaction matrix.
/// * `rhs` - The incident field vector (length $3N$).
/// * `tolerance` - Convergence tolerance on the relative residual.
/// * `max_iterations` - Maximum total number of iterations across all restarts.
///
/// # Returns
/// The solved dipole moment vector $\mathbf{p}$ (length $3N$).
pub fn solve_gmres(
    matrix: &Array2<Complex64>,
    rhs: &Array1<Complex64>,
    tolerance: f64,
    max_iterations: usize,
) -> Result<Array1<Complex64>, SolverError> {
    let n = rhs.len();
    let restart_dim = 30.min(n); // Krylov subspace dimension per restart cycle

    // Initial guess: zero vector
    let mut x = Array1::<Complex64>::zeros(n);

    let rhs_norm = vector_norm(rhs);
    if rhs_norm < 1e-30 {
        return Ok(x);
    }

    let abs_tol = tolerance * rhs_norm;
    let mut total_iters = 0;

    loop {
        // Compute residual: r = rhs - A*x
        let r = rhs - &matvec(matrix, &x);
        let beta = vector_norm(&r);

        if beta < abs_tol {
            return Ok(x);
        }

        if total_iters >= max_iterations {
            return Err(SolverError::ConvergenceFailure {
                max_iter: max_iterations,
                residual: beta / rhs_norm,
            });
        }

        // Arnoldi iteration with Givens rotations
        // V: orthonormal basis vectors (columns of V are v_0..v_m)
        let mut v_basis: Vec<Array1<Complex64>> = Vec::with_capacity(restart_dim + 1);
        v_basis.push(r.mapv(|c| c / Complex64::from(beta)));

        // Upper Hessenberg matrix H (stored as (m+1) x m)
        let mut h = Array2::<Complex64>::zeros((restart_dim + 1, restart_dim));

        // Givens rotation parameters
        let mut cs: Vec<f64> = Vec::with_capacity(restart_dim);
        let mut sn: Vec<Complex64> = Vec::with_capacity(restart_dim);

        // RHS of the least-squares problem: g = beta * e_1
        let mut g = Array1::<Complex64>::zeros(restart_dim + 1);
        g[0] = Complex64::from(beta);

        let mut j = 0;

        while j < restart_dim && total_iters < max_iterations {
            // Arnoldi step: w = A * v_j
            let w = matvec(matrix, &v_basis[j]);
            let mut w = w;

            // Modified Gram-Schmidt orthogonalisation
            for i in 0..=j {
                let h_ij = vector_dot(&v_basis[i], &w);
                h[[i, j]] = h_ij;
                w = w - &v_basis[i].mapv(|c| c * h_ij);
            }

            let h_jp1_j = vector_norm(&w);
            h[[j + 1, j]] = Complex64::from(h_jp1_j);

            if h_jp1_j > 1e-30 {
                v_basis.push(w.mapv(|c| c / Complex64::from(h_jp1_j)));
            } else {
                // Lucky breakdown — exact solution in the Krylov subspace
                v_basis.push(Array1::zeros(n));
            }

            // Apply previous Givens rotations to the new column of H
            for i in 0..j {
                let temp = cs[i] * h[[i, j]] + sn[i] * h[[i + 1, j]];
                h[[i + 1, j]] = -sn[i].conj() * h[[i, j]] + Complex64::from(cs[i]) * h[[i + 1, j]];
                h[[i, j]] = temp;
            }

            // Compute new Givens rotation to eliminate h[j+1, j]
            let (c, s) = givens_rotation(h[[j, j]], h[[j + 1, j]]);
            cs.push(c);
            sn.push(s);

            // Apply new rotation to H and g
            let temp = Complex64::from(c) * h[[j, j]] + s * h[[j + 1, j]];
            h[[j + 1, j]] = Complex64::from(0.0);
            h[[j, j]] = temp;

            let temp = Complex64::from(c) * g[j] + s * g[j + 1];
            g[j + 1] = -s.conj() * g[j] + Complex64::from(c) * g[j + 1];
            g[j] = temp;

            total_iters += 1;
            j += 1;

            // Check convergence
            if g[j].norm() < abs_tol {
                break;
            }
        }

        // Solve the upper triangular system H[0..j, 0..j] * y = g[0..j]
        let y = back_substitute(&h, &g, j);

        // Update solution: x = x + V * y
        for i in 0..j {
            x = x + &v_basis[i].mapv(|c| c * y[i]);
        }

        // Check if we've converged
        let residual = rhs - &matvec(matrix, &x);
        let res_norm = vector_norm(&residual);
        if res_norm < abs_tol {
            return Ok(x);
        }

        if total_iters >= max_iterations {
            return Err(SolverError::ConvergenceFailure {
                max_iter: max_iterations,
                residual: res_norm / rhs_norm,
            });
        }
    }
}

/// Matrix-vector product: y = A * x.
fn matvec(a: &Array2<Complex64>, x: &Array1<Complex64>) -> Array1<Complex64> {
    a.dot(x)
}

/// Complex inner product: <a, b> = sum(conj(a_i) * b_i).
fn vector_dot(a: &Array1<Complex64>, b: &Array1<Complex64>) -> Complex64 {
    a.iter()
        .zip(b.iter())
        .map(|(ai, bi)| ai.conj() * bi)
        .sum()
}

/// Euclidean norm of a complex vector.
fn vector_norm(v: &Array1<Complex64>) -> f64 {
    v.iter().map(|c| c.norm_sqr()).sum::<f64>().sqrt()
}

/// Compute Givens rotation parameters to zero out b in [a; b].
///
/// Returns (c, s) such that:
///   c * a + s * b = r (real, non-negative)
///   -conj(s) * a + c * b = 0
fn givens_rotation(a: Complex64, b: Complex64) -> (f64, Complex64) {
    if b.norm() < 1e-30 {
        return (1.0, Complex64::from(0.0));
    }
    if a.norm() < 1e-30 {
        return (0.0, b / Complex64::from(b.norm()));
    }
    let r = (a.norm_sqr() + b.norm_sqr()).sqrt();
    let c = a.norm() / r;
    let s = (a / Complex64::from(a.norm())) * b.conj() / Complex64::from(r);
    (c, s)
}

/// Back-substitution for an upper-triangular system H[0..m, 0..m] * y = g[0..m].
fn back_substitute(h: &Array2<Complex64>, g: &Array1<Complex64>, m: usize) -> Vec<Complex64> {
    let mut y = vec![Complex64::from(0.0); m];
    for i in (0..m).rev() {
        let mut sum = g[i];
        for j in (i + 1)..m {
            sum -= h[[i, j]] * y[j];
        }
        if h[[i, i]].norm() > 1e-30 {
            y[i] = sum / h[[i, i]];
        }
    }
    y
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    #[test]
    fn test_gmres_identity_system() {
        // Solve Ix = b, should give x = b exactly
        let n = 6;
        let mut matrix = Array2::<Complex64>::zeros((n, n));
        for i in 0..n {
            matrix[[i, i]] = Complex64::from(1.0);
        }
        let rhs = Array1::from_vec(
            (0..n)
                .map(|i| Complex64::new(i as f64 + 1.0, 0.5 * i as f64))
                .collect(),
        );

        let x = solve_gmres(&matrix, &rhs, 1e-10, 100).unwrap();

        for i in 0..n {
            assert!(
                (x[i] - rhs[i]).norm() < 1e-8,
                "x[{}] = {:?}, expected {:?}",
                i,
                x[i],
                rhs[i]
            );
        }
    }

    #[test]
    fn test_gmres_complex_system() {
        // Solve a 3x3 complex system and verify residual
        let matrix = array![
            [
                Complex64::new(4.0, 1.0),
                Complex64::new(-1.0, 0.5),
                Complex64::new(0.0, 0.0)
            ],
            [
                Complex64::new(-1.0, -0.5),
                Complex64::new(4.0, 0.0),
                Complex64::new(-1.0, 0.2)
            ],
            [
                Complex64::new(0.0, 0.0),
                Complex64::new(-1.0, -0.2),
                Complex64::new(4.0, -1.0)
            ],
        ];
        let rhs = array![
            Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 1.0),
            Complex64::new(-1.0, 0.5),
        ];

        let x = solve_gmres(&matrix, &rhs, 1e-12, 100).unwrap();

        // Verify: ||A*x - b|| / ||b|| < tolerance
        let residual = &rhs - &matrix.dot(&x);
        let res_norm = vector_norm(&residual);
        let rhs_norm = vector_norm(&rhs);
        assert!(
            res_norm / rhs_norm < 1e-10,
            "Relative residual {:.2e} too large",
            res_norm / rhs_norm
        );
    }

    #[test]
    fn test_gmres_agrees_with_known_solution() {
        // Diagonal system: trivial exact solution
        let n = 9;
        let mut matrix = Array2::<Complex64>::zeros((n, n));
        let mut rhs = Array1::<Complex64>::zeros(n);
        for i in 0..n {
            let alpha = Complex64::new(2.0 + i as f64, 0.3 * i as f64);
            matrix[[i, i]] = alpha;
            rhs[i] = alpha * Complex64::new(1.0, -0.5); // x_i = 1 - 0.5i for all i
        }

        let x = solve_gmres(&matrix, &rhs, 1e-12, 100).unwrap();

        for i in 0..n {
            let expected = Complex64::new(1.0, -0.5);
            assert!(
                (x[i] - expected).norm() < 1e-8,
                "x[{}] = {:?}, expected {:?}",
                i,
                x[i],
                expected
            );
        }
    }

    #[test]
    fn test_gmres_convergence_failure() {
        // Zero matrix should not converge (if rhs is non-zero)
        let matrix = Array2::<Complex64>::zeros((3, 3));
        let rhs = array![
            Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 1.0),
            Complex64::new(1.0, 1.0),
        ];

        let result = solve_gmres(&matrix, &rhs, 1e-10, 50);
        assert!(result.is_err(), "Should fail to converge for singular system");
    }
}
