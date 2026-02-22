//! Interaction matrix assembly for the CDA.
//!
//! Constructs the $3N \times 3N$ complex interaction matrix $\mathbf{A}$
//! where each $3 \times 3$ block $(i, j)$ is:
//!
//! - Diagonal ($i = j$): $\alpha_i^{-1}$ (inverse polarisability)
//! - Off-diagonal ($i \neq j$): $-\mathbf{G}(\mathbf{r}_i, \mathbf{r}_j)$
//!
//! The self-consistent dipole moments $\mathbf{p}$ satisfy:
//! $$\mathbf{A} \mathbf{p} = \mathbf{E}_{\text{inc}}$$

use ndarray::{Array1, Array2};
use num_complex::Complex64;
use rayon::prelude::*;

use crate::types::{Dipole, IncidentField};

/// Assemble the full $3N \times 3N$ interaction matrix.
///
/// Off-diagonal blocks (Green's tensor evaluations) are computed in parallel
/// using Rayon. Diagonal blocks (inverse polarisability) are filled sequentially
/// since they are cheap.
///
/// # Arguments
/// * `dipoles` - Slice of dipoles with positions and polarisabilities.
/// * `k` - Wavenumber in the medium (nm^{-1}).
/// * `use_fcd` - If true, use the Filtered Coupled Dipole (volume-averaged)
///   Green's function for near-field interactions. Requires `cell_size`.
/// * `cell_size` - Dipole lattice spacing (nm). Used by FCD for the
///   integration volume. Ignored if `use_fcd` is false.
///
/// # Returns
/// The interaction matrix $\mathbf{A}$ such that $\mathbf{A}\mathbf{p} = \mathbf{E}_{\text{inc}}$.
pub fn assemble_interaction_matrix(
    dipoles: &[Dipole],
    k: f64,
    use_fcd: bool,
    cell_size: f64,
) -> Array2<Complex64> {
    let n = dipoles.len();
    let dim = 3 * n;
    let mut matrix = Array2::<Complex64>::zeros((dim, dim));

    // Diagonal blocks: inverse polarisability (cheap, sequential)
    for i in 0..n {
        let inv_alpha = invert_3x3(&dipoles[i].polarisability);
        for row in 0..3 {
            for col in 0..3 {
                matrix[[3 * i + row, 3 * i + col]] = inv_alpha[3 * row + col];
            }
        }
    }

    // Off-diagonal blocks: -G(r_i, r_j), computed in parallel.
    // Each block is a 3×3 Green's tensor — compute all of them via Rayon,
    // then place into the matrix sequentially.
    let blocks: Vec<(usize, usize, [[Complex64; 3]; 3])> = (0..n)
        .into_par_iter()
        .flat_map(|i| {
            let dipoles_ref = dipoles;
            (0..n).into_par_iter().filter(move |&j| i != j).map(move |j| {
                let g = if use_fcd {
                    super::greens::dyadic_greens_tensor_filtered(
                        &dipoles_ref[i].position,
                        &dipoles_ref[j].position,
                        k,
                        cell_size,
                    )
                } else {
                    super::greens::dyadic_greens_tensor(
                        &dipoles_ref[i].position,
                        &dipoles_ref[j].position,
                        k,
                    )
                };
                let mut block = [[Complex64::from(0.0); 3]; 3];
                for row in 0..3 {
                    for col in 0..3 {
                        block[row][col] = g[[row, col]];
                    }
                }
                (i, j, block)
            })
        })
        .collect();

    for (i, j, block) in blocks {
        for row in 0..3 {
            for col in 0..3 {
                matrix[[3 * i + row, 3 * j + col]] = -block[row][col];
            }
        }
    }

    matrix
}

/// Construct the incident field vector (length 3N) for all dipoles.
pub fn build_incident_field_vector(
    dipoles: &[Dipole],
    incident: &IncidentField,
    k: f64,
) -> Array1<Complex64> {
    let n = dipoles.len();
    let mut rhs = Array1::<Complex64>::zeros(3 * n);

    for (i, dipole) in dipoles.iter().enumerate() {
        let e = incident.at_position(&dipole.position, k);
        rhs[3 * i] = e[0];
        rhs[3 * i + 1] = e[1];
        rhs[3 * i + 2] = e[2];
    }

    rhs
}

/// Public wrapper for 3x3 inversion, used by cross-section computation.
pub fn invert_3x3_pub(m: &[Complex64; 9]) -> [Complex64; 9] {
    invert_3x3(m)
}

/// Invert a 3x3 complex matrix stored as a flat [9] array (row-major).
fn invert_3x3(m: &[Complex64; 9]) -> [Complex64; 9] {
    // For an isotropic polarisability (diagonal), this simplifies to 1/alpha on the diagonal.
    // General case: compute cofactor matrix and determinant.
    let a = m[0]; let b = m[1]; let c = m[2];
    let d = m[3]; let e = m[4]; let f = m[5];
    let g = m[6]; let h = m[7]; let k = m[8];

    let det = a * (e * k - f * h) - b * (d * k - f * g) + c * (d * h - e * g);

    assert!(
        det.norm() > 1e-30,
        "Singular polarisability tensor (det = {:.2e})",
        det.norm()
    );

    let inv_det = Complex64::from(1.0) / det;

    [
        inv_det * (e * k - f * h),
        inv_det * (c * h - b * k),
        inv_det * (b * f - c * e),
        inv_det * (f * g - d * k),
        inv_det * (a * k - c * g),
        inv_det * (c * d - a * f),
        inv_det * (d * h - e * g),
        inv_det * (b * g - a * h),
        inv_det * (a * e - b * d),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_invert_3x3_identity() {
        let one = Complex64::from(1.0);
        let zero = Complex64::from(0.0);
        let id = [one, zero, zero, zero, one, zero, zero, zero, one];
        let inv = invert_3x3(&id);
        for i in 0..9 {
            let expected = if i % 4 == 0 { 1.0 } else { 0.0 };
            assert!((inv[i].re - expected).abs() < 1e-12);
            assert!(inv[i].im.abs() < 1e-12);
        }
    }

    #[test]
    fn test_invert_3x3_isotropic() {
        let alpha = Complex64::new(2.0, 0.5);
        let zero = Complex64::from(0.0);
        let m = [alpha, zero, zero, zero, alpha, zero, zero, zero, alpha];
        let inv = invert_3x3(&m);
        let expected = Complex64::from(1.0) / alpha;
        assert!((inv[0] - expected).norm() < 1e-12);
        assert!((inv[4] - expected).norm() < 1e-12);
        assert!((inv[8] - expected).norm() < 1e-12);
    }
}
