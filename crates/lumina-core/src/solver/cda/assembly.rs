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

use ndarray::Array2;
use num_complex::Complex64;

use crate::types::Dipole;

/// Assemble the full $3N \times 3N$ interaction matrix.
///
/// # Arguments
/// * `dipoles` - Slice of dipoles with positions and polarisabilities.
/// * `k` - Wavenumber in the medium (nm^{-1}).
///
/// # Returns
/// The interaction matrix $\mathbf{A}$ such that $\mathbf{A}\mathbf{p} = \mathbf{E}_{\text{inc}}$.
pub fn assemble_interaction_matrix(dipoles: &[Dipole], k: f64) -> Array2<Complex64> {
    let n = dipoles.len();
    let dim = 3 * n;
    let mut matrix = Array2::<Complex64>::zeros((dim, dim));

    for i in 0..n {
        // Diagonal block: inverse polarisability
        // The polarisability is stored as a flat 3x3 in row-major order
        for row in 0..3 {
            for col in 0..3 {
                matrix[[3 * i + row, 3 * i + col]] = dipoles[i].polarisability[3 * row + col];
            }
        }
        // TODO: invert the 3x3 polarisability block on the diagonal

        // Off-diagonal blocks: -G(r_i, r_j)
        for j in 0..n {
            if i == j {
                continue;
            }
            let g = super::greens::dyadic_greens_tensor(&dipoles[i].position, &dipoles[j].position, k);
            for row in 0..3 {
                for col in 0..3 {
                    matrix[[3 * i + row, 3 * j + col]] = -g[[row, col]];
                }
            }
        }
    }

    matrix
}
