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
use crate::solver::ewald::EwaldGreens;
use crate::solver::substrate::SubstrateRuntime;
use lumina_geometry::lattice::LatticeSpec;

/// Assemble the full $3N \times 3N$ interaction matrix.
///
/// Off-diagonal blocks (Green's tensor evaluations) are computed in parallel
/// using Rayon. Diagonal blocks (inverse polarisability) are filled sequentially
/// since they are cheap.
///
/// # Arguments
/// * `dipoles`   - Slice of dipoles with positions and polarisabilities.
/// * `k`         - Wavenumber in the medium (nm⁻¹).
/// * `use_fcd`   - If true, use the FCD volume-averaged Green's function for
///   near-field interactions. Ignored for periodic assemblies.
/// * `cell_size` - Dipole lattice spacing (nm). Used by FCD. Ignored if
///   `use_fcd` is false.
/// * `lattice`   - If `Some`, uses the periodic Ewald Green's function instead
///   of the free-space Green's function.
/// * `k_bloch`   - In-plane Bloch wavevector **k**∥ (nm⁻¹). Used only when
///   `lattice` is `Some`. Pass `[0;3]` for the Γ point.
/// * `substrate` - If `Some`, adds the substrate correction (image-dipole or
///   Sommerfeld integral) for each dipole pair. When both `lattice` and
///   `substrate` are present, the Fresnel image-dipole sum is also
///   Ewald-accelerated; the Sommerfeld path always uses free-space G.
///
/// # Returns
/// The interaction matrix $\mathbf{A}$ such that $\mathbf{A}\mathbf{p} = \mathbf{E}_{\text{inc}}$.
pub fn assemble_interaction_matrix(
    dipoles: &[Dipole],
    k: f64,
    use_fcd: bool,
    cell_size: f64,
    lattice: Option<&LatticeSpec>,
    k_bloch: [f64; 3],
    substrate: Option<&SubstrateRuntime>,
) -> Array2<Complex64> {
    let n = dipoles.len();
    let dim = 3 * n;
    let mut matrix = Array2::<Complex64>::zeros((dim, dim));

    // Diagonal blocks: inverse polarisability (cheap, sequential).
    // For periodic arrays, the Ewald self-image sum (R≠0) is added to the
    // diagonal in the off-diagonal loop below (i==j, R≠0 lattice images).
    for i in 0..n {
        let inv_alpha = invert_3x3(&dipoles[i].polarisability);
        for row in 0..3 {
            for col in 0..3 {
                matrix[[3 * i + row, 3 * i + col]] = inv_alpha[3 * row + col];
            }
        }
    }

    // Off-diagonal blocks: −G(r_i, r_j), computed in parallel.
    // For periodic systems, each off-diagonal block also includes the sum
    // over lattice images via EwaldGreens.
    //
    // SAFETY: each Rayon thread handles a unique row-block i, so writes to
    // rows [3i..3i+2] are data-race-free.
    let matrix_ptr = matrix.as_mut_ptr() as usize;
    let stride = dim;

    // Wrap the optional lattice in an Arc so the Rayon closure can share it.
    let ewald: Option<std::sync::Arc<EwaldGreens>> = lattice
        .map(|lat| std::sync::Arc::new(EwaldGreens::new(lat.clone())));

    (0..n).into_par_iter().for_each(|i| {
        // Clone the Arc for this thread (cheap — increments ref-count only).
        let ewald_ref = ewald.as_ref().map(std::sync::Arc::as_ref);

        for j in 0..n {
            // ── Direct G contribution ──────────────────────────────────────
            // In free-space (no Ewald), the self-block i==j is skipped:
            // the α⁻¹ diagonal is already set above and there is no
            // R=0 lattice-image contribution.
            let skip_direct = i == j && ewald_ref.is_none();
            if !skip_direct {
                let g = if let Some(eg) = ewald_ref {
                    // Periodic: Ewald sum over all lattice images of j.
                    // When i==j, the R=0 self-interaction is excluded inside
                    // EwaldGreens (only non-zero image lattice contributions).
                    let r = [
                        dipoles[i].position[0] - dipoles[j].position[0],
                        dipoles[i].position[1] - dipoles[j].position[1],
                        dipoles[i].position[2] - dipoles[j].position[2],
                    ];
                    eg.evaluate(r, k_bloch, k)
                } else if use_fcd {
                    super::greens::dyadic_greens_tensor_filtered(
                        &dipoles[i].position,
                        &dipoles[j].position,
                        k,
                        cell_size,
                    )
                } else {
                    super::greens::dyadic_greens_tensor(
                        &dipoles[i].position,
                        &dipoles[j].position,
                        k,
                    )
                };

                for (row, g_row) in g.iter().enumerate() {
                    for (col, &g_val) in g_row.iter().enumerate() {
                        unsafe {
                            let offset = (3 * i + row) * stride + (3 * j + col);
                            *(matrix_ptr as *mut Complex64).add(offset) -= g_val;
                        }
                    }
                }
            }

            // ── Substrate contribution ──────────────────────────────────────
            // Runs for ALL (i, j) pairs, including the self-image i == j.
            if let Some(sub_rt) = substrate {
                let rj = dipoles[j].position;
                let ri = dipoles[i].position;
                match sub_rt {
                    SubstrateRuntime::Fresnel { z_interface, delta_eps } => {
                        // Image position: mirror of r_j through the interface.
                        let rj_img = [rj[0], rj[1], 2.0 * z_interface - rj[2]];
                        let diff = [ri[0] - rj_img[0], ri[1] - rj_img[1], ri[2] - rj_img[2]];
                        let dist_sq = diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2];
                        if dist_sq >= 1e-20 {
                            let g_img = if let Some(eg) = ewald_ref {
                                eg.evaluate(diff, k_bloch, k)
                            } else {
                                super::greens::dyadic_greens_tensor(&ri, &rj_img, k)
                            };
                            let m = [-delta_eps, -delta_eps, *delta_eps];
                            for (row, g_row) in g_img.iter().enumerate() {
                                for (col, &g_val) in g_row.iter().enumerate() {
                                    unsafe {
                                        let offset = (3 * i + row) * stride + (3 * j + col);
                                        *(matrix_ptr as *mut Complex64).add(offset) -= m[col] * g_val;
                                    }
                                }
                            }
                        }
                    }
                    SubstrateRuntime::Sommerfeld { z_interface, greens } => {
                        // Full Sommerfeld reflected Green's tensor — no image position needed.
                        let g_refl = greens.evaluate(ri, rj, *z_interface);
                        for (row, g_row) in g_refl.iter().enumerate() {
                            for (col, &g_val) in g_row.iter().enumerate() {
                                unsafe {
                                    let offset = (3 * i + row) * stride + (3 * j + col);
                                    *(matrix_ptr as *mut Complex64).add(offset) -= g_val;
                                }
                            }
                        }
                    }
                }
            }
        }
    });

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

/// Compute the matrix-vector product **A · x** without explicitly assembling
/// the $3N \times 3N$ interaction matrix.
///
/// Functionally identical to `assemble_interaction_matrix(…).dot(x)` but
/// requires only $O(N)$ working memory instead of $O(N^2)$, making it safe
/// for large systems (N > 3 000) where the full matrix would exhaust RAM.
///
/// Each output row-block $i$ is independent, so the outer loop is
/// parallelised with Rayon. The cost per call is $O(N^2)$ — the same as one
/// explicit matrix–vector product.
///
/// # Arguments
/// Same as [`assemble_interaction_matrix`], plus `x` (the input vector, length $3N$).
pub fn matvec_on_the_fly(
    dipoles: &[Dipole],
    x: &Array1<Complex64>,
    k: f64,
    use_fcd: bool,
    cell_size: f64,
    lattice: Option<&LatticeSpec>,
    k_bloch: [f64; 3],
    substrate: Option<&SubstrateRuntime>,
) -> Array1<Complex64> {
    let n = dipoles.len();
    let ewald: Option<std::sync::Arc<EwaldGreens>> = lattice
        .map(|lat| std::sync::Arc::new(EwaldGreens::new(lat.clone())));

    // Compute each row-block i independently in parallel.
    let rows: Vec<[Complex64; 3]> = (0..n)
        .into_par_iter()
        .map(|i| {
            let ewald_ref = ewald.as_ref().map(std::sync::Arc::as_ref);
            let mut yi = [Complex64::from(0.0); 3];

            // Diagonal block: α_i^{-1} · x[3i..3i+2]
            let inv_alpha = invert_3x3(&dipoles[i].polarisability);
            for a in 0..3 {
                for b in 0..3 {
                    yi[a] += inv_alpha[3 * a + b] * x[3 * i + b];
                }
            }

            for j in 0..n {
                // ── Direct Green's function contribution ────────────────────
                // For free-space (no Ewald) the i==j self-block is skipped;
                // the α⁻¹ diagonal is already handled above.
                let skip_direct = i == j && ewald_ref.is_none();
                if !skip_direct {
                    let g = if let Some(eg) = ewald_ref {
                        let r = [
                            dipoles[i].position[0] - dipoles[j].position[0],
                            dipoles[i].position[1] - dipoles[j].position[1],
                            dipoles[i].position[2] - dipoles[j].position[2],
                        ];
                        eg.evaluate(r, k_bloch, k)
                    } else if use_fcd {
                        super::greens::dyadic_greens_tensor_filtered(
                            &dipoles[i].position,
                            &dipoles[j].position,
                            k,
                            cell_size,
                        )
                    } else {
                        super::greens::dyadic_greens_tensor(
                            &dipoles[i].position,
                            &dipoles[j].position,
                            k,
                        )
                    };
                    for a in 0..3 {
                        for b in 0..3 {
                            yi[a] -= g[a][b] * x[3 * j + b];
                        }
                    }
                }

                // ── Substrate contribution ───────────────────────────────────
                // Runs for ALL j (including i == j self-image).
                if let Some(sub_rt) = substrate {
                    let rj = dipoles[j].position;
                    let ri = dipoles[i].position;
                    match sub_rt {
                        SubstrateRuntime::Fresnel { z_interface, delta_eps } => {
                            let rj_img = [rj[0], rj[1], 2.0 * z_interface - rj[2]];
                            let diff = [ri[0]-rj_img[0], ri[1]-rj_img[1], ri[2]-rj_img[2]];
                            let dist_sq = diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2];
                            if dist_sq >= 1e-20 {
                                let g_img = if let Some(eg) = ewald_ref {
                                    eg.evaluate(diff, k_bloch, k)
                                } else {
                                    super::greens::dyadic_greens_tensor(&ri, &rj_img, k)
                                };
                                let m = [-delta_eps, -delta_eps, *delta_eps];
                                for a in 0..3 {
                                    for b in 0..3 {
                                        yi[a] -= m[b] * g_img[a][b] * x[3 * j + b];
                                    }
                                }
                            }
                        }
                        SubstrateRuntime::Sommerfeld { z_interface, greens } => {
                            let g_refl = greens.evaluate(ri, rj, *z_interface);
                            for a in 0..3 {
                                for b in 0..3 {
                                    yi[a] -= g_refl[a][b] * x[3 * j + b];
                                }
                            }
                        }
                    }
                }
            }
            yi
        })
        .collect();

    let mut y = Array1::<Complex64>::zeros(3 * n);
    for (i, yi) in rows.iter().enumerate() {
        for c in 0..3 {
            y[3 * i + c] = yi[c];
        }
    }
    y
}

/// Compute the 3×3 inverse polarisability block for each dipole.
///
/// Returns a `Vec` of length N, each a column-major `[Complex64; 9]`.
/// This is a pure inversion of each dipole's pre-built polarisability tensor;
/// no k, FCD, or cell_size dependency.
pub fn compute_diagonal_blocks(dipoles: &[Dipole]) -> Vec<[Complex64; 9]> {
    dipoles.iter().map(|d| invert_3x3(&d.polarisability)).collect()
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
