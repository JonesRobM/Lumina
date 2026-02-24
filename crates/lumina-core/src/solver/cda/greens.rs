//! Dyadic Green's function for free-space electromagnetic interaction.
//!
//! The interaction between dipoles $i$ and $j$ is governed by the free-space
//! dyadic Green's function:
//!
//! $$
//! \mathbf{G}(\mathbf{r}, \mathbf{r}') = \frac{k^2 e^{ikR}}{4\pi \epsilon_0 R}
//! \left[ \left(1 + \frac{ikR - 1}{k^2 R^2}\right)\mathbf{I}
//! + \frac{3 - 3ikR - k^2 R^2}{k^2 R^2} \frac{\mathbf{R}\mathbf{R}^T}{R^2} \right]
//! $$
//!
//! where $\mathbf{R} = \mathbf{r} - \mathbf{r}'$ and $R = |\mathbf{R}|$.

use num_complex::Complex64;

/// Stack-allocated 3×3 complex tensor (zero heap allocation).
pub type Tensor3x3 = [[Complex64; 3]; 3];

/// Compute the 3x3 dyadic Green's tensor between two points.
///
/// Returns a stack-allocated `[[Complex64; 3]; 3]` to avoid heap allocation
/// in the inner assembly loop (~N² calls per wavelength).
///
/// # Arguments
/// * `r1` - Position of the observation point (nm).
/// * `r2` - Position of the source point (nm).
/// * `k` - Wavenumber in the medium (nm^{-1}), i.e. $k = 2\pi n / \lambda$.
///
/// # Returns
/// A 3x3 complex tensor representing $\mathbf{G}(\mathbf{r}_1, \mathbf{r}_2)$.
///
/// # Panics
/// Panics if `r1 == r2` (self-interaction is not defined via this function).
pub fn dyadic_greens_tensor(r1: &[f64; 3], r2: &[f64; 3], k: f64) -> Tensor3x3 {
    let rx = r1[0] - r2[0];
    let ry = r1[1] - r2[1];
    let rz = r1[2] - r2[2];
    let r_sq = rx * rx + ry * ry + rz * rz;
    let r = r_sq.sqrt();

    assert!(r > 1e-15, "Self-interaction: r1 and r2 must not coincide");

    let kr = k * r;
    let kr_sq = kr * kr;
    let ikr = Complex64::new(0.0, kr);
    let exp_ikr = ikr.exp();

    // Prefactor: k^2 * exp(ikR) / (4 * pi * R)
    // Note: we absorb epsilon_0 into the polarisability definition.
    let prefactor = k * k * exp_ikr / (4.0 * std::f64::consts::PI * r);

    // Scalar coefficients for the identity and dyadic terms
    let a = Complex64::from(1.0) + (ikr - Complex64::from(1.0)) / Complex64::from(kr_sq);
    let b = (Complex64::from(3.0) - Complex64::from(3.0) * ikr - Complex64::from(kr_sq))
        / Complex64::from(kr_sq);

    // Unit displacement vector components
    let r_hat = [rx / r, ry / r, rz / r];

    // Build the 3x3 tensor: G_ij = prefactor * (a * delta_ij + b * r_hat_i * r_hat_j)
    let zero = Complex64::from(0.0);
    let mut g = [[zero; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let delta_ij = if i == j { 1.0 } else { 0.0 };
            g[i][j] = prefactor * (a * delta_ij + b * r_hat[i] * r_hat[j]);
        }
    }

    g
}

/// Compute the volume-averaged (filtered) Green's tensor using Gauss-Legendre
/// quadrature over the source cell.
///
/// For near-field interactions (|r_obs - r_src| <= 2*cell_size), the point-dipole
/// Green's function has a 1/R³ singularity that creates staircase artefacts at
/// material boundaries. The Integrated Green's Tensor (IGT) method (Yurkin &
/// Hoekstra, JQSRT 2007) replaces G(r_obs, r_src) with the volume average:
///
/// $G_{\text{IGT}}(\mathbf{r}_i, \mathbf{r}_j) = \frac{1}{V_j} \int_{V_j} G(\mathbf{r}_i, \mathbf{r}') \, d\mathbf{r}'$
///
/// This smooths the near-field singularity and dramatically improves accuracy
/// for metallic particles.
///
/// For distant pairs (|r_obs - r_src| > 2*cell_size), the standard point-dipole
/// Green's function is returned (the volume-averaging correction is negligible).
///
/// # Arguments
/// * `r_obs` - Position of the observation point (nm).
/// * `r_src` - Centre of the source cell (nm).
/// * `k` - Wavenumber in the medium (nm⁻¹).
/// * `cell_size` - Side length of the cubic dipole cell (nm).
pub fn dyadic_greens_tensor_filtered(
    r_obs: &[f64; 3],
    r_src: &[f64; 3],
    k: f64,
    cell_size: f64,
) -> Tensor3x3 {
    let dx = r_obs[0] - r_src[0];
    let dy = r_obs[1] - r_src[1];
    let dz = r_obs[2] - r_src[2];
    let dist = (dx * dx + dy * dy + dz * dz).sqrt();

    // For distant pairs, the point-dipole approximation is sufficient
    if dist > 2.0 * cell_size {
        return dyadic_greens_tensor(r_obs, r_src, k);
    }

    // 3-point Gauss-Legendre quadrature on [-1, 1]
    let gl_nodes: [f64; 3] = [
        -(3.0_f64 / 5.0).sqrt(),
        0.0,
        (3.0_f64 / 5.0).sqrt(),
    ];
    let gl_weights: [f64; 3] = [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0];

    let half_d = cell_size / 2.0;

    let zero = Complex64::from(0.0);
    let mut g_avg = [[zero; 3]; 3];
    let mut total_weight = 0.0;

    // 3D quadrature: 3^3 = 27 evaluation points per cell
    for (ix, &nx) in gl_nodes.iter().enumerate() {
        for (iy, &ny) in gl_nodes.iter().enumerate() {
            for (iz, &nz) in gl_nodes.iter().enumerate() {
                let r_prime = [
                    r_src[0] + half_d * nx,
                    r_src[1] + half_d * ny,
                    r_src[2] + half_d * nz,
                ];

                // Skip if quadrature point coincides with observation point
                let ddx = r_obs[0] - r_prime[0];
                let ddy = r_obs[1] - r_prime[1];
                let ddz = r_obs[2] - r_prime[2];
                let r_sq = ddx * ddx + ddy * ddy + ddz * ddz;
                if r_sq < 1e-20 {
                    continue;
                }

                let w = gl_weights[ix] * gl_weights[iy] * gl_weights[iz];
                let g = dyadic_greens_tensor(r_obs, &r_prime, k);

                for a in 0..3 {
                    for b in 0..3 {
                        g_avg[a][b] += Complex64::from(w) * g[a][b];
                    }
                }
                total_weight += w;
            }
        }
    }

    // Normalise by total weight (should be 8.0 for 3-point GL on [-1,1]^3)
    if total_weight > 1e-30 {
        let inv_w = Complex64::from(1.0 / total_weight);
        for row in &mut g_avg {
            for val in row.iter_mut() {
                *val *= inv_w;
            }
        }
    }

    g_avg
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_greens_reciprocity() {
        // G(r, r') should equal G^T(r', r)
        let r1 = [0.0, 0.0, 0.0];
        let r2 = [10.0, 5.0, 3.0];
        let k = 2.0 * std::f64::consts::PI / 500.0; // lambda = 500 nm

        let g_12 = dyadic_greens_tensor(&r1, &r2, k);
        let g_21 = dyadic_greens_tensor(&r2, &r1, k);

        for i in 0..3 {
            for j in 0..3 {
                assert_abs_diff_eq!(g_12[i][j].re, g_21[j][i].re, epsilon = 1e-12);
                assert_abs_diff_eq!(g_12[i][j].im, g_21[j][i].im, epsilon = 1e-12);
            }
        }
    }

    #[test]
    #[should_panic(expected = "Self-interaction")]
    fn test_greens_self_interaction_panics() {
        let r = [1.0, 2.0, 3.0];
        let k = 0.01;
        dyadic_greens_tensor(&r, &r, k);
    }

    #[test]
    fn test_fcd_converges_to_point_dipole_at_distance() {
        // For large separation (R >> d), FCD should agree with the point-dipole G
        let r1 = [0.0, 0.0, 0.0];
        let r2 = [30.0, 0.0, 0.0]; // 10 cell widths away
        let k = 2.0 * std::f64::consts::PI / 500.0;
        let cell_size = 3.0;

        let g_point = dyadic_greens_tensor(&r1, &r2, k);
        let g_fcd = dyadic_greens_tensor_filtered(&r1, &r2, k, cell_size);

        for i in 0..3 {
            for j in 0..3 {
                assert_abs_diff_eq!(g_point[i][j].re, g_fcd[i][j].re, epsilon = 1e-10);
                assert_abs_diff_eq!(g_point[i][j].im, g_fcd[i][j].im, epsilon = 1e-10);
            }
        }
    }

    #[test]
    fn test_fcd_is_finite_for_nearest_neighbours() {
        // For R = d (nearest neighbour), FCD should be smooth and finite
        let r1 = [0.0, 0.0, 0.0];
        let r2 = [3.0, 0.0, 0.0]; // one cell width apart
        let k = 2.0 * std::f64::consts::PI / 500.0;
        let cell_size = 3.0;

        let g_fcd = dyadic_greens_tensor_filtered(&r1, &r2, k, cell_size);

        for i in 0..3 {
            for j in 0..3 {
                assert!(g_fcd[i][j].re.is_finite(), "FCD G[{},{}] real part not finite", i, j);
                assert!(g_fcd[i][j].im.is_finite(), "FCD G[{},{}] imag part not finite", i, j);
            }
        }
    }

    #[test]
    fn test_fcd_differs_from_point_at_short_range() {
        // At R = d, FCD should differ meaningfully from point-dipole G
        let r1 = [0.0, 0.0, 0.0];
        let r2 = [3.0, 0.0, 0.0];
        let k = 2.0 * std::f64::consts::PI / 500.0;
        let cell_size = 3.0;

        let g_point = dyadic_greens_tensor(&r1, &r2, k);
        let g_fcd = dyadic_greens_tensor_filtered(&r1, &r2, k, cell_size);

        // At least one element should show a meaningful difference
        let mut max_diff = 0.0_f64;
        for i in 0..3 {
            for j in 0..3 {
                let diff = (g_fcd[i][j] - g_point[i][j]).norm();
                max_diff = max_diff.max(diff);
            }
        }
        assert!(
            max_diff > 1e-8,
            "FCD should differ from point-dipole G at short range, max diff = {:.2e}",
            max_diff
        );
    }
}
