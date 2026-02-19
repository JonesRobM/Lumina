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

use ndarray::Array2;
use num_complex::Complex64;

/// Compute the 3x3 dyadic Green's tensor between two points.
///
/// # Arguments
/// * `r1` - Position of the observation point (nm).
/// * `r2` - Position of the source point (nm).
/// * `k` - Wavenumber in the medium (nm^{-1}), i.e. $k = 2\pi n / \lambda$.
///
/// # Returns
/// A 3x3 complex matrix representing $\mathbf{G}(\mathbf{r}_1, \mathbf{r}_2)$.
///
/// # Panics
/// Panics if `r1 == r2` (self-interaction is not defined via this function).
pub fn dyadic_greens_tensor(r1: &[f64; 3], r2: &[f64; 3], k: f64) -> Array2<Complex64> {
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
    let mut g = Array2::<Complex64>::zeros((3, 3));
    for i in 0..3 {
        for j in 0..3 {
            let delta_ij = if i == j { 1.0 } else { 0.0 };
            g[[i, j]] = prefactor * (a * delta_ij + b * r_hat[i] * r_hat[j]);
        }
    }

    g
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
                assert_abs_diff_eq!(g_12[[i, j]].re, g_21[[j, i]].re, epsilon = 1e-12);
                assert_abs_diff_eq!(g_12[[i, j]].im, g_21[[j, i]].im, epsilon = 1e-12);
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
}
