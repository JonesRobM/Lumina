//! Analytical Mie theory for homogeneous spheres.
//!
//! Provides exact extinction, scattering, and absorption cross-sections
//! for a homogeneous sphere in a uniform medium. Used as the primary
//! validation benchmark for the CDA solver.
//!
//! # Reference
//! Bohren & Huffman, *Absorption and Scattering of Light by Small Particles* (1983).

use num_complex::Complex64;

/// Compute Mie theory cross-sections for a homogeneous sphere.
///
/// # Arguments
/// * `radius_nm` - Sphere radius in nanometres.
/// * `wavelength_nm` - Incident wavelength in nanometres (in the medium).
/// * `epsilon_sphere` - Complex dielectric function of the sphere at this wavelength.
/// * `n_medium` - Refractive index of the surrounding medium (real).
///
/// # Returns
/// Tuple of (extinction, scattering, absorption) cross-sections in nm^2.
pub fn mie_cross_sections(
    radius_nm: f64,
    wavelength_nm: f64,
    epsilon_sphere: Complex64,
    n_medium: f64,
) -> (f64, f64, f64) {
    let k = 2.0 * std::f64::consts::PI * n_medium / wavelength_nm;
    let x = k * radius_nm; // size parameter

    // Relative refractive index
    let n_sphere = epsilon_sphere.sqrt();
    let m = n_sphere / Complex64::from(n_medium);

    // Number of terms (Wiscombe's criterion)
    let n_max = (x + 4.0 * x.powf(1.0 / 3.0) + 2.0).ceil() as usize;
    let n_max = n_max.max(3);

    let (a_coeffs, b_coeffs) = mie_coefficients(x, m, n_max);

    let mut c_ext = 0.0;
    let mut c_sca = 0.0;

    for n in 0..a_coeffs.len() {
        let order = (n + 1) as f64; // n starts from 1 in the physics convention
        let weight = 2.0 * order + 1.0;

        c_ext += weight * (a_coeffs[n].re + b_coeffs[n].re);
        c_sca += weight * (a_coeffs[n].norm_sqr() + b_coeffs[n].norm_sqr());
    }

    let prefactor = 2.0 * std::f64::consts::PI / (k * k);
    c_ext *= prefactor;
    c_sca *= prefactor;
    let c_abs = c_ext - c_sca;

    (c_ext, c_sca, c_abs)
}

/// Compute Mie coefficients $a_n$ and $b_n$.
///
/// Uses the Riccati-Bessel functions:
/// - $\psi_n(z) = z j_n(z)$ (Riccati-Bessel of the first kind)
/// - $\xi_n(z) = z h_n^{(1)}(z)$ (Riccati-Bessel of the third kind)
fn mie_coefficients(
    x: f64,
    m: Complex64,
    n_max: usize,
) -> (Vec<Complex64>, Vec<Complex64>) {
    let mx = m * x;
    // Compute Riccati-Bessel functions via downward recurrence for psi(mx)
    // and upward recurrence for psi(x) and xi(x).
    let psi_mx = riccati_bessel_psi_complex(mx, n_max);
    let dpsi_mx = riccati_bessel_dpsi_complex(mx, n_max, &psi_mx);

    let psi_x = riccati_bessel_psi_real(x, n_max);
    let dpsi_x = riccati_bessel_dpsi_real(x, n_max, &psi_x);

    let xi_x = riccati_bessel_xi_real(x, n_max);
    let dxi_x = riccati_bessel_dxi_real(x, n_max, &xi_x);

    let mut a_coeffs = Vec::with_capacity(n_max);
    let mut b_coeffs = Vec::with_capacity(n_max);

    for n in 1..=n_max {
        // a_n = (m * psi_n(mx) * psi_n'(x) - psi_n(x) * psi_n'(mx))
        //     / (m * psi_n(mx) * xi_n'(x)  - xi_n(x)  * psi_n'(mx))
        let psi_n_mx = psi_mx[n];
        let dpsi_n_mx = dpsi_mx[n];
        let psi_n_x = Complex64::from(psi_x[n]);
        let dpsi_n_x = Complex64::from(dpsi_x[n]);
        let xi_n_x = xi_x[n];
        let dxi_n_x = dxi_x[n];

        let a_num = m * psi_n_mx * dpsi_n_x - psi_n_x * dpsi_n_mx;
        let a_den = m * psi_n_mx * dxi_n_x - xi_n_x * dpsi_n_mx;
        a_coeffs.push(a_num / a_den);

        // b_n = (psi_n(mx) * psi_n'(x) - m * psi_n(x) * psi_n'(mx))
        //     / (psi_n(mx) * xi_n'(x)  - m * xi_n(x)  * psi_n'(mx))
        let b_num = psi_n_mx * dpsi_n_x - m * psi_n_x * dpsi_n_mx;
        let b_den = psi_n_mx * dxi_n_x - m * xi_n_x * dpsi_n_mx;
        b_coeffs.push(b_num / b_den);
    }

    (a_coeffs, b_coeffs)
}

// --- Riccati-Bessel functions ---

/// Riccati-Bessel function psi_n(z) = z * j_n(z) for complex argument.
/// Computed via upward recurrence from psi_0 = sin(z), psi_1 = sin(z)/z - cos(z).
fn riccati_bessel_psi_complex(z: Complex64, n_max: usize) -> Vec<Complex64> {
    let mut psi = vec![Complex64::from(0.0); n_max + 1];
    psi[0] = z.sin();
    if n_max == 0 {
        return psi;
    }
    psi[1] = z.sin() / z - z.cos();

    for n in 1..n_max {
        let nf = Complex64::from((2 * n + 1) as f64);
        psi[n + 1] = nf / z * psi[n] - psi[n - 1];
    }

    psi
}

/// Derivative psi_n'(z) via the recurrence relation:
/// psi_n'(z) = psi_{n-1}(z) - n/z * psi_n(z)
fn riccati_bessel_dpsi_complex(z: Complex64, n_max: usize, psi: &[Complex64]) -> Vec<Complex64> {
    let mut dpsi = vec![Complex64::from(0.0); n_max + 1];
    dpsi[0] = z.cos(); // psi_0' = cos(z)

    for n in 1..=n_max {
        let nf = Complex64::from(n as f64);
        dpsi[n] = psi[n - 1] - nf / z * psi[n];
    }

    dpsi
}

/// Riccati-Bessel function psi_n(x) for real argument.
fn riccati_bessel_psi_real(x: f64, n_max: usize) -> Vec<f64> {
    let mut psi = vec![0.0; n_max + 1];
    psi[0] = x.sin();
    if n_max == 0 {
        return psi;
    }
    psi[1] = x.sin() / x - x.cos();

    for n in 1..n_max {
        psi[n + 1] = (2 * n + 1) as f64 / x * psi[n] - psi[n - 1];
    }

    psi
}

/// Derivative psi_n'(x) for real argument.
fn riccati_bessel_dpsi_real(x: f64, n_max: usize, psi: &[f64]) -> Vec<f64> {
    let mut dpsi = vec![0.0; n_max + 1];
    dpsi[0] = x.cos();

    for n in 1..=n_max {
        dpsi[n] = psi[n - 1] - n as f64 / x * psi[n];
    }

    dpsi
}

/// Riccati-Bessel function xi_n(x) = x * h_n^(1)(x) = psi_n(x) + i * chi_n(x)
/// for real argument. chi_n(x) = -x * y_n(x) where y_n is the spherical Neumann function.
fn riccati_bessel_xi_real(x: f64, n_max: usize) -> Vec<Complex64> {
    let psi = riccati_bessel_psi_real(x, n_max);
    let chi = riccati_bessel_chi_real(x, n_max);

    let mut xi = vec![Complex64::from(0.0); n_max + 1];
    for n in 0..=n_max {
        xi[n] = Complex64::new(psi[n], chi[n]);
    }
    xi
}

/// chi_n(x) = -x * y_n(x) (Riccati-Bessel of the second kind for real argument).
fn riccati_bessel_chi_real(x: f64, n_max: usize) -> Vec<f64> {
    let mut chi = vec![0.0; n_max + 1];
    chi[0] = -x.cos();
    if n_max == 0 {
        return chi;
    }
    chi[1] = -x.cos() / x - x.sin();

    for n in 1..n_max {
        chi[n + 1] = (2 * n + 1) as f64 / x * chi[n] - chi[n - 1];
    }

    chi
}

/// Derivative xi_n'(x) for real argument.
fn riccati_bessel_dxi_real(x: f64, n_max: usize, xi: &[Complex64]) -> Vec<Complex64> {
    let mut dxi = vec![Complex64::from(0.0); n_max + 1];
    // xi_0'(x) = cos(x) + i * sin(x) = exp(ix)
    dxi[0] = Complex64::new(x.cos(), x.sin());

    for n in 1..=n_max {
        let nf = Complex64::from(n as f64);
        let xc = Complex64::from(x);
        dxi[n] = xi[n - 1] - nf / xc * xi[n];
    }

    dxi
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mie_small_dielectric_sphere() {
        // A small dielectric sphere in the Rayleigh limit should have
        // C_ext ~ C_sca (negligible absorption) and C_sca ~ x^4.
        let radius = 10.0; // nm
        let wavelength = 600.0; // nm
        let epsilon = Complex64::new(4.0, 0.0); // lossless dielectric, n=2
        let n_medium = 1.0;

        let (c_ext, c_sca, c_abs) = mie_cross_sections(radius, wavelength, epsilon, n_medium);

        assert!(c_ext > 0.0, "Extinction must be positive");
        assert!(c_sca > 0.0, "Scattering must be positive");
        assert!(c_abs.abs() < c_ext * 0.01, "Lossless sphere should have negligible absorption");
        assert!(
            (c_ext - c_sca).abs() / c_ext < 0.01,
            "C_ext ≈ C_sca for lossless sphere"
        );
    }

    #[test]
    fn test_mie_metallic_sphere_has_absorption() {
        // Gold-like sphere: should have significant absorption near the plasmon resonance.
        let radius = 20.0;
        let wavelength = 520.0;
        let epsilon = Complex64::new(-4.0, 2.0);
        let n_medium = 1.0;

        let (c_ext, c_sca, c_abs) = mie_cross_sections(radius, wavelength, epsilon, n_medium);

        assert!(c_ext > 0.0);
        assert!(c_abs > 0.0, "Metallic sphere must have absorption");
        assert!(c_sca >= 0.0);
        assert!(
            (c_ext - c_sca - c_abs).abs() / c_ext < 1e-10,
            "C_ext = C_sca + C_abs"
        );
    }

    #[test]
    fn test_mie_known_rayleigh_limit() {
        // In the Rayleigh limit (x << 1), the extinction cross-section is:
        // C_ext = 8*pi/3 * k^4 * a^6 * |K|^2 where K = (eps-1)/(eps+2) for n_m=1
        // We check that Mie agrees with this for a very small sphere.
        let radius = 1.0; // nm, very small
        let wavelength = 500.0;
        let epsilon = Complex64::new(2.25, 0.0); // n = 1.5
        let n_medium = 1.0;

        let k = 2.0 * std::f64::consts::PI * n_medium / wavelength;
        let x = k * radius;

        let (c_ext, _, _) = mie_cross_sections(radius, wavelength, epsilon, n_medium);

        // Rayleigh formula for scattering (dominant for small lossless sphere)
        let k_factor = (epsilon - 1.0) / (epsilon + 2.0);
        let c_rayleigh = 8.0 / 3.0 * std::f64::consts::PI * k.powi(4) * radius.powi(6) * k_factor.norm_sqr();

        // Should agree to within a few percent for x << 1
        assert!(
            (c_ext - c_rayleigh).abs() / c_rayleigh < 0.05,
            "Mie ({:.6e}) should match Rayleigh ({:.6e}) for x={:.4e}",
            c_ext,
            c_rayleigh,
            x
        );
    }
}
