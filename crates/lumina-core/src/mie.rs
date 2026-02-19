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
/// * `wavelength_nm` - Incident wavelength in nanometres.
/// * `epsilon_sphere` - Complex dielectric function of the sphere at this wavelength.
/// * `n_medium` - Refractive index of the surrounding medium (real).
///
/// # Returns
/// Tuple of (extinction, scattering, absorption) cross-sections in nm^2.
pub fn mie_cross_sections(
    _radius_nm: f64,
    _wavelength_nm: f64,
    _epsilon_sphere: Complex64,
    _n_medium: f64,
) -> (f64, f64, f64) {
    // TODO: Implement Mie coefficients a_n, b_n via Bessel/Hankel functions.
    // Truncate the series when |a_n|, |b_n| < tolerance.
    todo!("Mie theory cross-sections")
}

/// Compute individual Mie coefficients $a_n$ and $b_n$.
///
/// These are the electric and magnetic multipole coefficients of order $n$.
fn _mie_coefficients(
    _size_parameter: f64,
    _relative_refractive_index: Complex64,
    _n_max: usize,
) -> (Vec<Complex64>, Vec<Complex64>) {
    todo!("Mie coefficients a_n, b_n")
}
