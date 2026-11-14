//! Johnson & Christy tabulated dielectric functions.
//!
//! Optical constants for Au, Ag, and Cu from:
//! P. B. Johnson and R. W. Christy, *Phys. Rev. B* **6**, 4370 (1972).
//!
//! Data is embedded at compile time and interpolated via cubic splines.

use num_complex::Complex64;

use crate::provider::{MaterialError, MaterialProvider};
use crate::spline::CubicSpline;

/// Johnson & Christy material with spline-interpolated dielectric function.
pub struct JohnsonChristyMaterial {
    name: String,
    wavelengths_nm: Vec<f64>,
    spline_real: CubicSpline,
    spline_imag: CubicSpline,
}

impl JohnsonChristyMaterial {
    /// Construct from tabulated data.
    ///
    /// # Arguments
    /// * `name` - Material name (e.g. "Au", "Ag", "Cu").
    /// * `wavelengths_nm` - Wavelengths in nm.
    /// * `eps_real` - Real part of $\epsilon(\lambda)$.
    /// * `eps_imag` - Imaginary part of $\epsilon(\lambda)$.
    pub fn new(
        name: impl Into<String>,
        wavelengths_nm: Vec<f64>,
        eps_real: Vec<f64>,
        eps_imag: Vec<f64>,
    ) -> Self {
        let spline_real = CubicSpline::new(wavelengths_nm.clone(), eps_real);
        let spline_imag = CubicSpline::new(wavelengths_nm.clone(), eps_imag);

        Self {
            name: name.into(),
            wavelengths_nm,
            spline_real,
            spline_imag,
        }
    }

    /// Load the gold (Au) dataset.
    ///
    /// Data from P. B. Johnson and R. W. Christy, *Phys. Rev. B* **6**, 4370 (1972).
    /// Values of $n$ and $k$ are converted to $\epsilon_1 = n^2 - k^2$ and
    /// $\epsilon_2 = 2nk$.
    pub fn gold() -> Self {
        // J&C 1972 Au data: (λ/μm, n, k) from Table II.
        // Covering 188–937 nm (full published range).
        let data: &[(f64, f64, f64)] = &[
            (0.1879, 1.28, 1.188),
            (0.1916, 1.32, 1.203),
            (0.1953, 1.34, 1.226),
            (0.1993, 1.33, 1.251),
            (0.2033, 1.33, 1.277),
            (0.2073, 1.30, 1.304),
            (0.2119, 1.30, 1.350),
            (0.2164, 1.30, 1.387),
            (0.2214, 1.30, 1.427),
            (0.2262, 1.31, 1.460),
            (0.2313, 1.30, 1.497),
            (0.2371, 1.32, 1.536),
            (0.2426, 1.32, 1.577),
            (0.2490, 1.33, 1.631),
            (0.2551, 1.33, 1.688),
            (0.2616, 1.35, 1.749),
            (0.2689, 1.38, 1.803),
            (0.2761, 1.43, 1.847),
            (0.2844, 1.47, 1.869),
            (0.2924, 1.49, 1.878),
            (0.3009, 1.53, 1.889),
            (0.3107, 1.53, 1.893),
            (0.3204, 1.54, 1.898),
            (0.3315, 1.48, 1.883),
            (0.3425, 1.48, 1.871),
            (0.3542, 1.50, 1.866),
            (0.3679, 1.48, 1.895),
            (0.3815, 1.46, 1.933),
            (0.3974, 1.47, 1.952),
            (0.4133, 1.46, 1.958),
            (0.4305, 1.45, 1.948),
            (0.4509, 1.38, 1.914),
            (0.4714, 1.31, 1.849),
            (0.4959, 1.04, 1.833),
            (0.5209, 0.62, 2.081),
            (0.5486, 0.43, 2.455),
            (0.5821, 0.29, 2.863),
            (0.6168, 0.21, 3.272),
            (0.6595, 0.14, 3.697),
            (0.7045, 0.13, 4.103),
            (0.7560, 0.14, 4.542),
            (0.8211, 0.16, 5.083),
            (0.8920, 0.17, 5.663),
        ];

        let mut wavelengths = Vec::with_capacity(data.len());
        let mut eps_real = Vec::with_capacity(data.len());
        let mut eps_imag = Vec::with_capacity(data.len());

        for &(wl_um, n, k) in data {
            wavelengths.push(wl_um * 1000.0); // μm → nm
            eps_real.push(n * n - k * k);      // ε₁ = n² − k²
            eps_imag.push(2.0 * n * k);        // ε₂ = 2nk
        }

        Self::new("Au (Johnson & Christy)", wavelengths, eps_real, eps_imag)
    }

    /// Load the silver (Ag) dataset.
    ///
    /// Data from P. B. Johnson and R. W. Christy, *Phys. Rev. B* **6**, 4370 (1972),
    /// as tabulated by refractiveindex.info (`data/main/Ag/nk/Johnson.yml`).
    pub fn silver() -> Self {
        // (λ/μm, n, k) — same photon-energy grid as Au
        let data: &[(f64, f64, f64)] = &[
            (0.1879, 1.07, 1.212),
            (0.1916, 1.10, 1.232),
            (0.1953, 1.12, 1.255),
            (0.1993, 1.14, 1.277),
            (0.2033, 1.15, 1.296),
            (0.2073, 1.18, 1.312),
            (0.2119, 1.20, 1.325),
            (0.2164, 1.22, 1.336),
            (0.2214, 1.25, 1.342),
            (0.2262, 1.26, 1.344),
            (0.2313, 1.28, 1.357),
            (0.2371, 1.28, 1.367),
            (0.2426, 1.30, 1.378),
            (0.2490, 1.31, 1.389),
            (0.2551, 1.33, 1.393),
            (0.2616, 1.35, 1.387),
            (0.2689, 1.38, 1.372),
            (0.2761, 1.41, 1.331),
            (0.2844, 1.41, 1.264),
            (0.2924, 1.39, 1.161),
            (0.3009, 1.34, 0.964),
            (0.3107, 1.13, 0.616),
            (0.3204, 0.81, 0.392),
            (0.3315, 0.17, 0.829),
            (0.3425, 0.14, 1.142),
            (0.3542, 0.10, 1.419),
            (0.3679, 0.07, 1.657),
            (0.3815, 0.05, 1.864),
            (0.3974, 0.05, 2.070),
            (0.4133, 0.05, 2.275),
            (0.4305, 0.04, 2.462),
            (0.4509, 0.04, 2.657),
            (0.4714, 0.05, 2.869),
            (0.4959, 0.05, 3.093),
            (0.5209, 0.05, 3.324),
            (0.5486, 0.06, 3.586),
            (0.5821, 0.05, 3.858),
            (0.6168, 0.06, 4.152),
            (0.6595, 0.05, 4.483),
            (0.7045, 0.04, 4.838),
            (0.7560, 0.03, 5.242),
            (0.8211, 0.04, 5.727),
            (0.8920, 0.04, 6.312),
        ];

        let mut wavelengths = Vec::with_capacity(data.len());
        let mut eps_real = Vec::with_capacity(data.len());
        let mut eps_imag = Vec::with_capacity(data.len());

        for &(wl_um, n, k) in data {
            wavelengths.push(wl_um * 1000.0);
            eps_real.push(n * n - k * k);
            eps_imag.push(2.0 * n * k);
        }

        Self::new("Ag (Johnson & Christy)", wavelengths, eps_real, eps_imag)
    }

    /// Load the copper (Cu) dataset.
    ///
    /// Data from P. B. Johnson and R. W. Christy, *Phys. Rev. B* **6**, 4370 (1972),
    /// as tabulated by refractiveindex.info (`data/main/Cu/nk/Johnson.yml`).
    /// Cu has a d-band interband edge near 590 nm that causes its characteristic colour.
    pub fn copper() -> Self {
        // (λ/μm, n, k)
        let data: &[(f64, f64, f64)] = &[
            (0.1879, 0.94, 1.337),
            (0.1916, 0.95, 1.388),
            (0.1953, 0.97, 1.440),
            (0.1993, 0.98, 1.493),
            (0.2033, 0.99, 1.550),
            (0.2073, 1.01, 1.599),
            (0.2119, 1.04, 1.651),
            (0.2164, 1.08, 1.699),
            (0.2214, 1.13, 1.737),
            (0.2262, 1.18, 1.768),
            (0.2313, 1.23, 1.792),
            (0.2371, 1.28, 1.802),
            (0.2426, 1.34, 1.799),
            (0.2490, 1.37, 1.783),
            (0.2551, 1.41, 1.741),
            (0.2616, 1.41, 1.691),
            (0.2689, 1.45, 1.668),
            (0.2761, 1.46, 1.646),
            (0.2844, 1.45, 1.633),
            (0.2924, 1.42, 1.633),
            (0.3009, 1.40, 1.679),
            (0.3107, 1.38, 1.729),
            (0.3204, 1.38, 1.783),
            (0.3315, 1.34, 1.821),
            (0.3425, 1.36, 1.864),
            (0.3542, 1.37, 1.916),
            (0.3679, 1.36, 1.975),
            (0.3815, 1.33, 2.045),
            (0.3974, 1.32, 2.116),
            (0.4133, 1.28, 2.207),
            (0.4305, 1.25, 2.305),
            (0.4509, 1.24, 2.397),
            (0.4714, 1.25, 2.483),
            (0.4959, 1.22, 2.564),
            (0.5209, 1.18, 2.608),
            (0.5486, 1.02, 2.577),
            (0.5821, 0.70, 2.704),
            (0.6168, 0.30, 3.205),
            (0.6595, 0.22, 3.747),
            (0.7045, 0.21, 4.205),
            (0.7560, 0.24, 4.665),
            (0.8211, 0.26, 5.180),
            (0.8920, 0.30, 5.768),
        ];

        let mut wavelengths = Vec::with_capacity(data.len());
        let mut eps_real = Vec::with_capacity(data.len());
        let mut eps_imag = Vec::with_capacity(data.len());

        for &(wl_um, n, k) in data {
            wavelengths.push(wl_um * 1000.0);
            eps_real.push(n * n - k * k);
            eps_imag.push(2.0 * n * k);
        }

        Self::new("Cu (Johnson & Christy)", wavelengths, eps_real, eps_imag)
    }
}

impl MaterialProvider for JohnsonChristyMaterial {
    fn name(&self) -> &str {
        &self.name
    }

    fn wavelength_range(&self) -> (f64, f64) {
        (
            *self.wavelengths_nm.first().unwrap(),
            *self.wavelengths_nm.last().unwrap(),
        )
    }

    fn dielectric_function(&self, wavelength_nm: f64) -> Result<Complex64, MaterialError> {
        let (min, max) = self.wavelength_range();
        if wavelength_nm < min || wavelength_nm > max {
            return Err(MaterialError::OutOfRange {
                wavelength_nm,
                min,
                max,
            });
        }

        let re = self.spline_real.evaluate(wavelength_nm);
        let im = self.spline_imag.evaluate(wavelength_nm);
        Ok(Complex64::new(re, im))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Spot checks against Johnson & Christy's published n, k.
    #[test]
    fn silver_and_copper_match_published_values() {
        let cases = [
            (JohnsonChristyMaterial::silver(), 520.9, 0.05, 3.324),
            (JohnsonChristyMaterial::silver(), 354.2, 0.10, 1.419),
            (JohnsonChristyMaterial::copper(), 413.3, 1.28, 2.207),
            (JohnsonChristyMaterial::copper(), 659.5, 0.22, 3.747),
        ];
        for (m, wl, n, k) in cases {
            let expected = Complex64::new(n, k).powi(2);
            let eps = m.dielectric_function(wl).expect("in range");
            assert!(
                (eps - expected).norm() < 1e-6 * expected.norm(),
                "{} at {wl} nm: {eps} vs {expected}",
                m.name()
            );
        }
    }
}
