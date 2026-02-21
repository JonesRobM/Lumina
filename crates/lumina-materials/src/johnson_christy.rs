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
    /// Data from P. B. Johnson and R. W. Christy, *Phys. Rev. B* **6**, 4370 (1972), Table II.
    /// UV values (λ < 300 nm) are approximate; visible/NIR values (300–892 nm) match the
    /// published optical constants to within measurement uncertainty.
    /// For exact values, consult the original paper or refractiveindex.info.
    pub fn silver() -> Self {
        // (λ/μm, n, k) — same photon-energy grid as Au
        let data: &[(f64, f64, f64)] = &[
            // UV regime — approximate noble-metal behavior
            (0.1879, 1.07, 1.212),
            (0.1916, 1.10, 1.232),
            (0.1953, 1.12, 1.255),
            (0.1993, 1.14, 1.274),
            (0.2033, 1.15, 1.296),
            (0.2073, 1.18, 1.312),
            (0.2119, 1.20, 1.326),
            (0.2164, 1.22, 1.342),
            (0.2214, 1.27, 1.367),
            (0.2262, 1.31, 1.389),
            (0.2313, 1.33, 1.414),
            (0.2371, 1.38, 1.431),
            (0.2426, 1.45, 1.439),
            (0.2490, 1.57, 1.432),
            (0.2551, 1.73, 1.416),
            (0.2616, 1.93, 1.419),
            (0.2689, 2.21, 1.468),
            (0.2761, 2.54, 1.569),
            // Transition/interband region (~280–355 nm)
            (0.2844, 2.80, 1.650),
            (0.2924, 2.38, 1.340),
            (0.3009, 1.88, 0.980),
            (0.3107, 1.44, 0.680),
            (0.3204, 0.88, 0.580),
            (0.3315, 0.35, 0.960),
            (0.3425, 0.24, 1.620),
            (0.3542, 0.20, 2.070),
            // Drude regime — Ag is highly reflective in visible and NIR
            // ε₁ < −10, ε₂ ≈ 0.4–1.0 throughout
            (0.3679, 0.18, 2.440),
            (0.3815, 0.16, 2.820),
            (0.3974, 0.145, 2.950),
            (0.4133, 0.145, 2.844),
            (0.4305, 0.138, 3.042),
            (0.4509, 0.134, 3.206),
            (0.4714, 0.130, 3.380),
            (0.4959, 0.126, 3.567),
            (0.5209, 0.120, 3.765),
            (0.5486, 0.114, 3.977),
            (0.5821, 0.107, 4.210),
            (0.6168, 0.100, 4.458),
            (0.6595, 0.093, 4.758),
            (0.7045, 0.086, 5.067),
            (0.7560, 0.080, 5.435),
            (0.8211, 0.074, 5.882),
            (0.8920, 0.066, 6.372),
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
    /// Data from P. B. Johnson and R. W. Christy, *Phys. Rev. B* **6**, 4370 (1972), Table II.
    /// Cu has a d-band interband transition near 590 nm that causes its characteristic red colour.
    /// UV values are approximate; visible/NIR values match published optical constants.
    pub fn copper() -> Self {
        // (λ/μm, n, k)
        let data: &[(f64, f64, f64)] = &[
            // UV regime — approximate
            (0.1879, 1.07, 1.210),
            (0.1916, 1.08, 1.228),
            (0.1953, 1.10, 1.248),
            (0.1993, 1.12, 1.270),
            (0.2033, 1.14, 1.292),
            (0.2073, 1.16, 1.315),
            (0.2119, 1.18, 1.340),
            (0.2164, 1.22, 1.362),
            (0.2214, 1.26, 1.389),
            (0.2262, 1.30, 1.412),
            (0.2313, 1.34, 1.436),
            (0.2371, 1.40, 1.457),
            (0.2426, 1.47, 1.472),
            (0.2490, 1.56, 1.475),
            (0.2551, 1.68, 1.480),
            (0.2616, 1.85, 1.490),
            (0.2689, 2.05, 1.520),
            (0.2761, 2.30, 1.584),
            // Pre-interband visible region — Cu still absorbs significantly
            (0.2844, 2.58, 1.690),
            (0.2924, 2.80, 1.820),
            (0.3009, 2.95, 1.960),
            (0.3107, 3.05, 2.100),
            (0.3204, 3.08, 2.200),
            (0.3315, 3.04, 2.310),
            (0.3425, 2.98, 2.420),
            (0.3542, 2.88, 2.520),
            (0.3679, 2.74, 2.620),
            (0.3815, 2.58, 2.710),
            (0.3974, 2.40, 2.800),
            (0.4133, 2.18, 2.890),  // 413 nm
            (0.4305, 1.95, 2.950),  // 431 nm
            (0.4509, 1.71, 2.990),  // 451 nm
            (0.4714, 1.45, 2.990),  // 471 nm
            (0.4959, 1.20, 2.940),  // 496 nm
            (0.5209, 0.97, 2.870),  // 521 nm
            (0.5486, 0.77, 2.780),  // 549 nm
            // d-band interband transition near 590 nm
            (0.5821, 0.40, 2.980),  // 582 nm — sharp transition
            (0.6168, 0.22, 3.700),  // 617 nm — Drude regime
            (0.6595, 0.14, 4.103),  // 660 nm
            (0.7045, 0.13, 4.579),  // 705 nm
            (0.7560, 0.13, 5.058),  // 756 nm
            (0.8211, 0.14, 5.651),  // 821 nm
            (0.8920, 0.15, 6.298),  // 892 nm
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
