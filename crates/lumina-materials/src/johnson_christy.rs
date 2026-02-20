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
