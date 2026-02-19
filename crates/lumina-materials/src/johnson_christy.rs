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
    pub fn gold() -> Self {
        // TODO: Embed the actual J&C data for gold.
        // Placeholder with a few representative points.
        let wavelengths = vec![400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0, 800.0];
        let eps_real = vec![-1.66, -1.95, -2.83, -5.42, -9.64, -14.05, -18.47, -23.08, -28.24];
        let eps_imag = vec![5.29, 4.01, 3.07, 2.26, 1.66, 1.36, 1.18, 1.10, 1.08];

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
