//! Material property provider trait.
//!
//! All material data sources implement [`MaterialProvider`], which returns
//! frequency-dependent complex dielectric functions, polarisabilities, or
//! hyperpolarisabilities.

use num_complex::Complex64;
use thiserror::Error;

/// Errors from material providers.
#[derive(Debug, Error)]
pub enum MaterialError {
    #[error("Wavelength {wavelength_nm} nm is outside the data range [{min}, {max}] nm")]
    OutOfRange {
        wavelength_nm: f64,
        min: f64,
        max: f64,
    },

    #[error("Material not found: {0}")]
    NotFound(String),

    #[error("Data error: {0}")]
    DataError(String),
}

/// Provides frequency-dependent material properties.
///
/// Implementations include tabulated experimental data (Johnson & Christy,
/// Palik) and ab initio results (DFT polarisabilities).
pub trait MaterialProvider: Send + Sync {
    /// Human-readable name of this material.
    fn name(&self) -> &str;

    /// Wavelength range over which data is available (nm).
    fn wavelength_range(&self) -> (f64, f64);

    /// Complex dielectric function $\epsilon(\lambda)$ at a given wavelength.
    fn dielectric_function(&self, wavelength_nm: f64) -> Result<Complex64, MaterialError>;

    /// Complex refractive index $\tilde{n} = n + ik$ at a given wavelength.
    ///
    /// Default implementation derives from $\epsilon = \tilde{n}^2$.
    fn refractive_index(&self, wavelength_nm: f64) -> Result<Complex64, MaterialError> {
        let eps = self.dielectric_function(wavelength_nm)?;
        Ok(eps.sqrt())
    }
}
