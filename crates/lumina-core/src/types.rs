//! Core types shared across the Lumina framework.
//!
//! This module defines the fundamental data structures used throughout the
//! simulation pipeline: dipoles, simulation parameters, and result containers.

use ndarray::Array2;
use num_complex::Complex64;
use serde::{Deserialize, Serialize};

/// A single point dipole in the simulation domain.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Dipole {
    /// Position in 3D space (nanometres).
    pub position: [f64; 3],
    /// Complex polarisability tensor (3x3) at the current frequency.
    /// Stored as a flat [9] array in row-major order.
    pub polarisability: [Complex64; 9],
}

impl Dipole {
    /// Create a dipole with isotropic (scalar) polarisability.
    pub fn isotropic(position: [f64; 3], alpha: Complex64) -> Self {
        let zero = Complex64::new(0.0, 0.0);
        Self {
            position,
            polarisability: [
                alpha, zero, zero,
                zero, alpha, zero,
                zero, zero, alpha,
            ],
        }
    }
}

/// Incident plane wave specification.
#[derive(Debug, Clone)]
pub struct IncidentField {
    /// Propagation direction (unit vector).
    pub direction: [f64; 3],
    /// Polarisation vector (unit vector, perpendicular to direction).
    pub polarisation: [f64; 3],
    /// Amplitude (V/m). Typically set to 1.0 for cross-section calculations.
    pub amplitude: f64,
}

impl Default for IncidentField {
    fn default() -> Self {
        Self {
            direction: [0.0, 0.0, 1.0],    // propagating along z
            polarisation: [1.0, 0.0, 0.0],  // x-polarised
            amplitude: 1.0,
        }
    }
}

impl IncidentField {
    /// Evaluate the incident field at a given position.
    ///
    /// $\mathbf{E}_{\text{inc}}(\mathbf{r}) = E_0 \hat{\mathbf{e}} \exp(i \mathbf{k} \cdot \mathbf{r})$
    pub fn at_position(&self, position: &[f64; 3], k: f64) -> [Complex64; 3] {
        let kdotr = k * (
            self.direction[0] * position[0]
            + self.direction[1] * position[1]
            + self.direction[2] * position[2]
        );
        let phase = Complex64::new(0.0, kdotr).exp();
        let e0 = Complex64::from(self.amplitude);
        [
            e0 * self.polarisation[0] * phase,
            e0 * self.polarisation[1] * phase,
            e0 * self.polarisation[2] * phase,
        ]
    }
}

/// Parameters defining a simulation run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationParams {
    /// Wavelength range [start, end] in nanometres.
    pub wavelength_range: [f64; 2],
    /// Number of wavelength points to evaluate.
    pub num_wavelengths: usize,
    /// Refractive index of the surrounding medium.
    pub environment_n: f64,
    /// Solver tolerance for iterative methods (e.g. GMRES).
    pub solver_tolerance: f64,
    /// Maximum iterations for iterative solvers.
    pub max_iterations: usize,
}

impl Default for SimulationParams {
    fn default() -> Self {
        Self {
            wavelength_range: [400.0, 900.0],
            num_wavelengths: 100,
            environment_n: 1.0,
            solver_tolerance: 1e-6,
            max_iterations: 1000,
        }
    }
}

/// Optical cross-sections at a single wavelength.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CrossSections {
    /// Wavelength (nm).
    pub wavelength_nm: f64,
    /// Extinction cross-section (nm^2).
    pub extinction: f64,
    /// Absorption cross-section (nm^2).
    pub absorption: f64,
    /// Scattering cross-section (nm^2).
    pub scattering: f64,
}

/// The solved dipole moments for a given excitation.
#[derive(Debug, Clone)]
pub struct DipoleResponse {
    /// Wavelength (nm) at which this was computed.
    pub wavelength_nm: f64,
    /// Complex dipole moments, shape (N, 3).
    pub moments: Array2<Complex64>,
    /// Local electric field at each dipole, shape (N, 3).
    pub local_fields: Array2<Complex64>,
}

/// Near-field map on a 2D plane.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NearFieldMap {
    /// Grid positions, shape (Nx * Ny, 3).
    pub positions: Vec<[f64; 3]>,
    /// |E|^2 values at each grid point.
    pub field_intensity: Vec<f64>,
    /// Number of points along x.
    pub nx: usize,
    /// Number of points along y.
    pub ny: usize,
    /// Spatial extent: [x_min, x_max, y_min, y_max] in nm.
    pub extent: [f64; 4],
}

/// Complete results from a simulation run.
#[derive(Debug, Clone)]
pub struct SimulationResult {
    /// Cross-section spectra across all wavelengths.
    pub spectra: Vec<CrossSections>,
    /// Dipole responses at each wavelength (optional, can be large).
    pub dipole_responses: Option<Vec<DipoleResponse>>,
    /// Near-field maps (optional).
    pub near_field_maps: Option<Vec<NearFieldMap>>,
}

/// Compute the Clausius-Mossotti polarisability for a small volume element.
///
/// $\alpha_{\text{CM}} = 3 V \epsilon_0 \frac{\epsilon - \epsilon_m}{\epsilon + 2\epsilon_m}$
///
/// Note: we work in Gaussian-like units where $\epsilon_0$ is absorbed, so the
/// returned polarisability has units of volume (nm^3).
///
/// # Arguments
/// * `volume_nm3` - Volume of the dipole element in nm^3 (typically $d^3$).
/// * `epsilon` - Complex dielectric function of the material.
/// * `epsilon_m` - Dielectric constant of the surrounding medium ($n_m^2$).
pub fn clausius_mossotti(volume_nm3: f64, epsilon: Complex64, epsilon_m: f64) -> Complex64 {
    let eps_m = Complex64::from(epsilon_m);
    3.0 * volume_nm3 * (epsilon - eps_m) / (epsilon + 2.0 * eps_m)
}

/// Apply the Draine radiative correction to the Clausius-Mossotti polarisability.
///
/// This ensures optical theorem consistency (energy conservation):
/// $\alpha_{\text{rad}} = \frac{\alpha_{\text{CM}}}{1 - \frac{i k^3}{6\pi} \alpha_{\text{CM}}}$
///
/// # Arguments
/// * `alpha_cm` - Clausius-Mossotti polarisability.
/// * `k` - Wavenumber in the medium (nm^{-1}).
pub fn radiative_correction(alpha_cm: Complex64, k: f64) -> Complex64 {
    let correction = Complex64::new(0.0, k.powi(3) / (6.0 * std::f64::consts::PI));
    alpha_cm / (Complex64::from(1.0) - correction * alpha_cm)
}
