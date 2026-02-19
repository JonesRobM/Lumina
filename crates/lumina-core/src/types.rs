//! Core types shared across the LuminaCDA framework.
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
