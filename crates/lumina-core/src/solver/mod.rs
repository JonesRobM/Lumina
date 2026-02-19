//! Optical solver abstraction and implementations.
//!
//! The [`OpticalSolver`] trait defines the interface that all simulation
//! methods must implement. The CDA (Coupled Dipole Approximation) is the
//! first implementation; future methods (BEM, FDTD, T-matrix) will implement
//! the same trait.

pub mod cda;

use crate::types::{CrossSections, Dipole, DipoleResponse, NearFieldMap, SimulationParams};
use thiserror::Error;

/// Errors that can occur during an optical solve.
#[derive(Debug, Error)]
pub enum SolverError {
    #[error("Solver failed to converge after {max_iter} iterations (residual: {residual:.2e})")]
    ConvergenceFailure { max_iter: usize, residual: f64 },

    #[error("Invalid geometry: {0}")]
    InvalidGeometry(String),

    #[error("Linear algebra error: {0}")]
    LinAlgError(String),

    #[error("Compute backend error: {0}")]
    ComputeError(String),
}

/// The core trait that all optical simulation methods must implement.
///
/// This abstraction allows the GUI and CLI to operate against any solver
/// without knowledge of the underlying numerical method. For v0.1 the
/// only implementation is the CDA; future versions will add BEM, FDTD,
/// T-matrix, and other approaches.
pub trait OpticalSolver {
    /// Compute optical cross-sections (extinction, absorption, scattering)
    /// for the given dipole configuration at a single wavelength.
    fn compute_cross_sections(
        &self,
        dipoles: &[Dipole],
        wavelength_nm: f64,
        params: &SimulationParams,
    ) -> Result<CrossSections, SolverError>;

    /// Solve for the self-consistent dipole moments at a single wavelength.
    fn solve_dipoles(
        &self,
        dipoles: &[Dipole],
        wavelength_nm: f64,
        params: &SimulationParams,
    ) -> Result<DipoleResponse, SolverError>;

    /// Compute near-field intensity on a 2D observation plane.
    fn compute_near_field(
        &self,
        dipoles: &[Dipole],
        response: &DipoleResponse,
        plane: &NearFieldPlane,
    ) -> Result<NearFieldMap, SolverError>;

    /// Human-readable name of the solver method.
    fn method_name(&self) -> &str;
}

/// Defines a 2D observation plane for near-field calculations.
#[derive(Debug, Clone)]
pub struct NearFieldPlane {
    /// Centre of the plane (nm).
    pub centre: [f64; 3],
    /// Normal vector to the plane.
    pub normal: [f64; 3],
    /// Half-width of the plane (nm).
    pub half_width: f64,
    /// Half-height of the plane (nm).
    pub half_height: f64,
    /// Number of grid points along the width.
    pub nx: usize,
    /// Number of grid points along the height.
    pub ny: usize,
}
