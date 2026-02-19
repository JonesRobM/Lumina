//! Coupled Dipole Approximation (CDA) solver.
//!
//! This module implements the CDA, which models the optical response of a
//! nanostructure as an ensemble of $N$ polarisable point dipoles interacting
//! via the full dyadic Green's function. The self-consistent dipole moments
//! are obtained by solving a $3N \times 3N$ complex linear system.
//!
//! # Method selection
//!
//! - **Direct solve** (LU decomposition via `faer`): Used when $N \leq 1000$.
//!   Exact but $O(N^3)$ in time and $O(N^2)$ in memory.
//! - **Iterative solve** (GMRES): Used when $N > 1000$. Requires only
//!   matrix-vector products, reducing memory to $O(N)$ per iteration.

pub mod assembly;
pub mod direct;
pub mod greens;
pub mod iterative;

use super::{NearFieldPlane, OpticalSolver, SolverError};
use crate::types::{CrossSections, Dipole, DipoleResponse, NearFieldMap};
use crate::types::SimulationParams;

/// The CDA solver, holding configuration for the numerical method.
pub struct CdaSolver {
    /// Threshold number of dipoles above which the iterative solver is used.
    pub iterative_threshold: usize,
}

impl Default for CdaSolver {
    fn default() -> Self {
        Self {
            iterative_threshold: 1000,
        }
    }
}

impl CdaSolver {
    pub fn new(iterative_threshold: usize) -> Self {
        Self {
            iterative_threshold,
        }
    }
}

impl OpticalSolver for CdaSolver {
    fn compute_cross_sections(
        &self,
        _dipoles: &[Dipole],
        _wavelength_nm: f64,
        _params: &SimulationParams,
    ) -> Result<CrossSections, SolverError> {
        // TODO: Implement — assemble interaction matrix, solve, compute
        // extinction/absorption/scattering from dipole moments.
        todo!("CDA cross-section computation")
    }

    fn solve_dipoles(
        &self,
        _dipoles: &[Dipole],
        _wavelength_nm: f64,
        _params: &SimulationParams,
    ) -> Result<DipoleResponse, SolverError> {
        // TODO: Implement — core CDA linear system solve.
        todo!("CDA dipole solve")
    }

    fn compute_near_field(
        &self,
        _dipoles: &[Dipole],
        _response: &DipoleResponse,
        _plane: &NearFieldPlane,
    ) -> Result<NearFieldMap, SolverError> {
        // TODO: Implement — sum dipole field contributions on observation grid.
        todo!("CDA near-field computation")
    }

    fn method_name(&self) -> &str {
        "Coupled Dipole Approximation (CDA)"
    }
}
