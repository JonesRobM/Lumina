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

use ndarray::Array2;
use num_complex::Complex64;

use super::{NearFieldPlane, OpticalSolver, SolverError};
use crate::types::{CrossSections, Dipole, DipoleResponse, IncidentField, NearFieldMap, SimulationParams};

/// The CDA solver, holding configuration for the numerical method.
pub struct CdaSolver {
    /// Threshold number of dipoles above which the iterative solver is used.
    pub iterative_threshold: usize,
    /// Incident field specification.
    pub incident_field: IncidentField,
    /// Whether to use the Filtered Coupled Dipole (FCD) Green's function.
    /// FCD volume-averages the Green's tensor over each source cell, smoothing
    /// the 1/R³ near-field singularity that causes staircase artefacts for
    /// metallic particles.
    pub use_fcd: bool,
    /// Dipole lattice spacing (nm). Required for FCD integration volume.
    pub cell_size: f64,
}

impl Default for CdaSolver {
    fn default() -> Self {
        Self {
            iterative_threshold: 1000,
            incident_field: IncidentField::default(),
            use_fcd: true,
            cell_size: 3.0,
        }
    }
}

impl CdaSolver {
    pub fn new(iterative_threshold: usize) -> Self {
        Self {
            iterative_threshold,
            incident_field: IncidentField::default(),
            use_fcd: true,
            cell_size: 3.0,
        }
    }

    /// Create a CDA solver with explicit FCD configuration.
    pub fn with_fcd(iterative_threshold: usize, use_fcd: bool, cell_size: f64) -> Self {
        Self {
            iterative_threshold,
            incident_field: IncidentField::default(),
            use_fcd,
            cell_size,
        }
    }

    /// Compute the wavenumber in the medium from wavelength and refractive index.
    fn wavenumber(wavelength_nm: f64, n_medium: f64) -> f64 {
        2.0 * std::f64::consts::PI * n_medium / wavelength_nm
    }
}

impl OpticalSolver for CdaSolver {
    fn compute_cross_sections(
        &self,
        dipoles: &[Dipole],
        wavelength_nm: f64,
        params: &SimulationParams,
    ) -> Result<CrossSections, SolverError> {
        let response = self.solve_dipoles(dipoles, wavelength_nm, params)?;
        let k = Self::wavenumber(wavelength_nm, params.environment_n);
        let e0_sq = self.incident_field.amplitude * self.incident_field.amplitude;

        let n = dipoles.len();

        // Extinction: C_ext = (k / |E0|^2) * Sum_i Im(E_inc,i* . p_i)
        // Note: the 4π is already absorbed into the Green's function prefactor,
        // so the cross-section formula uses k/|E₀|² rather than 4πk/|E₀|².
        let mut c_ext = 0.0;
        for i in 0..n {
            let e_inc = self.incident_field.at_position(&dipoles[i].position, k);
            for c in 0..3 {
                c_ext += (e_inc[c].conj() * response.moments[[i, c]]).im;
            }
        }
        c_ext *= k / e0_sq;

        // Absorption: C_abs = (k / |E0|^2) * Sum_i [ Im(p_i . (alpha_i^{-1})* . p_i*) - (2/3)k^3 |p_i|^2 ]
        // Same 4π convention as extinction (absorbed into Green's function).
        let mut c_abs = 0.0;
        let k3 = k.powi(3);
        for i in 0..n {
            // Get inverse polarisability (alpha^{-1}) for this dipole
            let inv_alpha = assembly::invert_3x3_pub(&dipoles[i].polarisability);

            // Compute p_i . (alpha_i^{-1})* . p_i*
            // This is sum_{a,b} p_a * conj(inv_alpha_{ab}) * conj(p_b)
            let mut pa_inv_alpha_conj_p = Complex64::from(0.0);
            for a in 0..3 {
                for b in 0..3 {
                    pa_inv_alpha_conj_p += response.moments[[i, a]]
                        * inv_alpha[3 * a + b].conj()
                        * response.moments[[i, b]].conj();
                }
            }

            // |p_i|^2
            let p_sq: f64 = (0..3)
                .map(|c| response.moments[[i, c]].norm_sqr())
                .sum();

            c_abs += pa_inv_alpha_conj_p.im - (2.0 / 3.0) * k3 * p_sq;
        }
        c_abs *= k / e0_sq;

        let c_sca = c_ext - c_abs;

        Ok(CrossSections {
            wavelength_nm,
            extinction: c_ext,
            absorption: c_abs,
            scattering: c_sca.max(0.0), // Numerical noise can make this slightly negative
        })
    }

    fn solve_dipoles(
        &self,
        dipoles: &[Dipole],
        wavelength_nm: f64,
        params: &SimulationParams,
    ) -> Result<DipoleResponse, SolverError> {
        if dipoles.is_empty() {
            return Err(SolverError::InvalidGeometry("No dipoles provided".into()));
        }

        let k = Self::wavenumber(wavelength_nm, params.environment_n);
        let n = dipoles.len();

        // Assemble the interaction matrix
        let matrix = assembly::assemble_interaction_matrix(dipoles, k, self.use_fcd, self.cell_size);

        // Build the RHS (incident field at each dipole)
        let rhs = assembly::build_incident_field_vector(dipoles, &self.incident_field, k);

        // Solve: use direct for small systems, GMRES for large
        let solution = if n <= self.iterative_threshold {
            direct::solve_direct(&matrix, &rhs)
                .map_err(|e| SolverError::LinAlgError(e.to_string()))?
        } else {
            iterative::solve_gmres(
                &matrix,
                &rhs,
                params.solver_tolerance,
                params.max_iterations,
            )?
        };

        // Reshape the solution into (N, 3) dipole moments
        let mut moments = Array2::<Complex64>::zeros((n, 3));
        let mut local_fields = Array2::<Complex64>::zeros((n, 3));

        for i in 0..n {
            for c in 0..3 {
                moments[[i, c]] = solution[3 * i + c];
            }
            // Local field = E_inc + sum_j G_ij . p_j  (which equals alpha^{-1} . p_i)
            // More simply: E_loc,i = alpha_i^{-1} . p_i
            let inv_alpha = assembly::invert_3x3_pub(&dipoles[i].polarisability);
            for a in 0..3 {
                let mut e_loc = Complex64::from(0.0);
                for b in 0..3 {
                    e_loc += inv_alpha[3 * a + b] * moments[[i, b]];
                }
                local_fields[[i, a]] = e_loc;
            }
        }

        Ok(DipoleResponse {
            wavelength_nm,
            moments,
            local_fields,
        })
    }

    fn compute_near_field(
        &self,
        dipoles: &[Dipole],
        response: &DipoleResponse,
        plane: &NearFieldPlane,
    ) -> Result<NearFieldMap, SolverError> {
        let k = Self::wavenumber(
            response.wavelength_nm,
            1.0, // TODO: pass params through
        );

        crate::fields::compute_near_field_map(dipoles, response, plane, k, &self.incident_field)
    }

    fn method_name(&self) -> &str {
        "Coupled Dipole Approximation (CDA)"
    }
}
