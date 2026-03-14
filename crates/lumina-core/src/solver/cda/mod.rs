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
//! - **Iterative solve** (GMRES): Used when $N > 1000$. The sub-strategy is
//!   chosen by comparing the matrix size $(3N)^2 \times 16$ bytes against
//!   `matrix_memory_budget` (default 2 GiB):
//!   - **Within budget + GPU**: assemble once on CPU, upload to GPU, GMRES
//!     matvec on GPU (fast).
//!   - **Within budget, CPU only**: assemble once, BLAS matvec reused across
//!     all GMRES iterations.
//!   - **Exceeds budget**: on-the-fly matvec — $O(N)$ peak memory, no
//!     assembly.  Krylov vectors only (≈ 90 × 3N complex doubles).
//!
//!   The budget gate is checked *before* attempting GPU assembly, so large
//!   systems always take the on-the-fly path regardless of GPU availability.

pub mod assembly;
pub mod direct;
pub mod fft_matvec;
pub mod greens;
pub mod iterative;

use std::sync::Arc;

use ndarray::Array2;
use num_complex::Complex64;

use super::{NearFieldPlane, OpticalSolver, SolverError};
use super::substrate::SubstrateSpec;
use crate::types::{CrossSections, Dipole, DipoleResponse, FarFieldMap, IncidentField, NearFieldMap, SimulationParams};
use lumina_compute::{ComputeBackend, CpuBackend};
use lumina_geometry::lattice::LatticeSpec;

/// The CDA solver, holding configuration for the numerical method.
pub struct CdaSolver {
    /// Threshold number of dipoles above which the iterative solver is used.
    pub iterative_threshold: usize,
    /// Incident field specification.
    pub incident_field: IncidentField,
    /// Whether to use the Filtered Coupled Dipole (FCD) Green's function.
    /// FCD volume-averages the Green's tensor over each source cell, smoothing
    /// the 1/R³ near-field singularity that causes staircase artefacts for
    /// metallic particles. Ignored for periodic assemblies.
    pub use_fcd: bool,
    /// Dipole lattice spacing (nm). Required for FCD integration volume.
    pub cell_size: f64,
    /// Compute backend for matrix-vector products and parallel operations.
    pub backend: Arc<dyn ComputeBackend>,
    /// Periodic lattice specification. When `Some`, the Ewald periodic
    /// Green's function replaces the free-space dyadic Green's function.
    pub lattice: Option<LatticeSpec>,
    /// Planar substrate specification. When `Some`, image-dipole corrections
    /// are added to the interaction matrix at each wavelength using the
    /// pre-computed Fresnel factor from `SimulationParams::substrate_delta_eps`.
    pub substrate: Option<SubstrateSpec>,
    /// Memory budget for the interaction matrix (bytes).
    ///
    /// For large systems on the CPU GMRES path, the solver compares the matrix
    /// size `(3N)² × 16` against this budget:
    ///
    /// - **Within budget**: assemble the full matrix once and use BLAS for
    ///   each matvec (fast — memory is O(N²), but reused across iterations).
    /// - **Exceeds budget**: compute A·x on-the-fly per GMRES iteration
    ///   (O(N) memory, same O(N²) work per call).
    ///
    /// Default: 2 GiB (fits matrices up to N ≈ 3 860).
    /// Set to `0` to always use on-the-fly.
    pub matrix_memory_budget: u64,
}

/// Default matrix memory budget: 2 GiB.
const DEFAULT_MATRIX_MEMORY_BUDGET: u64 = 2 * 1024 * 1024 * 1024;

impl Default for CdaSolver {
    fn default() -> Self {
        Self {
            iterative_threshold: 1000,
            incident_field: IncidentField::default(),
            use_fcd: true,
            cell_size: 3.0,
            backend: Arc::new(CpuBackend::new()),
            lattice: None,
            substrate: None,
            matrix_memory_budget: DEFAULT_MATRIX_MEMORY_BUDGET,
        }
    }
}

impl CdaSolver {
    pub fn new(iterative_threshold: usize) -> Self {
        Self { iterative_threshold, ..Default::default() }
    }

    /// Create a CDA solver with explicit FCD configuration.
    pub fn with_fcd(iterative_threshold: usize, use_fcd: bool, cell_size: f64) -> Self {
        Self { iterative_threshold, use_fcd, cell_size, ..Default::default() }
    }

    /// Create a CDA solver with a custom incident field (polarisation/direction).
    pub fn with_incident(iterative_threshold: usize, use_fcd: bool, cell_size: f64, incident_field: IncidentField) -> Self {
        Self { iterative_threshold, use_fcd, cell_size, incident_field, ..Default::default() }
    }

    /// Create a CDA solver with a specific compute backend.
    pub fn with_backend(iterative_threshold: usize, use_fcd: bool, cell_size: f64, backend: Arc<dyn ComputeBackend>) -> Self {
        Self { iterative_threshold, use_fcd, cell_size, backend, ..Default::default() }
    }

    /// Compute cross-sections from a pre-solved [`DipoleResponse`].
    ///
    /// Avoids re-solving the linear system when the dipole response is already
    /// available — for example when also computing the SHG response at the same
    /// wavelength.
    pub fn cross_sections_from_response(
        &self,
        dipoles: &[Dipole],
        response: &DipoleResponse,
        params: &SimulationParams,
    ) -> CrossSections {
        let k = Self::wavenumber(response.wavelength_nm, params.environment_n);
        let e0_sq = self.incident_field.amplitude * self.incident_field.amplitude;
        let n = dipoles.len();

        let mut c_ext = 0.0;
        for i in 0..n {
            let e_inc = self.incident_field.at_position(&dipoles[i].position, k);
            for c in 0..3 {
                c_ext += (e_inc[c].conj() * response.moments[[i, c]]).im;
            }
        }
        c_ext *= k / e0_sq;

        let mut c_abs = 0.0;
        let k3 = k.powi(3);
        for i in 0..n {
            let inv_alpha = assembly::invert_3x3_pub(&dipoles[i].polarisability);
            let mut pa_inv_alpha_conj_p = Complex64::from(0.0);
            for a in 0..3 {
                for b in 0..3 {
                    pa_inv_alpha_conj_p += response.moments[[i, a]]
                        * inv_alpha[3 * a + b].conj()
                        * response.moments[[i, b]].conj();
                }
            }
            let p_sq: f64 = (0..3).map(|c| response.moments[[i, c]].norm_sqr()).sum();
            c_abs += pa_inv_alpha_conj_p.im - (2.0 / 3.0) * k3 * p_sq;
        }
        c_abs *= k / e0_sq;

        CrossSections {
            wavelength_nm: response.wavelength_nm,
            extinction: c_ext,
            absorption: c_abs,
            scattering: (c_ext - c_abs).max(0.0),
            circular_dichroism: None,
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
            circular_dichroism: None,
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

        // Substrate runtime — set per-wavelength by runner/GUI in SimulationParams.
        let substrate = params.substrate_runtime.as_ref();

        // Build the RHS (incident field at each dipole)
        let rhs = assembly::build_incident_field_vector(dipoles, &self.incident_field, k);

        // Solve: use direct for small systems, GMRES for large.
        //
        // The iterative path has three sub-strategies, chosen by comparing the
        // matrix size (3N)²·16 bytes against `matrix_memory_budget`:
        //
        // - **Within budget, GPU available**: assemble once on CPU, upload to
        //   GPU, GMRES matvec on GPU (fast).
        // - **Within budget, CPU only**: assemble once on CPU, BLAS matvec
        //   (fast — reuses matrix across all GMRES iterations).
        // - **Exceeds budget**: on-the-fly matvec — evaluates G(rᵢ,rⱼ)·xⱼ
        //   per GMRES iteration without storing the matrix; O(N) peak memory.
        //
        // The budget gate is applied *before* attempting GPU to avoid a silent
        // O(N²) CPU allocation on the assembly path for large systems.
        let matrix_bytes = (3 * n as u64).pow(2) * 16;

        let solution = if n <= self.iterative_threshold {
            // Small system: assemble matrix and solve directly (LU).
            let matrix = assembly::assemble_interaction_matrix(
                dipoles,
                k,
                self.use_fcd,
                self.cell_size,
                self.lattice.as_ref(),
                params.k_bloch,
                substrate,
            );
            direct::solve_direct(&matrix, &rhs)
                .map_err(|e| SolverError::LinAlgError(e.to_string()))?
        } else if self.lattice.is_none()
            && params.substrate_runtime.is_none()
            && !self.use_fcd
        {
            // FFT path: O(N log N) for regular cubic grids; falls back to on-the-fly
            // if detect_regular_grid returns None (non-uniform or irregular spacing).
            if let Some(plan) = fft_matvec::FftMatvecPlan::try_build(dipoles, self.cell_size, k) {
                let matvec_fn = |x: &ndarray::Array1<Complex64>| -> Result<ndarray::Array1<Complex64>, SolverError> {
                    Ok(plan.matvec(dipoles, x.view()))
                };
                iterative::solve_gmres(&matvec_fn, &rhs, params.solver_tolerance, params.max_iterations)?
            } else {
                // Grid not regular — fall through to on-the-fly.
                let lattice_ref = self.lattice.as_ref();
                let matvec_fn = |x: &ndarray::Array1<Complex64>| -> Result<ndarray::Array1<Complex64>, SolverError> {
                    Ok(assembly::matvec_on_the_fly(
                        dipoles,
                        x,
                        k,
                        self.use_fcd,
                        self.cell_size,
                        lattice_ref,
                        params.k_bloch,
                        substrate,
                    ))
                };
                iterative::solve_gmres(&matvec_fn, &rhs, params.solver_tolerance, params.max_iterations)?
            }
        } else if matrix_bytes <= self.matrix_memory_budget {
            // Matrix fits in budget: assemble once and choose GPU or CPU matvec.
            let matrix = assembly::assemble_interaction_matrix(
                dipoles,
                k,
                self.use_fcd,
                self.cell_size,
                self.lattice.as_ref(),
                params.k_bloch,
                substrate,
            );
            if let Some(session) = self.backend.create_session(3 * n) {
                // GPU path: upload assembled matrix, GMRES matvec on GPU.
                session
                    .upload_matrix(&matrix)
                    .map_err(|e| SolverError::LinAlgError(e.to_string()))?;
                let matvec_fn = move |x: &ndarray::Array1<Complex64>| -> Result<ndarray::Array1<Complex64>, SolverError> {
                    session.matvec(x.view()).map_err(|e| SolverError::LinAlgError(e.to_string()))
                };
                iterative::solve_gmres(&matvec_fn, &rhs, params.solver_tolerance, params.max_iterations)?
            } else {
                // CPU cached path: BLAS matvec reuses matrix across iterations.
                let backend = Arc::clone(&self.backend);
                let matvec_fn = move |x: &ndarray::Array1<Complex64>| -> Result<ndarray::Array1<Complex64>, SolverError> {
                    backend.matvec(&matrix, x).map_err(|e| SolverError::LinAlgError(e.to_string()))
                };
                iterative::solve_gmres(&matvec_fn, &rhs, params.solver_tolerance, params.max_iterations)?
            }
        } else {
            // Matrix exceeds budget: on-the-fly matvec — O(N) memory, no assembly.
            let lattice_ref = self.lattice.as_ref();
            let matvec_fn = |x: &ndarray::Array1<Complex64>| -> Result<ndarray::Array1<Complex64>, SolverError> {
                Ok(assembly::matvec_on_the_fly(
                    dipoles,
                    x,
                    k,
                    self.use_fcd,
                    self.cell_size,
                    lattice_ref,
                    params.k_bloch,
                    substrate,
                ))
            };
            iterative::solve_gmres(&matvec_fn, &rhs, params.solver_tolerance, params.max_iterations)?
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

    fn compute_far_field(
        &self,
        dipoles: &[Dipole],
        response: &DipoleResponse,
        n_theta: usize,
        n_phi: usize,
    ) -> FarFieldMap {
        let k = Self::wavenumber(response.wavelength_nm, 1.0);
        crate::fields::compute_far_field(dipoles, response, k, n_theta, n_phi)
    }

    fn method_name(&self) -> &str {
        "Coupled Dipole Approximation (CDA)"
    }
}
