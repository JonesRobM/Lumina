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
//! - **FFT matvec** (block-Toeplitz convolution via `rustfft`): Used when
//!   $N > 1000$, dipoles are on a regular cubic grid, and no periodic BC,
//!   substrate, or FCD correction is active. $O(N \log N)$ per matvec.
//!   When the grid is irregular the solver falls through to GPU or CPU assembly.
//! - **Iterative solve** (GMRES): Used when $N > 1000$ and FFT path is
//!   not available. The sub-strategy is chosen by comparing the matrix
//!   size $(3N)^2 \times 16$ bytes against `matrix_memory_budget` (default 2 GiB):
//!   - **Within budget + GPU assembly** (no lattice/substrate/FCD): assemble
//!     on GPU via compute shader, GMRES matvec on GPU (fast).
//!   - **Within budget, CPU only**: assemble once, BLAS matvec reused across
//!     all GMRES iterations.
//!   - **Exceeds budget**: on-the-fly matvec -- $O(N)$ peak memory, no
//!     assembly.  Krylov vectors only (~= 90 x 3N complex doubles).

// pub(crate): `assemble_interaction_matrix` and `matvec_on_the_fly` must not be
// reachable from outside the crate, or a caller could route around the
// substrate/lattice compatibility guard applied in `solve_dipoles`.
pub(crate) mod assembly;
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
    /// the 1/R^3 near-field singularity that causes staircase artefacts for
    /// metallic particles. Ignored for periodic assemblies.
    pub use_fcd: bool,
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
    /// size `(3N)^2 x 16` against this budget:
    ///
    /// - **Within budget**: assemble the full matrix once and use BLAS for
    ///   each matvec (fast -- memory is O(N^2), but reused across iterations).
    /// - **Exceeds budget**: compute A*x on-the-fly per GMRES iteration
    ///   (O(N) memory, same O(N^2) work per call).
    ///
    /// Default: 2 GiB (fits matrices up to N ~= 3 860).
    /// Set to `0` to always use on-the-fly.
    pub matrix_memory_budget: u64,
}

/// Default matrix memory budget: 2 GiB.
const DEFAULT_MATRIX_MEMORY_BUDGET: u64 = 2 * 1024 * 1024 * 1024;

impl Default for CdaSolver {
    fn default() -> Self {
        Self {
            // 3000 dipoles => a 1.3 GiB matrix, inside the 2 GiB default budget. The
            // direct solver beats GMRES by ~9-14x across the measured range, so the
            // crossover belongs at the memory limit, not at an arbitrary N.
            iterative_threshold: 3000,
            incident_field: IncidentField::default(),
            use_fcd: true,
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
    pub fn with_fcd(iterative_threshold: usize, use_fcd: bool) -> Self {
        Self { iterative_threshold, use_fcd, ..Default::default() }
    }

    /// Create a CDA solver with a custom incident field (polarisation/direction).
    pub fn with_incident(iterative_threshold: usize, use_fcd: bool, incident_field: IncidentField) -> Self {
        Self { iterative_threshold, use_fcd, incident_field, ..Default::default() }
    }

    /// Create a CDA solver with a specific compute backend.
    pub fn with_backend(iterative_threshold: usize, use_fcd: bool, backend: Arc<dyn ComputeBackend>) -> Self {
        Self { iterative_threshold, use_fcd, backend, ..Default::default() }
    }

    /// Compute cross-sections from a pre-solved [`DipoleResponse`].
    ///
    /// Avoids re-solving the linear system when the dipole response is already
    /// available -- for example when also computing the SHG response at the same
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
        // Note: the 4pi is already absorbed into the Green's function prefactor,
        // so the cross-section formula uses k/|E0|^2 rather than 4pik/|E0|^2.
        let mut c_ext = 0.0;
        for i in 0..n {
            let e_inc = self.incident_field.at_position(&dipoles[i].position, k);
            for c in 0..3 {
                c_ext += (e_inc[c].conj() * response.moments[[i, c]]).im;
            }
        }
        c_ext *= k / e0_sq;

        // Absorption: C_abs = (k / |E0|^2) * Sum_i [ Im(p_i . (alpha_i^{-1})* . p_i*) - (2/3)k^3 |p_i|^2 ]
        // Same 4pi convention as extinction (absorbed into Green's function).
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

        // Substrate runtime -- set per-wavelength by runner/GUI in SimulationParams.
        let substrate = params.substrate_runtime.as_ref();
        crate::solver::substrate::check_substrate_lattice_compatible(
            self.lattice.as_ref(),
            substrate,
        )?;

        // Build the RHS (incident field at each dipole)
        let rhs = assembly::build_incident_field_vector(dipoles, &self.incident_field, k);

        // Block-Jacobi preconditioning for every iterative path. A_ii = alpha^-1, so
        // M^-1 = diag(alpha) and the preconditioner costs one 3x3 multiply per dipole.
        // GMRES then solves M^-1A*p = M^-1E for the same p. Unpreconditioned, this
        // needed ~1000 matvecs at N ~= 1000 -- slower than a direct LU solve.
        let precond = assembly::jacobi_preconditioner(dipoles);
        let rhs_p = assembly::apply_preconditioner(&precond, &rhs);

        // Solve: use direct for small systems, GMRES for large.
        //
        // The iterative path has three sub-strategies, chosen by comparing the
        // matrix size (3N)^2*16 bytes against `matrix_memory_budget`:
        //
        // - **Within budget, GPU available**: assemble once on CPU, upload to
        //   GPU, GMRES matvec on GPU (fast).
        // - **Within budget, CPU only**: assemble once on CPU, BLAS matvec
        //   (fast -- reuses matrix across all GMRES iterations).
        // - **Exceeds budget**: on-the-fly matvec -- evaluates G(r_i,r_j)*x_j
        //   per GMRES iteration without storing the matrix; O(N) peak memory.
        //
        // The budget gate is applied *before* attempting GPU to avoid a silent
        // O(N^2) CPU allocation on the assembly path for large systems.
        let matrix_bytes = (3 * n as u64).pow(2) * 16;

        // Direct LU is much faster than GMRES wherever it fits: measured 8.6x at
        // N = 1021, 10.9x at 1357 and 14.4x at 1791 on an Au sphere, and the gap widens
        // with N. It is gated on the memory budget too, because unlike the iterative
        // paths it has no on-the-fly fallback -- without the gate, raising
        // `iterative_threshold` would risk a hard OOM on a large system.
        let solution = if n <= self.iterative_threshold && matrix_bytes <= self.matrix_memory_budget {
            // Small system: assemble matrix and solve directly (LU).
            let matrix = assembly::assemble_interaction_matrix(
                dipoles,
                k,
                self.use_fcd,
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
            // FFT path: O(N log N) for regular cubic grids; falls back to GPU
            // assembly or on-the-fly if the grid is irregular.
            // Note: the default `CdaSolver` has `use_fcd = true`, so this path
            // is only active when the solver is constructed with `use_fcd = false`.
            if let Some(plan) = fft_matvec::FftMatvecPlan::try_build(dipoles, dipoles[0].cell_size, k) {
                let matvec_fn = |x: &ndarray::Array1<Complex64>| -> Result<ndarray::Array1<Complex64>, SolverError> {
                    Ok(plan.matvec(dipoles, x.view()))
                };
                {
                    let preconditioned = |x: &ndarray::Array1<Complex64>| -> Result<ndarray::Array1<Complex64>, SolverError> {
                        Ok(assembly::apply_preconditioner(&precond, &matvec_fn(x)?))
                    };
                    iterative::solve_gmres(&preconditioned, &rhs_p, params.solver_tolerance, params.max_iterations)?
                }
            } else if matrix_bytes <= self.matrix_memory_budget {
                // Grid not regular -- try GPU assembly (conditions match: no FCD,
                // no substrate, no Ewald), then CPU assembly.
                let diagonal_blocks = assembly::compute_diagonal_blocks(dipoles);
                let positions: Vec<[f64; 3]> = dipoles.iter().map(|d| d.position).collect();
                let gpu_session = self.backend.assemble_and_create_session(&positions, k, &diagonal_blocks);

                if let Some(Ok(session)) = gpu_session {
                    // GPU assembled matrix: GMRES matvec on GPU.
                    let matvec_fn = move |x: &ndarray::Array1<Complex64>| -> Result<ndarray::Array1<Complex64>, SolverError> {
                        session.matvec(x.view()).map_err(|e| SolverError::LinAlgError(e.to_string()))
                    };
                    {
                    let preconditioned = |x: &ndarray::Array1<Complex64>| -> Result<ndarray::Array1<Complex64>, SolverError> {
                        Ok(assembly::apply_preconditioner(&precond, &matvec_fn(x)?))
                    };
                    iterative::solve_gmres(&preconditioned, &rhs_p, params.solver_tolerance, params.max_iterations)?
                }
                } else {
                    // CPU assembly fallback (no GPU backend available).
                    let matrix = assembly::assemble_interaction_matrix(
                        dipoles,
                        k,
                        self.use_fcd,
                        self.lattice.as_ref(),
                        params.k_bloch,
                        substrate,
                    );
                    let backend = Arc::clone(&self.backend);
                    let matvec_fn = move |x: &ndarray::Array1<Complex64>| -> Result<ndarray::Array1<Complex64>, SolverError> {
                        backend.matvec(&matrix, x).map_err(|e| SolverError::LinAlgError(e.to_string()))
                    };
                    {
                    let preconditioned = |x: &ndarray::Array1<Complex64>| -> Result<ndarray::Array1<Complex64>, SolverError> {
                        Ok(assembly::apply_preconditioner(&precond, &matvec_fn(x)?))
                    };
                    iterative::solve_gmres(&preconditioned, &rhs_p, params.solver_tolerance, params.max_iterations)?
                }
                }
            } else {
                // Matrix exceeds budget: on-the-fly matvec -- O(N) memory, no assembly.
                let lattice_ref = self.lattice.as_ref();
                let matvec_fn = |x: &ndarray::Array1<Complex64>| -> Result<ndarray::Array1<Complex64>, SolverError> {
                    Ok(assembly::matvec_on_the_fly(
                        dipoles,
                        x,
                        k,
                        self.use_fcd,
                        lattice_ref,
                        params.k_bloch,
                        substrate,
                    ))
                };
                {
                    let preconditioned = |x: &ndarray::Array1<Complex64>| -> Result<ndarray::Array1<Complex64>, SolverError> {
                        Ok(assembly::apply_preconditioner(&precond, &matvec_fn(x)?))
                    };
                    iterative::solve_gmres(&preconditioned, &rhs_p, params.solver_tolerance, params.max_iterations)?
                }
            }
        } else if matrix_bytes <= self.matrix_memory_budget {
            // Lattice, substrate, or FCD active: GPU assembly not supported.
            // Matrix fits in budget: assemble once on CPU and use BLAS matvec.
            let matrix = assembly::assemble_interaction_matrix(
                dipoles,
                k,
                self.use_fcd,
                self.lattice.as_ref(),
                params.k_bloch,
                substrate,
            );
            if let Some(session) = self.backend.create_session(3 * n) {
                // GPU matvec with CPU-assembled matrix.
                session
                    .upload_matrix(&matrix)
                    .map_err(|e| SolverError::LinAlgError(e.to_string()))?;
                let matvec_fn = move |x: &ndarray::Array1<Complex64>| -> Result<ndarray::Array1<Complex64>, SolverError> {
                    session.matvec(x.view()).map_err(|e| SolverError::LinAlgError(e.to_string()))
                };
                {
                    let preconditioned = |x: &ndarray::Array1<Complex64>| -> Result<ndarray::Array1<Complex64>, SolverError> {
                        Ok(assembly::apply_preconditioner(&precond, &matvec_fn(x)?))
                    };
                    iterative::solve_gmres(&preconditioned, &rhs_p, params.solver_tolerance, params.max_iterations)?
                }
            } else {
                // Pure CPU cached matvec.
                let backend = Arc::clone(&self.backend);
                let matvec_fn = move |x: &ndarray::Array1<Complex64>| -> Result<ndarray::Array1<Complex64>, SolverError> {
                    backend.matvec(&matrix, x).map_err(|e| SolverError::LinAlgError(e.to_string()))
                };
                {
                    let preconditioned = |x: &ndarray::Array1<Complex64>| -> Result<ndarray::Array1<Complex64>, SolverError> {
                        Ok(assembly::apply_preconditioner(&precond, &matvec_fn(x)?))
                    };
                    iterative::solve_gmres(&preconditioned, &rhs_p, params.solver_tolerance, params.max_iterations)?
                }
            }
        } else {
            // Matrix exceeds budget: on-the-fly matvec -- O(N) memory, no assembly.
            let lattice_ref = self.lattice.as_ref();
            let matvec_fn = |x: &ndarray::Array1<Complex64>| -> Result<ndarray::Array1<Complex64>, SolverError> {
                Ok(assembly::matvec_on_the_fly(
                    dipoles,
                    x,
                    k,
                    self.use_fcd,
                    lattice_ref,
                    params.k_bloch,
                    substrate,
                ))
            };
            {
                    let preconditioned = |x: &ndarray::Array1<Complex64>| -> Result<ndarray::Array1<Complex64>, SolverError> {
                        Ok(assembly::apply_preconditioner(&precond, &matvec_fn(x)?))
                    };
                    iterative::solve_gmres(&preconditioned, &rhs_p, params.solver_tolerance, params.max_iterations)?
                }
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
            environment_n: params.environment_n,
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
        let k = Self::wavenumber(response.wavelength_nm, response.environment_n);

        crate::fields::compute_near_field_map(dipoles, response, plane, k, &self.incident_field)
    }

    fn compute_far_field(
        &self,
        dipoles: &[Dipole],
        response: &DipoleResponse,
        n_theta: usize,
        n_phi: usize,
    ) -> FarFieldMap {
        let k = Self::wavenumber(response.wavelength_nm, response.environment_n);
        crate::fields::compute_far_field(dipoles, response, k, n_theta, n_phi)
    }

    fn method_name(&self) -> &str {
        "Coupled Dipole Approximation (CDA)"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::solver::sommerfeld::SommerfeldGreens;
    use crate::solver::substrate::SubstrateRuntime;
    use crate::types::SimulationParams;
    use lumina_geometry::lattice::LatticeSpec;

    fn one_dipole() -> Vec<Dipole> {
        vec![Dipole::isotropic([0.0, 0.0, 30.0], Complex64::new(1000.0, 10.0), 3.0)]
    }

    fn sommerfeld_runtime() -> SubstrateRuntime {
        SubstrateRuntime::Sommerfeld {
            z_interface: 0.0,
            greens: Arc::new(SommerfeldGreens::new(600.0, Complex64::new(2.25, 0.0), 1.0)),
        }
    }

    fn params_with(substrate: Option<SubstrateRuntime>) -> SimulationParams {
        SimulationParams { substrate_runtime: substrate, ..Default::default() }
    }

    /// C2: a Sommerfeld substrate combined with a periodic lattice used to be
    /// accepted silently, with the reflected tensor evaluated pair-by-pair and
    /// never summed over lattice images -- an inconsistent mixture of a periodic
    /// direct term and a non-periodic substrate term. It must be refused.
    #[test]
    fn sommerfeld_with_lattice_is_rejected() {
        let mut solver = CdaSolver::new(1000);
        solver.lattice = Some(LatticeSpec::chain([100.0, 0.0, 0.0]));

        let err = solver
            .solve_dipoles(&one_dipole(), 600.0, &params_with(Some(sommerfeld_runtime())))
            .expect_err("Sommerfeld + lattice must be refused, not silently approximated");

        match err {
            SolverError::UnsupportedConfiguration(msg) => {
                assert!(
                    msg.contains("Sommerfeld") && msg.contains("periodic"),
                    "error should name both the substrate and the lattice: {msg}"
                );
            }
            other => panic!("expected UnsupportedConfiguration, got {other:?}"),
        }
    }

    /// The guard must not be over-broad: Sommerfeld alone is a supported case.
    #[test]
    fn sommerfeld_without_lattice_is_accepted() {
        let solver = CdaSolver::new(1000);
        solver
            .solve_dipoles(&one_dipole(), 600.0, &params_with(Some(sommerfeld_runtime())))
            .expect("Sommerfeld without a lattice is supported");
    }

    /// ...nor may it catch the Fresnel path, which *is* Ewald-summed correctly.
    #[test]
    fn fresnel_with_lattice_is_accepted() {
        let mut solver = CdaSolver::new(1000);
        solver.lattice = Some(LatticeSpec::chain([100.0, 0.0, 0.0]));
        let fresnel = SubstrateRuntime::Fresnel {
            z_interface: 0.0,
            delta_eps: Complex64::new(0.2, 0.0),
        };
        solver
            .solve_dipoles(&one_dipole(), 600.0, &params_with(Some(fresnel)))
            .expect("Fresnel + lattice is a supported combination");
    }

    /// D11: the iterative path is preconditioned (M^-1 = diag(alpha), since A_ii = alpha^-1),
    /// so it solves M^-1A*p = M^-1E. That must give the same p as the direct solve --
    /// a misapplied preconditioner would silently corrupt every large run.
    ///
    /// Note the preconditioner only *helps convergence* when alpha varies between
    /// dipoles. For a single material M^-1 = alpha*I is a scalar multiple of the
    /// identity, which leaves the Krylov subspace and GMRES convergence exactly
    /// unchanged -- measured residuals matched to 12 digits. It is kept because it
    /// is free and does help mixed-material scenes (residual after 60 iterations
    /// 0.146 -> 0.086 on an alternating-alpha cluster).
    #[test]
    fn preconditioned_gmres_agrees_with_the_direct_solve() {
        let mut dipoles = Vec::new();
        for i in 0..5 {
            for j in 0..5 {
                for k in 0..5 {
                    let alpha = if (i + j + k) % 2 == 0 {
                        Complex64::new(900.0, 120.0)
                    } else {
                        Complex64::new(40.0, 3.0)
                    };
                    dipoles.push(Dipole::isotropic(
                        [i as f64 * 4.0, j as f64 * 4.0, k as f64 * 4.0],
                        alpha,
                        4.0,
                    ));
                }
            }
        }
        let params = SimulationParams {
            max_iterations: 2000,
            solver_tolerance: 1e-10,
            ..Default::default()
        };

        // n = 125: threshold above it takes the direct path, below it takes GMRES.
        let direct = CdaSolver { iterative_threshold: 1000, use_fcd: false, ..Default::default() }
            .solve_dipoles(&dipoles, 600.0, &params)
            .expect("direct solve");
        let iterative = CdaSolver { iterative_threshold: 10, use_fcd: false, ..Default::default() }
            .solve_dipoles(&dipoles, 600.0, &params)
            .expect("preconditioned GMRES must converge");

        let scale = direct
            .moments
            .iter()
            .fold(0.0_f64, |m, v| m.max(v.norm()));
        assert!(scale > 1e-30, "direct solve produced a zero response");

        for (a, b) in direct.moments.iter().zip(iterative.moments.iter()) {
            assert!(
                (a - b).norm() / scale < 1e-6,
                "preconditioned GMRES disagrees with direct: {a:?} vs {b:?}"
            );
        }
    }


    /// C3, at the level a user would notice: an on-axis chain must actually couple
    /// through the lattice.
    ///
    /// Until the collinear branch was restored, `EwaldGreens::evaluate` returned an
    /// exactly zero tensor for every on-axis pair. For this geometry *every* pair is
    /// collinear, so every off-diagonal block vanished and the periodic system became
    /// completely **decoupled** -- each dipole reduced to an isolated one, p = alpha*E_inc.
    ///
    /// Note comparing against a free-space solve would NOT catch that: a decoupled
    /// system also differs from free space. The assertion has to be against the
    /// decoupled result itself.
    #[test]
    fn on_axis_chain_is_coupled_through_the_lattice() {
        use lumina_geometry::lattice::LatticeSpec;

        let alpha = Complex64::new(900.0, 80.0);
        let dipoles: Vec<Dipole> = (0..5)
            .map(|i| Dipole::isotropic([i as f64 * 12.0, 0.0, 0.0], alpha, 4.0))
            .collect();
        let params = SimulationParams { environment_n: 1.0, ..Default::default() };
        let wavelength = 600.0;
        let k = 2.0 * std::f64::consts::PI * params.environment_n / wavelength;

        let mut solver = CdaSolver { use_fcd: false, ..Default::default() };
        solver.lattice = Some(LatticeSpec::chain([80.0, 0.0, 0.0]));
        let periodic = solver
            .solve_dipoles(&dipoles, wavelength, &params)
            .expect("periodic solve");

        // The decoupled answer: with every off-diagonal block zero, A = diag(alpha^-1)
        // and the solve collapses to p_i = alpha*E_inc(r_i).
        let mut max_dev = 0.0_f64;
        let mut scale = 0.0_f64;
        for (i, d) in dipoles.iter().enumerate() {
            let e_inc = solver.incident_field.at_position(&d.position, k);
            for c in 0..3 {
                let isolated = alpha * e_inc[c];
                scale = scale.max(isolated.norm());
                max_dev = max_dev.max((periodic.moments[[i, c]] - isolated).norm());
            }
        }
        assert!(scale > 1e-30, "incident field is zero -- fixture is broken");

        assert!(
            max_dev / scale > 1e-3,
            "periodic solve equals the DECOUPLED result (rel dev {:.2e}) -- every on-axis pair's Green's tensor is zero, so the lattice does nothing",
            max_dev / scale
        );
    }

}
