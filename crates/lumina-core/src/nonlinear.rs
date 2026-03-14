//! Nonlinear optical response (SHG/THG).
//!
//! Implements the second- and third-order nonlinear optical response within
//! the CDA framework.
//!
//! For SHG, the nonlinear dipole moment at $2\omega$ is driven by the local
//! field at $\omega$:
//!
//! $$\mathbf{p}_i(2\omega) = \epsilon_0 \boldsymbol{\chi}^{(2)} :
//!   \mathbf{E}_{\text{loc}}(\omega)\, \mathbf{E}_{\text{loc}}(\omega)$$
//!
//! For THG, the nonlinear dipole moment at $3\omega$ is driven by the cube of
//! the local field:
//!
//! $$p_{i,a}(3\omega) = \sum_{b,c,d} \chi^{(3)}_{abcd}\,
//!   E_{\text{loc},b}(\omega)\, E_{\text{loc},c}(\omega)\, E_{\text{loc},d}(\omega)$$
//!
//! In both cases, the sources are re-radiated self-consistently through the CDA
//! interaction matrix at the harmonic frequency, with **no** incident field at
//! the harmonic.
//!
//! # SHG Workflow
//!
//! 1. Solve the linear CDA system at $\omega$ — obtain [`DipoleResponse`]
//!    (which stores `local_fields`).
//! 2. Assign a [`Chi2Tensor`] to each dipole. For a uniform structure use
//!    `vec![chi2; n_dipoles]`.
//! 3. Build dipoles with polarisability evaluated at $2\omega = \lambda/2$.
//! 4. Call [`compute_shg_response`].
//!
//! # THG Workflow
//!
//! 1. Solve the linear CDA system at $\omega$ — obtain [`DipoleResponse`].
//! 2. Assign a [`Chi3Tensor`] to each dipole.
//! 3. Build dipoles with polarisability evaluated at $3\omega = \lambda/3$.
//! 4. Call [`compute_thg_response`].

use std::sync::Arc;

use ndarray::{Array1, Array2};
use num_complex::Complex64;

use crate::solver::cda::{assembly, direct, iterative, CdaSolver};
use crate::solver::SolverError;
use crate::types::{
    Chi2Tensor, Chi3Tensor, Dipole, DipoleResponse, FarFieldMap, ShgResult, SimulationParams,
    ThgResult,
};

/// Compute nonlinear dipole source terms at $2\omega$ from local fields at $\omega$.
///
/// For each dipole $i$:
/// $$p^{\text{NL}}_{i,a}(2\omega) = \sum_{b,c} \chi^{(2)}_{abc}\,
///   E_{\text{loc},b,i}(\omega)\, E_{\text{loc},c,i}(\omega)$$
///
/// # Arguments
/// * `local_fields`  — Local electric fields at $\omega$, shape $(N, 3)$.
/// * `chi2_tensors`  — One [`Chi2Tensor`] per dipole; length must equal $N$.
///
/// # Returns
/// Source moment array of shape $(N, 3)$.
pub fn compute_shg_sources(
    local_fields: &Array2<Complex64>,
    chi2_tensors: &[Chi2Tensor],
) -> Array2<Complex64> {
    let n = local_fields.nrows();
    assert_eq!(
        chi2_tensors.len(),
        n,
        "chi2_tensors must have one entry per dipole"
    );

    let mut sources = Array2::zeros((n, 3));
    for i in 0..n {
        let e: [Complex64; 3] = [
            local_fields[[i, 0]],
            local_fields[[i, 1]],
            local_fields[[i, 2]],
        ];
        let p_nl = chi2_tensors[i].contract(&e);
        for a in 0..3 {
            sources[[i, a]] = p_nl[a];
        }
    }
    sources
}

/// Compute the self-consistent SHG response.
///
/// Drives the CDA interaction matrix at $2\omega$ with nonlinear source
/// moments computed from the linear local fields, then solves for the
/// re-radiated dipole moments at the harmonic frequency.
///
/// # Arguments
/// * `solver`         — CDA solver instance (for configuration and backend).
/// * `omega_response` — Linear dipole response at $\omega$ (contains local fields).
/// * `dipoles_2omega` — Dipoles with polarisability evaluated at $2\omega$.
///                      Positions must match those used for `omega_response`.
/// * `chi2_tensors`   — One [`Chi2Tensor`] per dipole; length $N$.
/// * `params_2omega`  — Simulation parameters for the $2\omega$ solve.
///                      The harmonic wavelength is taken from
///                      `omega_response.wavelength_nm / 2`.
/// * `calc_far_field` — If `true`, compute the far-field pattern at $2\omega$.
///
/// # Returns
/// [`ShgResult`] with source moments, driven moments, SHG intensity, and
/// optionally the far-field pattern.
pub fn compute_shg_response(
    solver: &CdaSolver,
    omega_response: &DipoleResponse,
    dipoles_2omega: &[Dipole],
    chi2_tensors: &[Chi2Tensor],
    params_2omega: &SimulationParams,
    calc_far_field: bool,
) -> Result<ShgResult, SolverError> {
    let n = dipoles_2omega.len();
    if n == 0 {
        return Err(SolverError::InvalidGeometry("No dipoles provided".into()));
    }
    assert_eq!(
        omega_response.moments.nrows(),
        n,
        "dipoles_2omega length must match omega_response dipole count"
    );

    let fundamental_nm = omega_response.wavelength_nm;
    let harmonic_nm = fundamental_nm / 2.0;
    let k_2omega = 2.0 * std::f64::consts::PI * params_2omega.environment_n / harmonic_nm;

    // Step 1: Nonlinear source moments p_NL from local fields at ω.
    let source_moments = compute_shg_sources(&omega_response.local_fields, chi2_tensors);

    // Substrate runtime for the 2ω solve.
    let substrate = params_2omega.substrate_runtime.as_ref();

    // Step 2: Interaction matrix at 2ω.
    let matrix = assembly::assemble_interaction_matrix(
        dipoles_2omega,
        k_2omega,
        solver.use_fcd,
        solver.cell_size,
        solver.lattice.as_ref(),
        params_2omega.k_bloch,
        substrate,
    );

    // Step 3: RHS = flat 3N vector of nonlinear sources (no incident field at 2ω).
    let mut rhs = Array1::zeros(3 * n);
    for i in 0..n {
        for c in 0..3 {
            rhs[3 * i + c] = source_moments[[i, c]];
        }
    }

    // Step 4: Solve A(2ω) · p(2ω) = p_NL.
    let solution = if n <= solver.iterative_threshold {
        direct::solve_direct(&matrix, &rhs)
            .map_err(|e| SolverError::LinAlgError(e.to_string()))?
    } else if let Some(session) = solver.backend.create_session(3 * n) {
        session
            .upload_matrix(&matrix)
            .map_err(|e| SolverError::LinAlgError(e.to_string()))?;
        let matvec_fn =
            move |x: &Array1<Complex64>| -> Result<Array1<Complex64>, SolverError> {
                session
                    .matvec(x.view())
                    .map_err(|e| SolverError::LinAlgError(e.to_string()))
            };
        iterative::solve_gmres(
            &matvec_fn,
            &rhs,
            params_2omega.solver_tolerance,
            params_2omega.max_iterations,
        )?
    } else {
        let backend = Arc::clone(&solver.backend);
        let matvec_fn =
            move |x: &Array1<Complex64>| -> Result<Array1<Complex64>, SolverError> {
                backend
                    .matvec(&matrix, x)
                    .map_err(|e| SolverError::LinAlgError(e.to_string()))
            };
        iterative::solve_gmres(
            &matvec_fn,
            &rhs,
            params_2omega.solver_tolerance,
            params_2omega.max_iterations,
        )?
    };

    // Step 5: Reshape into (N, 3) driven moments.
    let mut driven_moments = Array2::zeros((n, 3));
    for i in 0..n {
        for c in 0..3 {
            driven_moments[[i, c]] = solution[3 * i + c];
        }
    }

    // Step 6: Total SHG intensity ∝ Σᵢ |p_i(2ω)|².
    let mut shg_intensity = 0.0_f64;
    for i in 0..n {
        for c in 0..3 {
            shg_intensity += driven_moments[[i, c]].norm_sqr();
        }
    }

    // Step 7: Optional far-field pattern at 2ω.
    let far_field: Option<FarFieldMap> = if calc_far_field {
        let response_2omega = DipoleResponse {
            wavelength_nm: harmonic_nm,
            moments: driven_moments.clone(),
            local_fields: Array2::zeros((n, 3)),
        };
        Some(crate::fields::compute_far_field(
            dipoles_2omega,
            &response_2omega,
            k_2omega,
            36,
            72,
        ))
    } else {
        None
    };

    Ok(ShgResult {
        fundamental_nm,
        harmonic_nm,
        source_moments,
        driven_moments,
        shg_intensity,
        far_field,
    })
}

/// Compute nonlinear dipole source terms at $3\omega$ from local fields at $\omega$.
///
/// For each dipole $i$:
/// $$p^{\text{NL}}_{i,a}(3\omega) = \sum_{b,c,d} \chi^{(3)}_{abcd}\,
///   E_{\text{loc},b,i}(\omega)\, E_{\text{loc},c,i}(\omega)\, E_{\text{loc},d,i}(\omega)$$
///
/// # Arguments
/// * `local_fields`  — Local electric fields at $\omega$, shape $(N, 3)$.
/// * `chi3_tensors`  — One [`Chi3Tensor`] per dipole; length must equal $N$.
///
/// # Returns
/// Source moment array of shape $(N, 3)$.
pub fn compute_thg_sources(
    local_fields: &Array2<Complex64>,
    chi3_tensors: &[Chi3Tensor],
) -> Array2<Complex64> {
    let n = local_fields.nrows();
    assert_eq!(
        chi3_tensors.len(),
        n,
        "chi3_tensors must have one entry per dipole"
    );

    let mut sources = Array2::zeros((n, 3));
    for i in 0..n {
        let e: [Complex64; 3] = [
            local_fields[[i, 0]],
            local_fields[[i, 1]],
            local_fields[[i, 2]],
        ];
        let p_nl = chi3_tensors[i].contract(&e);
        for a in 0..3 {
            sources[[i, a]] = p_nl[a];
        }
    }
    sources
}

/// Compute the self-consistent THG response.
///
/// Drives the CDA interaction matrix at $3\omega$ with nonlinear source
/// moments computed from the linear local fields, then solves for the
/// re-radiated dipole moments at the harmonic frequency.
///
/// # Arguments
/// * `solver`         — CDA solver instance (for configuration and backend).
/// * `omega_response` — Linear dipole response at $\omega$ (contains local fields).
/// * `dipoles_3omega` — Dipoles with polarisability evaluated at $3\omega$.
///                      Positions must match those used for `omega_response`.
/// * `chi3_tensors`   — One [`Chi3Tensor`] per dipole; length $N$.
/// * `params_3omega`  — Simulation parameters for the $3\omega$ solve.
///                      The harmonic wavelength is taken from
///                      `omega_response.wavelength_nm / 3`.
/// * `calc_far_field` — If `true`, compute the far-field pattern at $3\omega$.
///
/// # Returns
/// [`ThgResult`] with source moments, driven moments, THG intensity, and
/// optionally the far-field pattern.
pub fn compute_thg_response(
    solver: &CdaSolver,
    omega_response: &DipoleResponse,
    dipoles_3omega: &[Dipole],
    chi3_tensors: &[Chi3Tensor],
    params_3omega: &SimulationParams,
    calc_far_field: bool,
) -> Result<ThgResult, SolverError> {
    let n = dipoles_3omega.len();
    if n == 0 {
        return Err(SolverError::InvalidGeometry("No dipoles provided".into()));
    }
    assert_eq!(
        omega_response.moments.nrows(),
        n,
        "dipoles_3omega length must match omega_response dipole count"
    );

    let fundamental_nm = omega_response.wavelength_nm;
    let harmonic_nm = fundamental_nm / 3.0;
    let k_3omega = 2.0 * std::f64::consts::PI * params_3omega.environment_n / harmonic_nm;

    // Step 1: Nonlinear source moments p_NL from local fields at ω.
    let source_moments = compute_thg_sources(&omega_response.local_fields, chi3_tensors);

    // Substrate runtime for the 3ω solve.
    let substrate = params_3omega.substrate_runtime.as_ref();

    // Step 2: Interaction matrix at 3ω.
    let matrix = assembly::assemble_interaction_matrix(
        dipoles_3omega,
        k_3omega,
        solver.use_fcd,
        solver.cell_size,
        solver.lattice.as_ref(),
        params_3omega.k_bloch,
        substrate,
    );

    // Step 3: RHS = flat 3N vector of nonlinear sources (no incident field at 3ω).
    let mut rhs = Array1::zeros(3 * n);
    for i in 0..n {
        for c in 0..3 {
            rhs[3 * i + c] = source_moments[[i, c]];
        }
    }

    // Step 4: Solve A(3ω) · p(3ω) = p_NL.
    let solution = if n <= solver.iterative_threshold {
        direct::solve_direct(&matrix, &rhs)
            .map_err(|e| SolverError::LinAlgError(e.to_string()))?
    } else if let Some(session) = solver.backend.create_session(3 * n) {
        session
            .upload_matrix(&matrix)
            .map_err(|e| SolverError::LinAlgError(e.to_string()))?;
        let matvec_fn =
            move |x: &Array1<Complex64>| -> Result<Array1<Complex64>, SolverError> {
                session
                    .matvec(x.view())
                    .map_err(|e| SolverError::LinAlgError(e.to_string()))
            };
        iterative::solve_gmres(
            &matvec_fn,
            &rhs,
            params_3omega.solver_tolerance,
            params_3omega.max_iterations,
        )?
    } else {
        let backend = Arc::clone(&solver.backend);
        let matvec_fn =
            move |x: &Array1<Complex64>| -> Result<Array1<Complex64>, SolverError> {
                backend
                    .matvec(&matrix, x)
                    .map_err(|e| SolverError::LinAlgError(e.to_string()))
            };
        iterative::solve_gmres(
            &matvec_fn,
            &rhs,
            params_3omega.solver_tolerance,
            params_3omega.max_iterations,
        )?
    };

    // Step 5: Reshape into (N, 3) driven moments.
    let mut driven_moments = Array2::zeros((n, 3));
    for i in 0..n {
        for c in 0..3 {
            driven_moments[[i, c]] = solution[3 * i + c];
        }
    }

    // Step 6: Total THG intensity ∝ Σᵢ |p_i(3ω)|².
    let mut thg_intensity = 0.0_f64;
    for i in 0..n {
        for c in 0..3 {
            thg_intensity += driven_moments[[i, c]].norm_sqr();
        }
    }

    // Step 7: Optional far-field pattern at 3ω.
    let far_field: Option<FarFieldMap> = if calc_far_field {
        let response_3omega = DipoleResponse {
            wavelength_nm: harmonic_nm,
            moments: driven_moments.clone(),
            local_fields: Array2::zeros((n, 3)),
        };
        Some(crate::fields::compute_far_field(
            dipoles_3omega,
            &response_3omega,
            k_3omega,
            36,
            72,
        ))
    } else {
        None
    };

    Ok(ThgResult {
        fundamental_nm,
        harmonic_nm,
        source_moments,
        driven_moments,
        thg_intensity,
        far_field,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;
    use num_complex::Complex64;
    use crate::types::{Chi2Tensor, Chi3Tensor, Dipole, DipoleResponse, SimulationParams};
    use crate::solver::cda::CdaSolver;

    fn make_response(wavelength_nm: f64, fields: Vec<[f64; 3]>) -> DipoleResponse {
        let n = fields.len();
        let mut local_fields = Array2::zeros((n, 3));
        for (i, f) in fields.iter().enumerate() {
            for c in 0..3 {
                local_fields[[i, c]] = Complex64::from(f[c]);
            }
        }
        DipoleResponse {
            wavelength_nm,
            moments: Array2::zeros((n, 3)),
            local_fields,
        }
    }

    /// Zero χ^(2) → zero source moments, regardless of local field.
    #[test]
    fn test_zero_chi2_gives_zero_sources() {
        let response = make_response(600.0, vec![[1.0, 0.5, 0.2]; 4]);
        let chi2 = vec![Chi2Tensor::zero(); 4];
        let sources = compute_shg_sources(&response.local_fields, &chi2);
        for i in 0..4 {
            for c in 0..3 {
                assert!(
                    sources[[i, c]].norm() < 1e-14,
                    "Expected zero source for zero chi2 at dipole {i}, component {c}"
                );
            }
        }
    }

    /// Single dipole: χ_xxx = 1, E = [E₀, 0, 0].
    /// Expected: p_x = E₀², p_y = p_z = 0.
    #[test]
    fn test_shg_source_xxx_contraction() {
        let e0 = 2.5_f64;
        let response = make_response(600.0, vec![[e0, 0.0, 0.0]]);
        let mut chi2 = Chi2Tensor::zero();
        chi2.components[0 * 9 + 0 * 3 + 0] = Complex64::new(1.0, 0.0); // χ_xxx = 1
        let sources = compute_shg_sources(&response.local_fields, &[chi2]);

        let expected = Complex64::new(e0 * e0, 0.0);
        assert!(
            (sources[[0, 0]] - expected).norm() < 1e-12,
            "p_x should equal chi_xxx * E₀²"
        );
        assert!(sources[[0, 1]].norm() < 1e-12, "p_y should be zero");
        assert!(sources[[0, 2]].norm() < 1e-12, "p_z should be zero");
    }

    /// Isotropic surface tensor: verify non-zero element positions and count.
    #[test]
    fn test_isotropic_surface_tensor_structure() {
        let chi_zzz = Complex64::new(1.0, 0.0);
        let chi_zxx = Complex64::new(0.5, 0.0);
        let t = Chi2Tensor::isotropic_surface(chi_zzz, chi_zxx);

        // χ_zzz at (2,2,2) → index 26
        assert!((t.components[26] - chi_zzz).norm() < 1e-14, "χ_zzz");
        // χ_zxx at (2,0,0) → index 18
        assert!((t.components[18] - chi_zxx).norm() < 1e-14, "χ_zxx");
        // χ_zyy at (2,1,1) → index 22
        assert!((t.components[22] - chi_zxx).norm() < 1e-14, "χ_zyy");
        // χ_xxz at (0,0,2) → index 2
        assert!((t.components[2] - chi_zxx).norm() < 1e-14, "χ_xxz");
        // χ_xzx at (0,2,0) → index 6
        assert!((t.components[6] - chi_zxx).norm() < 1e-14, "χ_xzx");
        // χ_yyz at (1,1,2) → index 14
        assert!((t.components[14] - chi_zxx).norm() < 1e-14, "χ_yyz");
        // χ_yzy at (1,2,1) → index 16
        assert!((t.components[16] - chi_zxx).norm() < 1e-14, "χ_yzy");

        let nonzero = t.components.iter().filter(|c| c.norm() > 1e-14).count();
        assert_eq!(nonzero, 7, "C_∞v tensor has exactly 7 non-zero elements");
    }

    /// End-to-end: zero χ^(2) → zero SHG intensity from the full solver.
    #[test]
    fn test_shg_response_zero_chi2() {
        let solver = CdaSolver {
            use_fcd: false,
            ..CdaSolver::default()
        };

        // Single isotropic dipole at origin with unit polarisability at 2ω.
        let dipole_2omega =
            Dipole::isotropic([0.0, 0.0, 0.0], Complex64::new(1.0, 0.0));

        // Linear response at 600 nm with a non-zero local field.
        let omega_response = make_response(600.0, vec![[3.0, 1.0, 0.5]]);

        let chi2 = vec![Chi2Tensor::zero()];
        let params = SimulationParams::default();

        let result = compute_shg_response(
            &solver,
            &omega_response,
            &[dipole_2omega],
            &chi2,
            &params,
            false,
        )
        .unwrap();

        assert!(
            result.shg_intensity < 1e-20,
            "Zero chi2 must give zero SHG intensity, got {}",
            result.shg_intensity
        );
        assert_eq!(result.fundamental_nm, 600.0);
        assert_eq!(result.harmonic_nm, 300.0);
    }

    /// Non-zero χ^(2): SHG intensity scales as |E_loc|⁴.
    #[test]
    fn test_shg_intensity_scales_with_field() {
        let solver = CdaSolver {
            use_fcd: false,
            ..CdaSolver::default()
        };
        let dipole = Dipole::isotropic([0.0, 0.0, 0.0], Complex64::new(1.0, 0.0));
        let mut chi2 = Chi2Tensor::zero();
        chi2.components[0 * 9 + 0 * 3 + 0] = Complex64::new(1.0, 0.0); // χ_xxx = 1
        let params = SimulationParams::default();

        // Field amplitude 1
        let r1 = compute_shg_response(
            &solver,
            &make_response(600.0, vec![[1.0, 0.0, 0.0]]),
            &[dipole.clone()],
            &[chi2.clone()],
            &params,
            false,
        )
        .unwrap();

        // Field amplitude 2 — SHG intensity should be 2⁴ = 16× larger.
        let r2 = compute_shg_response(
            &solver,
            &make_response(600.0, vec![[2.0, 0.0, 0.0]]),
            &[dipole],
            &[chi2],
            &params,
            false,
        )
        .unwrap();

        let ratio = r2.shg_intensity / r1.shg_intensity;
        assert!(
            (ratio - 16.0).abs() < 1e-8,
            "SHG intensity should scale as E⁴; expected ratio 16, got {ratio:.6}"
        );
    }

    // ── THG tests ─────────────────────────────────────────────────────────────

    /// Zero χ^(3) → zero source moments, regardless of local field.
    #[test]
    fn test_zero_chi3_gives_zero_sources() {
        let response = make_response(600.0, vec![[1.0, 0.5, 0.2]; 4]);
        let chi3 = vec![Chi3Tensor::zero(); 4];
        let sources = compute_thg_sources(&response.local_fields, &chi3);
        for i in 0..4 {
            for c in 0..3 {
                assert!(
                    sources[[i, c]].norm() < 1e-14,
                    "Expected zero source for zero chi3 at dipole {i}, component {c}"
                );
            }
        }
    }

    /// Single dipole: χ_xxxx = 1, E = [E₀, 0, 0].
    /// Expected: p_x = E₀³, p_y = p_z = 0.
    #[test]
    fn test_thg_source_xxxx_contraction() {
        let e0 = 2.0_f64;
        let response = make_response(600.0, vec![[e0, 0.0, 0.0]]);
        let mut chi3 = Chi3Tensor::zero();
        chi3.components[0 * 27 + 0 * 9 + 0 * 3 + 0] = Complex64::new(1.0, 0.0); // χ_xxxx = 1
        let sources = compute_thg_sources(&response.local_fields, &[chi3]);

        let expected = Complex64::new(e0 * e0 * e0, 0.0);
        assert!(
            (sources[[0, 0]] - expected).norm() < 1e-12,
            "p_x should equal chi_xxxx * E₀³, got {}",
            sources[[0, 0]]
        );
        assert!(sources[[0, 1]].norm() < 1e-12, "p_y should be zero");
        assert!(sources[[0, 2]].norm() < 1e-12, "p_z should be zero");
    }

    /// Isotropic bulk tensor: verify known non-zero element positions.
    #[test]
    fn test_isotropic_bulk_tensor_structure() {
        let chi_xxxx = Complex64::new(1.0, 0.0);
        let chi_xxyy = Complex64::new(0.25, 0.0);
        let t = Chi3Tensor::isotropic_bulk(chi_xxxx, chi_xxyy);

        // χ_xxxx at (0,0,0,0) → index 0
        assert!(
            (t.components[0] - chi_xxxx).norm() < 1e-14,
            "χ_xxxx should equal chi_xxxx"
        );
        // χ_xxyy at (0,0,1,1) → index 0*27+0*9+1*3+1 = 4
        assert!(
            (t.components[4] - chi_xxyy).norm() < 1e-14,
            "χ_xxyy should equal chi_xxyy"
        );
        // χ_xyxy at (0,1,0,1) → index 0*27+1*9+0*3+1 = 10
        assert!(
            (t.components[10] - chi_xxyy).norm() < 1e-14,
            "χ_xyxy should equal chi_xxyy"
        );
        // χ_xyyx at (0,1,1,0) → index 0*27+1*9+1*3+0 = 12
        assert!(
            (t.components[12] - chi_xxyy).norm() < 1e-14,
            "χ_xyyx should equal chi_xxyy"
        );
    }

    /// End-to-end: zero χ^(3) → zero THG intensity from the full solver.
    #[test]
    fn test_thg_response_zero_chi3() {
        let solver = CdaSolver {
            use_fcd: false,
            ..CdaSolver::default()
        };

        // Single isotropic dipole at origin with unit polarisability at 3ω.
        let dipole_3omega =
            Dipole::isotropic([0.0, 0.0, 0.0], Complex64::new(1.0, 0.0));

        // Linear response at 600 nm with a non-zero local field.
        let omega_response = make_response(600.0, vec![[3.0, 1.0, 0.5]]);

        let chi3 = vec![Chi3Tensor::zero()];
        let params = SimulationParams::default();

        let result = compute_thg_response(
            &solver,
            &omega_response,
            &[dipole_3omega],
            &chi3,
            &params,
            false,
        )
        .unwrap();

        assert!(
            result.thg_intensity < 1e-20,
            "Zero chi3 must give zero THG intensity, got {}",
            result.thg_intensity
        );
        assert_eq!(result.fundamental_nm, 600.0);
        assert!((result.harmonic_nm - 200.0).abs() < 1e-10);
    }

    /// Non-zero χ^(3): THG intensity scales as |E_loc|⁶.
    #[test]
    fn test_thg_intensity_scales_with_field() {
        let solver = CdaSolver {
            use_fcd: false,
            ..CdaSolver::default()
        };
        let dipole = Dipole::isotropic([0.0, 0.0, 0.0], Complex64::new(1.0, 0.0));
        let mut chi3 = Chi3Tensor::zero();
        chi3.components[0 * 27 + 0 * 9 + 0 * 3 + 0] = Complex64::new(1.0, 0.0); // χ_xxxx = 1
        let params = SimulationParams::default();

        // Field amplitude 1
        let r1 = compute_thg_response(
            &solver,
            &make_response(600.0, vec![[1.0, 0.0, 0.0]]),
            &[dipole.clone()],
            &[chi3.clone()],
            &params,
            false,
        )
        .unwrap();

        // Field amplitude 2 — THG intensity should be 2⁶ = 64× larger.
        let r2 = compute_thg_response(
            &solver,
            &make_response(600.0, vec![[2.0, 0.0, 0.0]]),
            &[dipole],
            &[chi3],
            &params,
            false,
        )
        .unwrap();

        let ratio = r2.thg_intensity / r1.thg_intensity;
        assert!(
            (ratio - 64.0).abs() < 1e-6,
            "THG intensity should scale as E⁶; expected ratio 64, got {ratio:.6}"
        );
    }
}
