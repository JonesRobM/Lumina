//! Simulation runner: ties together geometry, materials, and solver.

use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use anyhow::Result;
use num_complex::Complex64;
use rayon::prelude::*;

use lumina_compute::ComputeBackend;
use lumina_core::nonlinear::{compute_shg_response, compute_thg_response};
use lumina_core::solver::cda::CdaSolver;
use lumina_core::solver::substrate::substrate_reflection_factor;
use lumina_core::solver::{NearFieldPlane, OpticalSolver, SolverError};
use lumina_core::types::{
    clausius_mossotti, ellipsoid_polarisability_tensor, radiative_correction, Chi2Tensor,
    Chi3Tensor, CrossSections, Dipole, NearFieldMap, SimulationParams,
};
use lumina_materials::johnson_christy::JohnsonChristyMaterial;
use lumina_materials::palik::PalikMaterial;
use lumina_materials::provider::MaterialProvider;

use crate::config::{JobConfig, NonlinearConfig, WavelengthSpec};

/// SHG intensity at a single fundamental wavelength.
#[derive(Debug, Clone)]
pub struct ShgPoint {
    /// Fundamental wavelength (nm).
    pub fundamental_nm: f64,
    /// Second-harmonic wavelength λ/2 (nm).
    pub harmonic_nm: f64,
    /// Total SHG intensity ∝ Σᵢ |p_i(2ω)|² (nm⁶, unnormalised).
    pub shg_intensity: f64,
}

/// THG intensity at a single fundamental wavelength.
#[derive(Debug, Clone)]
pub struct ThgPoint {
    /// Fundamental wavelength (nm).
    pub fundamental_nm: f64,
    /// Third-harmonic wavelength λ/3 (nm).
    pub harmonic_nm: f64,
    /// Total THG intensity ∝ Σᵢ |p_i(3ω)|² (nm⁹, unnormalised).
    pub thg_intensity: f64,
}

/// Results from a simulation run.
pub struct SimulationOutput {
    pub spectra: Vec<CrossSections>,
    pub near_field: Option<NearFieldMap>,
    /// SHG spectrum, one point per fundamental wavelength. Empty if SHG was
    /// not requested or every wavelength failed material resolution at 2ω.
    pub shg_spectra: Vec<ShgPoint>,
    /// THG spectrum, one point per fundamental wavelength. Empty if THG was
    /// not requested or every wavelength failed material resolution at λ/3.
    pub thg_spectra: Vec<ThgPoint>,
}

/// Run a full simulation from a parsed job configuration.
pub fn run_simulation(job: &JobConfig) -> Result<SimulationOutput> {
    // Build wavelength grid
    let wavelengths = match &job.simulation.wavelengths {
        WavelengthSpec::Range { range, points } => {
            let start = range[0];
            let end = range[1];
            (0..*points)
                .map(|i| start + (end - start) * i as f64 / (*points - 1).max(1) as f64)
                .collect::<Vec<_>>()
        }
        WavelengthSpec::List { values } => values.clone(),
    };

    let params = SimulationParams {
        wavelength_range: [
            *wavelengths.first().unwrap_or(&400.0),
            *wavelengths.last().unwrap_or(&900.0),
        ],
        num_wavelengths: wavelengths.len(),
        environment_n: job.simulation.environment_n,
        solver_tolerance: job.simulation.solver_tolerance,
        max_iterations: job.simulation.max_iterations,
        k_bloch: [0.0, 0.0, 0.0],
        substrate_z_interface: job.simulation.substrate
            .as_ref()
            .map(|s| s.z_interface)
            .unwrap_or(0.0),
        substrate_delta_eps: None, // resolved per-wavelength below
    };

    // Build dipoles from the scene specification.
    let geom = job.geometry.build_geometry()
        .map_err(|e| anyhow::anyhow!("Geometry error: {}", e))?;

    for obj in &job.geometry.objects {
        println!(
            "  Object '{}': material={}, spacing={} nm",
            obj.name, obj.material, obj.dipole_spacing
        );
    }
    let total_dipoles = geom.positions.len();
    println!("Total dipoles: {}", total_dipoles);
    if total_dipoles == 0 {
        anyhow::bail!("No dipoles generated — check geometry configuration");
    }

    // Load all material providers
    let gold   = JohnsonChristyMaterial::gold();
    let silver = JohnsonChristyMaterial::silver();
    let copper = JohnsonChristyMaterial::copper();
    let tio2   = PalikMaterial::tio2();
    let sio2   = PalikMaterial::sio2();

    // Select compute backend based on config.
    let backend: Arc<dyn ComputeBackend> = create_backend(&job.simulation.backend);

    // Run simulation across wavelengths — parallel via Rayon.
    let solver = CdaSolver {
        backend,
        substrate: job.simulation.substrate.clone(),
        ..Default::default()
    };

    let n_wl = wavelengths.len();
    let progress = AtomicUsize::new(0);

    println!("Solving {} wavelengths in parallel (rayon)...", n_wl);

    let all_spectra: Vec<CrossSections> = wavelengths
        .par_iter()
        .map(|&wl| {
            let k = 2.0 * std::f64::consts::PI * params.environment_n / wl;
            let epsilon_m = params.environment_n * params.environment_n;

            // Resolve substrate Fresnel factor for this wavelength.
            let mut params_wl = params.clone();
            if let Some(sub) = &job.simulation.substrate {
                let eps_sub = resolve_material_epsilon(
                    &sub.material, wl, &gold, &silver, &copper, &tio2, &sio2,
                ).map_err(|e| anyhow::anyhow!("{}", e))?;
                params_wl.substrate_delta_eps = Some(substrate_reflection_factor(eps_sub, epsilon_m, sub.use_retarded));
            }

            let mut dipoles = Vec::with_capacity(total_dipoles);
            for i in 0..total_dipoles {
                let epsilon = resolve_material_epsilon(
                    &geom.materials[i], wl, &gold, &silver, &copper, &tio2, &sio2,
                )
                .map_err(|e| anyhow::anyhow!("{}", e))?;
                let volume = geom.spacings[i].powi(3);
                let dipole = if let Some(hint) = &geom.hints[i] {
                    // Anisotropic path: ellipsoidal depolarisation factors.
                    let tensor = ellipsoid_polarisability_tensor(
                        volume, hint.depol_factors, epsilon, epsilon_m, k, hint.rotation,
                    );
                    Dipole { position: geom.positions[i], polarisability: tensor }
                } else {
                    // Isotropic path: standard spherical CM.
                    let alpha_cm = clausius_mossotti(volume, epsilon, epsilon_m);
                    let alpha = radiative_correction(alpha_cm, k);
                    Dipole::isotropic(geom.positions[i], alpha)
                };
                dipoles.push(dipole);
            }

            let cs = match solver.compute_cross_sections(&dipoles, wl, &params_wl) {
                Ok(cs) => cs,
                Err(SolverError::ConvergenceFailure { max_iter, residual }) => {
                    eprintln!(
                        "  WARNING: λ={:.1} nm did not converge after {} iterations \
                         (residual: {:.2e}). Consider increasing max_iterations, \
                         refining dipole spacing, or checking the material at this wavelength.",
                        wl, max_iter, residual
                    );
                    CrossSections {
                        wavelength_nm: wl,
                        extinction: f64::NAN,
                        absorption: f64::NAN,
                        scattering: f64::NAN,
                        circular_dichroism: None,
                    }
                }
                Err(e) => return Err(anyhow::anyhow!("Solver error at λ={:.1} nm: {}", wl, e)),
            };

            let done = progress.fetch_add(1, Ordering::Relaxed) + 1;
            if done.is_multiple_of(10) || done == 1 || done == n_wl {
                println!(
                    "  [{}/{}] λ={:.1} nm: C_ext={:.2e}, C_abs={:.2e}, C_sca={:.2e}",
                    done, n_wl, wl, cs.extinction, cs.absorption, cs.scattering
                );
            }

            Ok(cs)
        })
        .collect::<Result<Vec<_>>>()?;

    // Report convergence failures
    let n_failed = all_spectra.iter().filter(|cs| cs.extinction.is_nan()).count();
    if n_failed > 0 {
        eprintln!(
            "WARNING: {}/{} wavelengths did not converge. \
             These points appear as NaN in the output. \
             Try reducing dipole spacing or increasing solver iterations.",
            n_failed, n_wl
        );
    }

    // Find peak extinction wavelength for near-field computation
    let (peak_wl, peak_ext) = all_spectra
        .iter()
        .map(|cs| (cs.wavelength_nm, cs.extinction))
        .fold((wavelengths[0], 0.0_f64), |(pw, pe), (w, e)| {
            if e > pe { (w, e) } else { (pw, pe) }
        });

    // Optionally compute near-field map at peak wavelength
    let near_field = if job.output.save_near_field && peak_ext > 0.0 {
        println!("Computing near-field map at λ={:.1} nm (peak extinction)...", peak_wl);

        // Rebuild dipoles at the peak wavelength
        let k_peak = 2.0 * std::f64::consts::PI * params.environment_n / peak_wl;
        let epsilon_m = params.environment_n * params.environment_n;
        let mut peak_params = params.clone();
        if let Some(sub) = &job.simulation.substrate {
            let eps_sub = resolve_material_epsilon(
                &sub.material, peak_wl, &gold, &silver, &copper, &tio2, &sio2,
            )?;
            peak_params.substrate_delta_eps = Some(substrate_reflection_factor(eps_sub, epsilon_m, sub.use_retarded));
        }
        let mut peak_dipoles = Vec::with_capacity(total_dipoles);
        for i in 0..total_dipoles {
            let epsilon = resolve_material_epsilon(
                &geom.materials[i], peak_wl, &gold, &silver, &copper, &tio2, &sio2,
            )?;
            let volume = geom.spacings[i].powi(3);
            let dipole = if let Some(hint) = &geom.hints[i] {
                let tensor = ellipsoid_polarisability_tensor(
                    volume, hint.depol_factors, epsilon, epsilon_m, k_peak, hint.rotation,
                );
                Dipole { position: geom.positions[i], polarisability: tensor }
            } else {
                let alpha_cm = clausius_mossotti(volume, epsilon, epsilon_m);
                let alpha = radiative_correction(alpha_cm, k_peak);
                Dipole::isotropic(geom.positions[i], alpha)
            };
            peak_dipoles.push(dipole);
        }

        match solver.solve_dipoles(&peak_dipoles, peak_wl, &peak_params) {
            Ok(response) => {
                // Estimate bounding radius of the structure for the observation plane
                let half_extent = geom.positions
                    .iter()
                    .map(|p| (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt())
                    .fold(0.0_f64, f64::max);
                let plane_half = (half_extent * 2.0).max(20.0);

                let plane = NearFieldPlane {
                    centre: [0.0, 0.0, 0.0],
                    normal: [0.0, 0.0, 1.0],
                    half_width: plane_half,
                    half_height: plane_half,
                    nx: 80,
                    ny: 80,
                };
                match solver.compute_near_field(&peak_dipoles, &response, &plane) {
                    Ok(nf) => Some(nf),
                    Err(e) => {
                        eprintln!("Warning: near-field computation failed: {}", e);
                        None
                    }
                }
            }
            Err(e) => {
                eprintln!("Warning: dipole solve for near-field failed: {}", e);
                None
            }
        }
    } else {
        None
    };

    // ── SHG sweep (optional) ──────────────────────────────────────────────────
    let shg_spectra = if let Some(nl) = &job.nonlinear {
        if nl.enable_shg {
            run_shg_sweep(
                nl,
                &wavelengths,
                &geom,
                &params,
                &solver,
                &gold,
                &silver,
                &copper,
                &tio2,
                &sio2,
            )
        } else {
            Vec::new()
        }
    } else {
        Vec::new()
    };

    // ── THG sweep (optional) ──────────────────────────────────────────────────
    let thg_spectra = if let Some(nl) = &job.nonlinear {
        if nl.enable_thg {
            run_thg_sweep(
                nl,
                &wavelengths,
                &geom,
                &params,
                &solver,
                &gold,
                &silver,
                &copper,
                &tio2,
                &sio2,
            )
        } else {
            Vec::new()
        }
    } else {
        Vec::new()
    };

    Ok(SimulationOutput { spectra: all_spectra, near_field, shg_spectra, thg_spectra })
}


/// Run the SHG sweep over all fundamental wavelengths.
///
/// For each wavelength λ:
/// 1. Build dipoles at λ (ω) and at λ/2 (2ω).
/// 2. Solve the linear CDA system at ω to get local fields.
/// 3. Drive the 2ω system with the χ^(2) source terms.
///
/// Wavelengths where 2ω material resolution fails are silently skipped.
#[allow(clippy::too_many_arguments)]
fn run_shg_sweep(
    nl: &NonlinearConfig,
    wavelengths: &[f64],
    geom: &lumina_geometry::scene::SceneDipoles,
    params: &SimulationParams,
    solver: &CdaSolver,
    gold: &lumina_materials::johnson_christy::JohnsonChristyMaterial,
    silver: &lumina_materials::johnson_christy::JohnsonChristyMaterial,
    copper: &lumina_materials::johnson_christy::JohnsonChristyMaterial,
    tio2: &lumina_materials::palik::PalikMaterial,
    sio2: &lumina_materials::palik::PalikMaterial,
) -> Vec<ShgPoint> {
    let chi2 = build_chi2_tensor(nl);
    let total_dipoles = geom.positions.len();
    let chi2_tensors = vec![chi2; total_dipoles];

    let n_wl = wavelengths.len();
    let progress = AtomicUsize::new(0);
    println!("SHG sweep: {} wavelengths...", n_wl);

    let results: Vec<Option<ShgPoint>> = wavelengths
        .par_iter()
        .map(|&wl| {
            let wl_2omega = wl / 2.0;
            let k_omega = 2.0 * std::f64::consts::PI * params.environment_n / wl;
            let k_2omega = 2.0 * std::f64::consts::PI * params.environment_n / wl_2omega;
            let epsilon_m = params.environment_n * params.environment_n;

            // Resolve substrate Fresnel factor at ω and 2ω.
            let mut params_omega = params.clone();
            let mut params_2omega = params.clone();
            if let Some(sub) = &solver.substrate {
                if let Ok(eps_sub) =
                    resolve_material_epsilon(&sub.material, wl, gold, silver, copper, tio2, sio2)
                {
                    params_omega.substrate_delta_eps =
                        Some(substrate_reflection_factor(eps_sub, epsilon_m, sub.use_retarded));
                }
                if let Ok(eps_sub) = resolve_material_epsilon(
                    &sub.material,
                    wl_2omega,
                    gold,
                    silver,
                    copper,
                    tio2,
                    sio2,
                ) {
                    params_2omega.substrate_delta_eps =
                        Some(substrate_reflection_factor(eps_sub, epsilon_m, sub.use_retarded));
                }
            }

            // Build dipoles at ω.
            let mut dipoles_omega = Vec::with_capacity(total_dipoles);
            for i in 0..total_dipoles {
                let epsilon = match resolve_material_epsilon(
                    &geom.materials[i],
                    wl,
                    gold,
                    silver,
                    copper,
                    tio2,
                    sio2,
                ) {
                    Ok(e) => e,
                    Err(_) => return None,
                };
                let volume = geom.spacings[i].powi(3);
                let dipole = if let Some(hint) = &geom.hints[i] {
                    let tensor = ellipsoid_polarisability_tensor(
                        volume,
                        hint.depol_factors,
                        epsilon,
                        epsilon_m,
                        k_omega,
                        hint.rotation,
                    );
                    Dipole { position: geom.positions[i], polarisability: tensor }
                } else {
                    let alpha = radiative_correction(
                        clausius_mossotti(volume, epsilon, epsilon_m),
                        k_omega,
                    );
                    Dipole::isotropic(geom.positions[i], alpha)
                };
                dipoles_omega.push(dipole);
            }

            // Build dipoles at 2ω (material evaluated at λ/2).
            let mut dipoles_2omega = Vec::with_capacity(total_dipoles);
            for i in 0..total_dipoles {
                let epsilon = match resolve_material_epsilon(
                    &geom.materials[i],
                    wl_2omega,
                    gold,
                    silver,
                    copper,
                    tio2,
                    sio2,
                ) {
                    Ok(e) => e,
                    Err(_) => return None, // 2ω out of material range — skip
                };
                let volume = geom.spacings[i].powi(3);
                let dipole = if let Some(hint) = &geom.hints[i] {
                    let tensor = ellipsoid_polarisability_tensor(
                        volume,
                        hint.depol_factors,
                        epsilon,
                        epsilon_m,
                        k_2omega,
                        hint.rotation,
                    );
                    Dipole { position: geom.positions[i], polarisability: tensor }
                } else {
                    let alpha = radiative_correction(
                        clausius_mossotti(volume, epsilon, epsilon_m),
                        k_2omega,
                    );
                    Dipole::isotropic(geom.positions[i], alpha)
                };
                dipoles_2omega.push(dipole);
            }

            // Linear solve at ω (need local fields for SHG sources).
            let omega_response = match solver.solve_dipoles(&dipoles_omega, wl, &params_omega) {
                Ok(r) => r,
                Err(e) => {
                    eprintln!("  SHG: linear solve failed at λ={:.1} nm: {}", wl, e);
                    return None;
                }
            };

            // SHG solve.
            let shg = match compute_shg_response(
                solver,
                &omega_response,
                &dipoles_2omega,
                &chi2_tensors,
                &params_2omega,
                nl.far_field,
            ) {
                Ok(r) => r,
                Err(e) => {
                    eprintln!("  SHG: harmonic solve failed at λ={:.1} nm: {}", wl, e);
                    return None;
                }
            };

            let done = progress.fetch_add(1, Ordering::Relaxed) + 1;
            if done.is_multiple_of(10) || done == 1 || done == n_wl {
                println!(
                    "  SHG [{}/{}] λ={:.1} nm → {:.1} nm: I_SHG={:.3e}",
                    done, n_wl, wl, shg.harmonic_nm, shg.shg_intensity
                );
            }

            Some(ShgPoint {
                fundamental_nm: shg.fundamental_nm,
                harmonic_nm: shg.harmonic_nm,
                shg_intensity: shg.shg_intensity,
            })
        })
        .collect();

    results.into_iter().flatten().collect()
}

/// Build a [`Chi2Tensor`] from the nonlinear config.
fn build_chi2_tensor(nl: &NonlinearConfig) -> Chi2Tensor {
    match nl.symmetry.as_str() {
        "isotropic_surface" => {
            let chi_zzz = nl
                .chi_zzz
                .map(|[r, i]| num_complex::Complex64::new(r, i))
                .unwrap_or(num_complex::Complex64::new(1.0, 0.0));
            let chi_zxx = nl
                .chi_zxx
                .map(|[r, i]| num_complex::Complex64::new(r, i))
                .unwrap_or(num_complex::Complex64::new(0.0, 0.0));
            Chi2Tensor::isotropic_surface(chi_zzz, chi_zxx)
        }
        _ => Chi2Tensor::zero(),
    }
}

/// Write SHG spectrum to a CSV file.
pub fn write_shg_csv(shg: &[ShgPoint], path: &Path, job: &JobConfig) -> Result<()> {
    use std::io::Write;

    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }

    let mut file = std::fs::File::create(path)?;
    writeln!(file, "# Lumina CDA Solver — SHG Spectrum")?;
    writeln!(file, "# Version: {}", env!("CARGO_PKG_VERSION"))?;
    writeln!(file, "# environment_n: {}", job.simulation.environment_n)?;
    if let Some(nl) = &job.nonlinear {
        writeln!(file, "# chi2_symmetry: {}", nl.symmetry)?;
        if let Some([r, i]) = nl.chi_zzz {
            writeln!(file, "# chi_zzz: {} + {}i nm³", r, i)?;
        }
        if let Some([r, i]) = nl.chi_zxx {
            writeln!(file, "# chi_zxx: {} + {}i nm³", r, i)?;
        }
    }
    writeln!(file, "#")?;
    writeln!(file, "fundamental_nm,harmonic_nm,shg_intensity_nm6")?;

    for pt in shg {
        writeln!(
            file,
            "{:.2},{:.2},{:.6e}",
            pt.fundamental_nm, pt.harmonic_nm, pt.shg_intensity
        )?;
    }

    println!("SHG spectrum written to: {}", path.display());
    Ok(())
}

/// Run the THG sweep over all fundamental wavelengths.
///
/// For each wavelength λ:
/// 1. Build dipoles at λ (ω) and at λ/3 (3ω).
/// 2. Solve the linear CDA system at ω to get local fields.
/// 3. Drive the 3ω system with the χ^(3) source terms.
///
/// Wavelengths where 3ω material resolution fails are silently skipped.
#[allow(clippy::too_many_arguments)]
fn run_thg_sweep(
    nl: &NonlinearConfig,
    wavelengths: &[f64],
    geom: &lumina_geometry::scene::SceneDipoles,
    params: &SimulationParams,
    solver: &CdaSolver,
    gold: &lumina_materials::johnson_christy::JohnsonChristyMaterial,
    silver: &lumina_materials::johnson_christy::JohnsonChristyMaterial,
    copper: &lumina_materials::johnson_christy::JohnsonChristyMaterial,
    tio2: &lumina_materials::palik::PalikMaterial,
    sio2: &lumina_materials::palik::PalikMaterial,
) -> Vec<ThgPoint> {
    let chi3 = build_chi3_tensor(nl);
    let total_dipoles = geom.positions.len();
    let chi3_tensors = vec![chi3; total_dipoles];

    let n_wl = wavelengths.len();
    let progress = AtomicUsize::new(0);
    println!("THG sweep: {} wavelengths...", n_wl);

    let results: Vec<Option<ThgPoint>> = wavelengths
        .par_iter()
        .map(|&wl| {
            let wl_3omega = wl / 3.0;
            let k_omega = 2.0 * std::f64::consts::PI * params.environment_n / wl;
            let k_3omega = 2.0 * std::f64::consts::PI * params.environment_n / wl_3omega;
            let epsilon_m = params.environment_n * params.environment_n;

            // Resolve substrate Fresnel factor at ω and 3ω.
            let mut params_omega = params.clone();
            let mut params_3omega = params.clone();
            if let Some(sub) = &solver.substrate {
                if let Ok(eps_sub) =
                    resolve_material_epsilon(&sub.material, wl, gold, silver, copper, tio2, sio2)
                {
                    params_omega.substrate_delta_eps =
                        Some(lumina_core::solver::substrate::substrate_reflection_factor(
                            eps_sub, epsilon_m, sub.use_retarded,
                        ));
                }
                if let Ok(eps_sub) = resolve_material_epsilon(
                    &sub.material,
                    wl_3omega,
                    gold,
                    silver,
                    copper,
                    tio2,
                    sio2,
                ) {
                    params_3omega.substrate_delta_eps =
                        Some(lumina_core::solver::substrate::substrate_reflection_factor(
                            eps_sub, epsilon_m, sub.use_retarded,
                        ));
                }
            }

            // Build dipoles at ω.
            let mut dipoles_omega = Vec::with_capacity(total_dipoles);
            for i in 0..total_dipoles {
                let epsilon = match resolve_material_epsilon(
                    &geom.materials[i],
                    wl,
                    gold,
                    silver,
                    copper,
                    tio2,
                    sio2,
                ) {
                    Ok(e) => e,
                    Err(_) => return None,
                };
                let volume = geom.spacings[i].powi(3);
                let dipole = if let Some(hint) = &geom.hints[i] {
                    let tensor = ellipsoid_polarisability_tensor(
                        volume,
                        hint.depol_factors,
                        epsilon,
                        epsilon_m,
                        k_omega,
                        hint.rotation,
                    );
                    Dipole { position: geom.positions[i], polarisability: tensor }
                } else {
                    let alpha = radiative_correction(
                        clausius_mossotti(volume, epsilon, epsilon_m),
                        k_omega,
                    );
                    Dipole::isotropic(geom.positions[i], alpha)
                };
                dipoles_omega.push(dipole);
            }

            // Build dipoles at 3ω (material evaluated at λ/3).
            let mut dipoles_3omega = Vec::with_capacity(total_dipoles);
            for i in 0..total_dipoles {
                let epsilon = match resolve_material_epsilon(
                    &geom.materials[i],
                    wl_3omega,
                    gold,
                    silver,
                    copper,
                    tio2,
                    sio2,
                ) {
                    Ok(e) => e,
                    Err(_) => return None, // 3ω out of material range — skip
                };
                let volume = geom.spacings[i].powi(3);
                let dipole = if let Some(hint) = &geom.hints[i] {
                    let tensor = ellipsoid_polarisability_tensor(
                        volume,
                        hint.depol_factors,
                        epsilon,
                        epsilon_m,
                        k_3omega,
                        hint.rotation,
                    );
                    Dipole { position: geom.positions[i], polarisability: tensor }
                } else {
                    let alpha = radiative_correction(
                        clausius_mossotti(volume, epsilon, epsilon_m),
                        k_3omega,
                    );
                    Dipole::isotropic(geom.positions[i], alpha)
                };
                dipoles_3omega.push(dipole);
            }

            // Linear solve at ω (need local fields for THG sources).
            let omega_response = match solver.solve_dipoles(&dipoles_omega, wl, &params_omega) {
                Ok(r) => r,
                Err(e) => {
                    eprintln!("  THG: linear solve failed at λ={:.1} nm: {}", wl, e);
                    return None;
                }
            };

            // THG solve.
            let thg = match compute_thg_response(
                solver,
                &omega_response,
                &dipoles_3omega,
                &chi3_tensors,
                &params_3omega,
                false,
            ) {
                Ok(r) => r,
                Err(e) => {
                    eprintln!("  THG: harmonic solve failed at λ={:.1} nm: {}", wl, e);
                    return None;
                }
            };

            let done = progress.fetch_add(1, Ordering::Relaxed) + 1;
            if done.is_multiple_of(10) || done == 1 || done == n_wl {
                println!(
                    "  THG [{}/{}] λ={:.1} nm → {:.1} nm: I_THG={:.3e}",
                    done, n_wl, wl, thg.harmonic_nm, thg.thg_intensity
                );
            }

            Some(ThgPoint {
                fundamental_nm: thg.fundamental_nm,
                harmonic_nm: thg.harmonic_nm,
                thg_intensity: thg.thg_intensity,
            })
        })
        .collect();

    results.into_iter().flatten().collect()
}

/// Build a [`Chi3Tensor`] from the nonlinear config.
fn build_chi3_tensor(nl: &NonlinearConfig) -> Chi3Tensor {
    match nl.chi3_symmetry.as_str() {
        "isotropic_bulk" => {
            let chi_xxxx = nl
                .chi3_xxxx
                .map(|[r, i]| num_complex::Complex64::new(r, i))
                .unwrap_or(num_complex::Complex64::new(1.0, 0.0));
            let chi_xxyy = nl
                .chi3_xxyy
                .map(|[r, i]| num_complex::Complex64::new(r, i))
                .unwrap_or(num_complex::Complex64::new(0.0, 0.0));
            Chi3Tensor::isotropic_bulk(chi_xxxx, chi_xxyy)
        }
        _ => Chi3Tensor::zero(),
    }
}

/// Write THG spectrum to a CSV file.
pub fn write_thg_csv(thg: &[ThgPoint], path: &Path, job: &JobConfig) -> Result<()> {
    use std::io::Write;

    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }

    let mut file = std::fs::File::create(path)?;
    writeln!(file, "# Lumina CDA Solver — THG Spectrum")?;
    writeln!(file, "# Version: {}", env!("CARGO_PKG_VERSION"))?;
    writeln!(file, "# environment_n: {}", job.simulation.environment_n)?;
    if let Some(nl) = &job.nonlinear {
        writeln!(file, "# chi3_symmetry: {}", nl.chi3_symmetry)?;
        if let Some([r, i]) = nl.chi3_xxxx {
            writeln!(file, "# chi3_xxxx: {} + {}i nm⁶", r, i)?;
        }
        if let Some([r, i]) = nl.chi3_xxyy {
            writeln!(file, "# chi3_xxyy: {} + {}i nm⁶", r, i)?;
        }
    }
    writeln!(file, "#")?;
    writeln!(file, "fundamental_nm,harmonic_nm,thg_intensity_nm9")?;

    for pt in thg {
        writeln!(
            file,
            "{:.2},{:.2},{:.6e}",
            pt.fundamental_nm, pt.harmonic_nm, pt.thg_intensity
        )?;
    }

    println!("THG spectrum written to: {}", path.display());
    Ok(())
}

/// Resolve the dielectric function for a material identifier at a given wavelength.
fn resolve_material_epsilon(
    material_id: &str,
    wavelength_nm: f64,
    gold: &JohnsonChristyMaterial,
    silver: &JohnsonChristyMaterial,
    copper: &JohnsonChristyMaterial,
    tio2: &PalikMaterial,
    sio2: &PalikMaterial,
) -> Result<Complex64> {
    let provider: &dyn MaterialProvider = match material_id {
        "Au_JC"       => gold,
        "Ag_JC"       => silver,
        "Cu_JC"       => copper,
        "TiO2_Palik"  => tio2,
        "SiO2_Palik"  => sio2,
        _ => anyhow::bail!(
            "Unknown material '{}'. Valid identifiers: Au_JC, Ag_JC, Cu_JC, TiO2_Palik, SiO2_Palik",
            material_id
        ),
    };
    provider
        .dielectric_function(wavelength_nm)
        .map_err(|e| anyhow::anyhow!("Material '{}' at {:.1} nm: {}", material_id, wavelength_nm, e))
}

/// Write cross-section spectra to a CSV file with a metadata header.
pub fn write_spectra_csv(
    spectra: &[CrossSections],
    path: &Path,
    job: &JobConfig,
) -> Result<()> {
    use std::io::Write;

    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }

    let mut file = std::fs::File::create(path)?;

    // Metadata header
    writeln!(file, "# Lumina CDA Solver — Cross-Section Spectra")?;
    writeln!(file, "# Version: {}", env!("CARGO_PKG_VERSION"))?;
    writeln!(file, "# environment_n: {}", job.simulation.environment_n)?;

    // List object materials and spacings
    for obj in &job.geometry.objects {
        writeln!(
            file,
            "# object '{}': material={}, dipole_spacing={} nm",
            obj.name, obj.material, obj.dipole_spacing
        )?;
    }
    writeln!(file, "#")?;

    // Add CD column only if present
    let has_cd = spectra.iter().any(|cs| cs.circular_dichroism.is_some());
    if has_cd {
        writeln!(file, "wavelength_nm,extinction_nm2,absorption_nm2,scattering_nm2,circular_dichroism_nm2")?;
    } else {
        writeln!(file, "wavelength_nm,extinction_nm2,absorption_nm2,scattering_nm2")?;
    }

    for cs in spectra {
        if has_cd {
            writeln!(
                file,
                "{:.2},{:.6e},{:.6e},{:.6e},{:.6e}",
                cs.wavelength_nm,
                cs.extinction,
                cs.absorption,
                cs.scattering,
                cs.circular_dichroism.unwrap_or(0.0),
            )?;
        } else {
            writeln!(
                file,
                "{:.2},{:.6e},{:.6e},{:.6e}",
                cs.wavelength_nm, cs.extinction, cs.absorption, cs.scattering
            )?;
        }
    }

    println!("Spectra written to: {}", path.display());
    Ok(())
}

/// Write cross-section spectra to a JSON file.
pub fn write_spectra_json(
    spectra: &[CrossSections],
    path: &Path,
) -> Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }

    let json = serde_json::to_string_pretty(spectra)
        .map_err(|e| anyhow::anyhow!("JSON serialisation error: {}", e))?;
    std::fs::write(path, json)?;

    println!("Spectra (JSON) written to: {}", path.display());
    Ok(())
}

/// Write a near-field map to a CSV file.
pub fn write_near_field_csv(
    near_field: &NearFieldMap,
    path: &Path,
) -> Result<()> {
    use std::io::Write;

    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }

    let mut file = std::fs::File::create(path)?;
    writeln!(file, "# Lumina CDA Solver — Near-Field Intensity Map")?;
    writeln!(file, "# Grid: {}x{}", near_field.nx, near_field.ny)?;
    writeln!(
        file,
        "# Extent: x=[{:.2}, {:.2}] y=[{:.2}, {:.2}] nm",
        near_field.extent[0], near_field.extent[1],
        near_field.extent[2], near_field.extent[3],
    )?;
    writeln!(file, "#")?;
    writeln!(file, "x_nm,y_nm,z_nm,|E|_sq")?;

    for (pos, intensity) in near_field.positions.iter().zip(near_field.field_intensity.iter()) {
        writeln!(
            file,
            "{:.4},{:.4},{:.4},{:.6e}",
            pos[0], pos[1], pos[2], intensity
        )?;
    }

    println!("Near-field map written to: {}", path.display());
    Ok(())
}

/// Create a compute backend based on the user's preference string.
///
/// - `"gpu"` — attempt GPU, fail if unavailable.
/// - `"cpu"` — always use CPU.
/// - `"auto"` (default) — try GPU, fall back to CPU.
fn create_backend(preference: &str) -> Arc<dyn ComputeBackend> {
    match preference {
        "cpu" => {
            println!("Backend: CPU");
            Arc::new(lumina_compute::CpuBackend::new())
        }
        "gpu" => {
            #[cfg(feature = "gpu")]
            {
                match lumina_compute::GpuBackend::new_blocking() {
                    Ok(gpu) => {
                        println!("Backend: {}", gpu.device_info().name);
                        Arc::new(gpu)
                    }
                    Err(e) => {
                        eprintln!("GPU requested but unavailable: {}. Aborting.", e);
                        std::process::exit(1);
                    }
                }
            }
            #[cfg(not(feature = "gpu"))]
            {
                eprintln!("GPU requested but binary was built without --features gpu. Aborting.");
                std::process::exit(1);
            }
        }
        _ => {
            // "auto" or any unrecognised value
            #[cfg(feature = "gpu")]
            {
                match lumina_compute::GpuBackend::new_blocking() {
                    Ok(gpu) => {
                        println!("Backend: {} (auto-detected)", gpu.device_info().name);
                        return Arc::new(gpu);
                    }
                    Err(e) => {
                        println!("GPU not available ({}), using CPU", e);
                    }
                }
            }
            Arc::new(lumina_compute::CpuBackend::new())
        }
    }
}
