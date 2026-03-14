//! Main application state and egui integration.

use std::sync::mpsc;
use eframe::egui;
use num_complex::Complex64;

use lumina_core::types::{CrossSections, FarFieldMap, NearFieldMap};
use lumina_geometry::scene::SceneSpec;

use crate::panels;

/// The main Lumina application.
pub struct LuminaApp {
    /// Which panel is currently selected in the sidebar.
    active_panel: Panel,
    /// Unified scene editor (replaces separate geometry + materials panels).
    pub scene_panel: panels::scene::ScenePanel,
    /// State for the simulation panel.
    pub simulation_state: panels::simulation::SimulationPanel,
    /// State for the results panel.
    pub results_state: panels::results::ResultsPanel,

    /// Channel receiver for simulation results from background thread.
    result_rx: Option<mpsc::Receiver<SimulationMessage>>,
}

/// Messages sent from the simulation background thread.
pub enum SimulationMessage {
    Progress { done: usize, total: usize },
    Complete {
        spectra: Vec<CrossSections>,
        near_field: Option<NearFieldMap>,
        far_field: Option<FarFieldMap>,
        /// SHG spectrum: `(fundamental_nm, shg_intensity)` pairs. Empty if SHG was not computed.
        shg_spectra: Vec<(f64, f64)>,
        /// THG spectrum: `(fundamental_nm, thg_intensity)` pairs. Empty if THG was not computed.
        thg_spectra: Vec<(f64, f64)>,
    },
    Error(String),
    DebugLog(String),
}

/// Sidebar navigation panels.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Panel {
    Scene,
    Simulation,
    Results,
}

impl Default for LuminaApp {
    fn default() -> Self {
        Self {
            active_panel: Panel::Scene,
            scene_panel: panels::scene::ScenePanel::default(),
            simulation_state: panels::simulation::SimulationPanel::default(),
            results_state: panels::results::ResultsPanel::default(),
            result_rx: None,
        }
    }
}

impl LuminaApp {
    /// Launch a simulation on a background thread using the current scene.
    fn launch_simulation(&mut self) {
        let (tx, rx) = mpsc::channel();
        self.result_rx = Some(rx);
        self.simulation_state.is_running = true;
        self.simulation_state.progress = 0.0;

        // Capture everything the thread needs by value.
        let scene: SceneSpec = self.scene_panel.scene.clone();

        let wl_start         = self.simulation_state.wavelength_start;
        let wl_end           = self.simulation_state.wavelength_end;
        let num_wl           = self.simulation_state.num_wavelengths;
        let env_n            = self.simulation_state.environment_n;
        let polarisation     = self.simulation_state.polarisation;
        let compute_cd       = self.simulation_state.compute_cd;
        let use_gpu          = self.simulation_state.use_gpu;
        let debug            = self.simulation_state.debug_output;
        let enable_substrate = self.simulation_state.enable_substrate;
        let substrate_z      = self.simulation_state.substrate_z_interface;
        let substrate_mat    = self.simulation_state.substrate_material.clone();
        let enable_shg       = self.simulation_state.enable_shg;
        let shg_symmetry     = self.simulation_state.shg_symmetry.clone();
        let chi_zzz_re       = self.simulation_state.chi_zzz_re;
        let chi_zzz_im       = self.simulation_state.chi_zzz_im;
        let chi_zxx_re       = self.simulation_state.chi_zxx_re;
        let chi_zxx_im       = self.simulation_state.chi_zxx_im;
        let enable_thg       = self.simulation_state.enable_thg;
        let thg_symmetry     = self.simulation_state.thg_symmetry.clone();
        let chi3_xxxx_re     = self.simulation_state.chi3_xxxx_re;
        let chi3_xxxx_im     = self.simulation_state.chi3_xxxx_im;
        let chi3_xxyy_re     = self.simulation_state.chi3_xxyy_re;
        let chi3_xxyy_im     = self.simulation_state.chi3_xxyy_im;

        std::thread::spawn(move || {
            use std::sync::Arc;
            use std::sync::atomic::{AtomicUsize, Ordering};
            use std::sync::Mutex;
            use rayon::prelude::*;

            use lumina_compute::ComputeBackend;
            use lumina_core::fields::compute_circular_dichroism;
            use lumina_core::solver::cda::CdaSolver;
            use lumina_core::solver::substrate::{SubstrateRuntime, SubstrateSpec};
            use lumina_core::solver::{NearFieldPlane, OpticalSolver, SolverError};
            use lumina_core::types::{
                clausius_mossotti, ellipsoid_polarisability_tensor, radiative_correction,
                Dipole, IncidentField, SimulationParams,
            };
            use lumina_materials::johnson_christy::JohnsonChristyMaterial;
            use lumina_materials::palik::PalikMaterial;
            use lumina_materials::provider::MaterialProvider;
            use crate::panels::simulation::IncidentPolarisation;

            // Debug logger helper.
            let tx_debug = tx.clone();
            let dbg = move |msg: String| {
                if debug {
                    let _ = tx_debug.send(SimulationMessage::DebugLog(msg));
                }
            };

            // ── Compute backend ──────────────────────────────────────────────
            let backend: Arc<dyn ComputeBackend> = if use_gpu {
                #[cfg(feature = "gpu")]
                {
                    match lumina_compute::GpuBackend::new_blocking() {
                        Ok(gpu) => {
                            dbg(format!("Backend: GPU ({})", gpu.device_info().name));
                            Arc::new(gpu)
                        }
                        Err(e) => {
                            dbg(format!("GPU unavailable ({}), falling back to CPU", e));
                            Arc::new(lumina_compute::CpuBackend::new())
                        }
                    }
                }
                #[cfg(not(feature = "gpu"))]
                {
                    dbg("Backend: CPU (GPU feature not compiled)".to_string());
                    Arc::new(lumina_compute::CpuBackend::new())
                }
            } else {
                dbg("Backend: CPU (selected)".to_string());
                Arc::new(lumina_compute::CpuBackend::new())
            };

            // ── Build geometry ───────────────────────────────────────────────
            let geom = match scene.build_geometry() {
                Ok(g) => g,
                Err(e) => {
                    let _ = tx.send(SimulationMessage::Error(format!("Geometry error: {e}")));
                    return;
                }
            };

            let n = geom.positions.len();
            if n == 0 {
                let _ = tx.send(SimulationMessage::Error(
                    "Scene produced no dipoles — check shape parameters and dipole spacing.".into(),
                ));
                return;
            }

            // Half-extent for near-field plane sizing.
            let half_extent = geom.positions
                .iter()
                .map(|p| (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt())
                .fold(0.0_f64, f64::max);

            dbg(format!(
                "Scene: {} objects, {} dipoles total, bounding radius {:.1} nm",
                scene.objects.len(), n, half_extent
            ));

            // ── Material providers ───────────────────────────────────────────
            // Created once here; borrowed by the Rayon closure below.
            let gold   = JohnsonChristyMaterial::gold();
            let silver = JohnsonChristyMaterial::silver();
            let copper = JohnsonChristyMaterial::copper();
            let tio2   = PalikMaterial::tio2();
            let sio2   = PalikMaterial::sio2();
            let al2o3  = PalikMaterial::al2o3();
            let si     = PalikMaterial::silicon();
            let gaas   = PalikMaterial::gaas();

            let resolve_eps = |mat: &str, wl: f64| -> Result<Complex64, String> {
                let provider: &dyn MaterialProvider = match mat {
                    "Au_JC"        => &gold,
                    "Ag_JC"        => &silver,
                    "Cu_JC"        => &copper,
                    "TiO2_Palik"   => &tio2,
                    "SiO2_Palik"   => &sio2,
                    "Al2O3_Palik"  => &al2o3,
                    "Si_Palik"     => &si,
                    "GaAs_Palik"   => &gaas,
                    other => return Err(format!(
                        "Unknown material '{other}'. Valid: Au_JC, Ag_JC, Cu_JC, \
                         TiO2_Palik, SiO2_Palik, Al2O3_Palik, Si_Palik, GaAs_Palik"
                    )),
                };
                provider.dielectric_function(wl).map_err(|e| format!("{e}"))
            };

            // ── Simulation parameters ────────────────────────────────────────
            let params = SimulationParams {
                wavelength_range: [wl_start, wl_end],
                num_wavelengths: num_wl,
                environment_n: env_n,
                ..Default::default()
            };

            let epsilon_m = env_n * env_n;
            let k_of = |wl: f64| 2.0 * std::f64::consts::PI * env_n / wl;

            dbg(format!(
                "Wavelengths: {wl_start:.1}–{wl_end:.1} nm, {num_wl} points, \
                 n_env = {env_n:.3}, CD = {compute_cd}"
            ));

            // ── Incident fields ──────────────────────────────────────────────
            let incident_x = IncidentField {
                direction: [0.0, 0.0, 1.0],
                polarisation: [1.0, 0.0, 0.0],
                amplitude: 1.0,
            };
            let incident_y = IncidentField {
                direction: [0.0, 0.0, 1.0],
                polarisation: [0.0, 1.0, 0.0],
                amplitude: 1.0,
            };
            let primary = match polarisation {
                IncidentPolarisation::X | IncidentPolarisation::Circular => incident_x.clone(),
                IncidentPolarisation::Y => incident_y.clone(),
            };

            // Each solver carries its own Arc<backend>.
            let substrate_spec = if enable_substrate {
                Some(SubstrateSpec { z_interface: substrate_z, material: substrate_mat.clone(), use_retarded: true, use_sommerfeld: false })
            } else {
                None
            };

            let mut solver_x = CdaSolver::with_incident(
                1000, true, geom.spacings[0], primary,
            );
            solver_x.backend = Arc::clone(&backend);
            solver_x.substrate = substrate_spec.clone();

            let mut solver_y = CdaSolver::with_incident(
                1000, true, geom.spacings[0], incident_y.clone(),
            );
            solver_y.backend = Arc::clone(&backend);
            solver_y.substrate = substrate_spec.clone();

            // ── Wavelength sweep ─────────────────────────────────────────────
            let wavelengths: Vec<f64> = (0..num_wl)
                .map(|i| wl_start + (wl_end - wl_start) * i as f64 / (num_wl - 1).max(1) as f64)
                .collect();

            let tx_par = Mutex::new(tx.clone());
            let progress_counter = AtomicUsize::new(0);
            let passes = if compute_cd { 2 } else { 1 };
            let sim_start = std::time::Instant::now();

            // Cap concurrent wavelength threads to avoid OOM.
            // On the CPU GMRES path the matrix is never assembled (on-the-fly
            // matvec), so memory per thread is O(N).  On the GPU path the
            // matrix is assembled per thread (3N×3N×16 bytes); cap threads to
            // keep that within an 8 GiB budget.
            let matrix_bytes = (3u64 * n as u64).pow(2) * 16;
            const MEMORY_BUDGET: u64 = 8 * 1024 * 1024 * 1024; // 8 GiB
            let max_par = ((MEMORY_BUDGET / matrix_bytes.max(1)) as usize)
                .max(1)
                .min(rayon::current_num_threads());
            let wl_pool = match rayon::ThreadPoolBuilder::new()
                .num_threads(max_par)
                .build()
            {
                Ok(p) => p,
                Err(e) => {
                    let _ = tx.send(SimulationMessage::Error(format!("Thread pool error: {e}")));
                    return;
                }
            };

            dbg(format!("Solving wavelengths in parallel ({max_par} threads)…"));

            let results: Result<Vec<(CrossSections, Vec<Dipole>)>, String> = wl_pool.install(|| wavelengths
                .par_iter()
                .map(|&wl| {
                    let t0 = std::time::Instant::now();
                    let k = k_of(wl);

                    // Resolve substrate runtime for this wavelength.
                    let mut params_wl = params.clone();
                    if enable_substrate {
                        match resolve_eps(&substrate_mat, wl) {
                            Ok(eps_sub) => {
                                if let Some(sub) = &substrate_spec {
                                    params_wl.substrate_runtime =
                                        Some(SubstrateRuntime::from_spec(sub, eps_sub, epsilon_m, wl));
                                }
                            }
                            Err(e) => {
                                if let Ok(lock) = tx_par.lock() {
                                    let _ = lock.send(SimulationMessage::DebugLog(format!(
                                        "Substrate material error at λ={wl:.1} nm: {e}"
                                    )));
                                }
                            }
                        }
                    }

                    // Build per-dipole polarisabilities for this wavelength.
                    let mut dipoles = Vec::with_capacity(n);
                    for i in 0..n {
                        let eps = resolve_eps(&geom.materials[i], wl).map_err(|e| {
                            format!("Material error at λ={wl:.1} nm: {e}")
                        })?;
                        let volume = geom.spacings[i].powi(3);
                        let dipole = if let Some(hint) = &geom.hints[i] {
                            // Anisotropic path: ellipsoidal depolarisation factors.
                            let tensor = ellipsoid_polarisability_tensor(
                                volume, hint.depol_factors, eps, epsilon_m, k, hint.rotation,
                            );
                            Dipole { position: geom.positions[i], polarisability: tensor }
                        } else {
                            // Isotropic path: standard spherical CM.
                            let alpha_cm = clausius_mossotti(volume, eps, epsilon_m);
                            let alpha = radiative_correction(alpha_cm, k);
                            Dipole::isotropic(geom.positions[i], alpha)
                        };
                        dipoles.push(dipole);
                    }

                    let mut cs = match solver_x.compute_cross_sections(&dipoles, wl, &params_wl) {
                        Ok(cs) => cs,
                        Err(SolverError::ConvergenceFailure { max_iter, residual }) => {
                            if let Ok(lock) = tx_par.lock() {
                                let _ = lock.send(SimulationMessage::DebugLog(format!(
                                    "WARNING: λ={wl:.1} nm did not converge after {max_iter} iterations \
                                     (residual: {residual:.2e}).",
                                )));
                            }
                            CrossSections {
                                wavelength_nm: wl,
                                extinction: f64::NAN,
                                absorption: f64::NAN,
                                scattering: f64::NAN,
                                circular_dichroism: None,
                            }
                        }
                        Err(e) => return Err(format!("Solver error at λ={wl:.1} nm: {e}")),
                    };

                    // Optional CD pass.
                    if compute_cd && !cs.extinction.is_nan() {
                        let rx = solver_x.solve_dipoles(&dipoles, wl, &params_wl);
                        let ry = solver_y.solve_dipoles(&dipoles, wl, &params_wl);
                        if let (Ok(rx), Ok(ry)) = (rx, ry) {
                            cs.circular_dichroism = Some(compute_circular_dichroism(&rx, &ry, k));
                        }
                    }

                    let done = progress_counter.fetch_add(1, Ordering::Relaxed) + 1;
                    if let Ok(lock) = tx_par.lock() {
                        let _ = lock.send(SimulationMessage::Progress {
                            done: done * passes,
                            total: num_wl * passes,
                        });
                        if debug {
                            let elapsed_ms = t0.elapsed().as_secs_f64() * 1e3;
                            let _ = lock.send(SimulationMessage::DebugLog(format!(
                                "[{done}/{num_wl}] λ={wl:.1} nm  C_ext={:.3e}  {elapsed_ms:.1} ms",
                                cs.extinction,
                            )));
                        }
                    }

                    Ok((cs, dipoles))
                })
                .collect()
            );

            // ── Collect results ──────────────────────────────────────────────
            let (spectra, peak_dipoles, peak_wl) = match results {
                Ok(pairs) => {
                    let mut peak_ext = 0.0_f64;
                    let mut peak_idx = 0;
                    let mut collected = Vec::with_capacity(pairs.len());
                    for (i, (cs, _)) in pairs.iter().enumerate() {
                        if cs.extinction > peak_ext {
                            peak_ext = cs.extinction;
                            peak_idx = i;
                        }
                        collected.push(cs.clone());
                    }
                    let n_nan = collected.iter().filter(|cs| cs.extinction.is_nan()).count();
                    if n_nan > 0 {
                        let _ = tx.send(SimulationMessage::DebugLog(format!(
                            "WARNING: {n_nan}/{num_wl} wavelengths did not converge (NaN in output).",
                        )));
                    }
                    let pw = wavelengths[peak_idx];
                    let pd = pairs.into_iter().nth(peak_idx).map(|(_, d)| d).unwrap_or_default();
                    (collected, pd, pw)
                }
                Err(e) => {
                    let _ = tx.send(SimulationMessage::Error(e));
                    return;
                }
            };

            dbg(format!(
                "Spectra complete: {:.2} s, peak C_ext at λ={peak_wl:.1} nm",
                sim_start.elapsed().as_secs_f64(),
            ));

            // ── Near-field + far-field at peak wavelength ────────────────────
            let mut near_field = None;
            let mut far_field = None;

            if !peak_dipoles.is_empty() {
                let mut peak_params = params.clone();
                if enable_substrate {
                    if let Ok(eps_sub) = resolve_eps(&substrate_mat, peak_wl) {
                        if let Some(sub) = &substrate_spec {
                            peak_params.substrate_runtime =
                                Some(SubstrateRuntime::from_spec(sub, eps_sub, epsilon_m, peak_wl));
                        }
                    }
                }
                if let Ok(response) = solver_x.solve_dipoles(&peak_dipoles, peak_wl, &peak_params) {
                    let plane_half = (half_extent * 2.0).max(20.0);
                    let plane = NearFieldPlane {
                        centre: [0.0, 0.0, 0.0],
                        normal: [0.0, 0.0, 1.0],
                        half_width: plane_half,
                        half_height: plane_half,
                        nx: 40,
                        ny: 40,
                    };
                    if let Ok(nf) = solver_x.compute_near_field(&peak_dipoles, &response, &plane) {
                        near_field = Some(nf);
                    }
                    far_field = Some(solver_x.compute_far_field(&peak_dipoles, &response, 36, 72));
                }
            }

            // ── SHG sweep (optional) ────────────────────────────────────────
            let shg_spectra: Vec<(f64, f64)> = if enable_shg {
                use lumina_core::nonlinear::compute_shg_response;
                use lumina_core::types::Chi2Tensor;

                let chi2 = match shg_symmetry.as_str() {
                    "isotropic_surface" => Chi2Tensor::isotropic_surface(
                        num_complex::Complex64::new(chi_zzz_re, chi_zzz_im),
                        num_complex::Complex64::new(chi_zxx_re, chi_zxx_im),
                    ),
                    _ => Chi2Tensor::zero(),
                };
                let chi2_tensors = vec![chi2; n];

                dbg("SHG sweep…".to_string());

                let shg_results: Vec<Option<(f64, f64)>> = wl_pool.install(|| wavelengths
                    .par_iter()
                    .map(|&wl| {
                        let wl_2omega = wl / 2.0;
                        let k_omega = k_of(wl);
                        let k_2omega = k_of(wl_2omega);

                        // Build dipoles at ω.
                        let mut dipoles_omega = Vec::with_capacity(n);
                        for i in 0..n {
                            let eps = resolve_eps(&geom.materials[i], wl).ok()?;
                            let volume = geom.spacings[i].powi(3);
                            let dipole = if let Some(hint) = &geom.hints[i] {
                                let tensor = ellipsoid_polarisability_tensor(
                                    volume, hint.depol_factors, eps, epsilon_m, k_omega,
                                    hint.rotation,
                                );
                                Dipole { position: geom.positions[i], polarisability: tensor }
                            } else {
                                let alpha = radiative_correction(
                                    clausius_mossotti(volume, eps, epsilon_m),
                                    k_omega,
                                );
                                Dipole::isotropic(geom.positions[i], alpha)
                            };
                            dipoles_omega.push(dipole);
                        }

                        // Build dipoles at 2ω (skip wavelength if material out of range).
                        let mut dipoles_2omega = Vec::with_capacity(n);
                        for i in 0..n {
                            let eps = resolve_eps(&geom.materials[i], wl_2omega).ok()?;
                            let volume = geom.spacings[i].powi(3);
                            let dipole = if let Some(hint) = &geom.hints[i] {
                                let tensor = ellipsoid_polarisability_tensor(
                                    volume, hint.depol_factors, eps, epsilon_m, k_2omega,
                                    hint.rotation,
                                );
                                Dipole { position: geom.positions[i], polarisability: tensor }
                            } else {
                                let alpha = radiative_correction(
                                    clausius_mossotti(volume, eps, epsilon_m),
                                    k_2omega,
                                );
                                Dipole::isotropic(geom.positions[i], alpha)
                            };
                            dipoles_2omega.push(dipole);
                        }

                        // Linear solve at ω.
                        let omega_response =
                            solver_x.solve_dipoles(&dipoles_omega, wl, &params).ok()?;

                        // SHG solve at 2ω.
                        let shg = compute_shg_response(
                            &solver_x,
                            &omega_response,
                            &dipoles_2omega,
                            &chi2_tensors,
                            &params,
                            false,
                        )
                        .ok()?;

                        Some((wl, shg.shg_intensity))
                    })
                    .collect()
                );

                shg_results.into_iter().flatten().collect()
            } else {
                Vec::new()
            };

            // ── THG sweep (optional) ────────────────────────────────────────
            let thg_spectra: Vec<(f64, f64)> = if enable_thg {
                use lumina_core::nonlinear::compute_thg_response;
                use lumina_core::types::Chi3Tensor;

                let chi3 = match thg_symmetry.as_str() {
                    "isotropic_bulk" => Chi3Tensor::isotropic_bulk(
                        num_complex::Complex64::new(chi3_xxxx_re, chi3_xxxx_im),
                        num_complex::Complex64::new(chi3_xxyy_re, chi3_xxyy_im),
                    ),
                    _ => Chi3Tensor::zero(),
                };
                let chi3_tensors = vec![chi3; n];

                dbg("THG sweep…".to_string());

                let thg_results: Vec<Option<(f64, f64)>> = wl_pool.install(|| wavelengths
                    .par_iter()
                    .map(|&wl| {
                        let wl_3omega = wl / 3.0;
                        let k_omega = k_of(wl);
                        let k_3omega = k_of(wl_3omega);

                        // Build dipoles at ω.
                        let mut dipoles_omega = Vec::with_capacity(n);
                        for i in 0..n {
                            let eps = resolve_eps(&geom.materials[i], wl).ok()?;
                            let volume = geom.spacings[i].powi(3);
                            let dipole = if let Some(hint) = &geom.hints[i] {
                                let tensor = ellipsoid_polarisability_tensor(
                                    volume, hint.depol_factors, eps, epsilon_m, k_omega,
                                    hint.rotation,
                                );
                                Dipole { position: geom.positions[i], polarisability: tensor }
                            } else {
                                let alpha = radiative_correction(
                                    clausius_mossotti(volume, eps, epsilon_m),
                                    k_omega,
                                );
                                Dipole::isotropic(geom.positions[i], alpha)
                            };
                            dipoles_omega.push(dipole);
                        }

                        // Build dipoles at 3ω (skip wavelength if material out of range).
                        let mut dipoles_3omega = Vec::with_capacity(n);
                        for i in 0..n {
                            let eps = resolve_eps(&geom.materials[i], wl_3omega).ok()?;
                            let volume = geom.spacings[i].powi(3);
                            let dipole = if let Some(hint) = &geom.hints[i] {
                                let tensor = ellipsoid_polarisability_tensor(
                                    volume, hint.depol_factors, eps, epsilon_m, k_3omega,
                                    hint.rotation,
                                );
                                Dipole { position: geom.positions[i], polarisability: tensor }
                            } else {
                                let alpha = radiative_correction(
                                    clausius_mossotti(volume, eps, epsilon_m),
                                    k_3omega,
                                );
                                Dipole::isotropic(geom.positions[i], alpha)
                            };
                            dipoles_3omega.push(dipole);
                        }

                        // Linear solve at ω.
                        let omega_response =
                            solver_x.solve_dipoles(&dipoles_omega, wl, &params).ok()?;

                        // THG solve at 3ω.
                        let thg = compute_thg_response(
                            &solver_x,
                            &omega_response,
                            &dipoles_3omega,
                            &chi3_tensors,
                            &params,
                            false,
                        )
                        .ok()?;

                        Some((wl, thg.thg_intensity))
                    })
                    .collect()
                );

                thg_results.into_iter().flatten().collect()
            } else {
                Vec::new()
            };

            let _ = tx.send(SimulationMessage::Complete { spectra, near_field, far_field, shg_spectra, thg_spectra });
        });
    }
}

impl eframe::App for LuminaApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Poll for simulation results.
        if let Some(rx) = &self.result_rx {
            while let Ok(msg) = rx.try_recv() {
                match msg {
                    SimulationMessage::Progress { done, total } => {
                        self.simulation_state.progress = done as f32 / total as f32;
                        ctx.request_repaint();
                    }
                    SimulationMessage::Complete { spectra, near_field, far_field, shg_spectra, thg_spectra } => {
                        self.results_state.spectra = Some(spectra);
                        self.results_state.near_field = near_field;
                        self.results_state.far_field = far_field;
                        self.results_state.shg_spectra =
                            if shg_spectra.is_empty() { None } else { Some(shg_spectra) };
                        self.results_state.thg_spectra =
                            if thg_spectra.is_empty() { None } else { Some(thg_spectra) };
                        self.results_state.has_results = true;
                        self.simulation_state.is_running = false;
                        self.simulation_state.progress = 1.0;
                        self.active_panel = Panel::Results;
                    }
                    SimulationMessage::Error(e) => {
                        self.simulation_state.is_running = false;
                        self.simulation_state.error_message = Some(e);
                    }
                    SimulationMessage::DebugLog(msg) => {
                        self.simulation_state.debug_log.push(msg);
                    }
                }
            }
        }

        // Request repaint while simulation is running.
        if self.simulation_state.is_running {
            ctx.request_repaint();
        }

        // Sidebar navigation.
        egui::SidePanel::left("nav_panel")
            .resizable(false)
            .default_width(140.0)
            .show(ctx, |ui| {
                ui.heading("Lumina");
                ui.separator();
                ui.selectable_value(&mut self.active_panel, Panel::Scene, "Scene");
                ui.selectable_value(&mut self.active_panel, Panel::Simulation, "Simulation");
                ui.selectable_value(&mut self.active_panel, Panel::Results, "Results");
            });

        // Main content area.
        let should_launch = std::cell::Cell::new(false);
        egui::CentralPanel::default().show(ctx, |ui| match self.active_panel {
            Panel::Scene => self.scene_panel.ui(ui),
            Panel::Simulation => {
                self.simulation_state.ui(ui);
                if self.simulation_state.launch_requested {
                    self.simulation_state.launch_requested = false;
                    should_launch.set(true);
                }
            }
            Panel::Results => self.results_state.ui(ui),
        });

        if should_launch.get() {
            self.launch_simulation();
        }
    }
}
