//! Main application state and egui integration.

use std::sync::mpsc;
use eframe::egui;
use num_complex::Complex64;

use lumina_core::types::{CrossSections, FarFieldMap, NearFieldMap};

use crate::panels;
use crate::panels::geometry::ShapeType;
use crate::panels::materials::MaterialChoice;

/// The main Lumina application.
pub struct LuminaApp {
    /// Which panel is currently selected in the sidebar.
    active_panel: Panel,
    /// State for the geometry panel.
    pub geometry_state: panels::geometry::GeometryPanel,
    /// State for the materials panel.
    pub materials_state: panels::materials::MaterialsPanel,
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
    },
    Error(String),
    DebugLog(String),
}

/// Sidebar navigation panels.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Panel {
    Geometry,
    Materials,
    Simulation,
    Results,
}

impl Default for LuminaApp {
    fn default() -> Self {
        Self {
            active_panel: Panel::Geometry,
            geometry_state: panels::geometry::GeometryPanel::default(),
            materials_state: panels::materials::MaterialsPanel::default(),
            simulation_state: panels::simulation::SimulationPanel::default(),
            results_state: panels::results::ResultsPanel::default(),
            result_rx: None,
        }
    }
}

impl LuminaApp {
    /// Launch a simulation on a background thread using the current panel state.
    fn launch_simulation(&mut self) {
        let (tx, rx) = mpsc::channel();
        self.result_rx = Some(rx);
        self.simulation_state.is_running = true;
        self.simulation_state.progress = 0.0;

        // --- Capture geometry parameters ---
        let shape       = self.geometry_state.selected_shape;
        let spacing     = self.geometry_state.dipole_spacing;

        // For ImportFile, capture the pre-parsed positions
        let xyz_positions: Option<Vec<[f64; 3]>> = if shape == ShapeType::ImportFile {
            self.geometry_state.cached_positions.clone()
        } else {
            None
        };

        let sphere_radius   = self.geometry_state.sphere_radius;
        let cyl_radius      = self.geometry_state.cylinder_radius;
        let cyl_length      = self.geometry_state.cylinder_length;
        let cub_hx          = self.geometry_state.cuboid_hx;
        let cub_hy          = self.geometry_state.cuboid_hy;
        let cub_hz          = self.geometry_state.cuboid_hz;
        let ell_a           = self.geometry_state.ellipsoid_a;
        let ell_b           = self.geometry_state.ellipsoid_b;
        let ell_c           = self.geometry_state.ellipsoid_c;
        let helix_r         = self.geometry_state.helix_radius;
        let helix_wr        = self.geometry_state.helix_wire_radius;
        let helix_pitch     = self.geometry_state.helix_pitch;
        let helix_turns     = self.geometry_state.helix_turns;

        // Estimate bounding half-extent for the near-field observation plane
        let shape_half_extent = match shape {
            ShapeType::Sphere    => sphere_radius,
            ShapeType::Cylinder  => cyl_radius.max(cyl_length / 2.0),
            ShapeType::Cuboid    => cub_hx.max(cub_hy).max(cub_hz),
            ShapeType::Ellipsoid => ell_a.max(ell_b).max(ell_c),
            ShapeType::Helix     => helix_r + helix_wr + helix_pitch * helix_turns / 2.0,
            ShapeType::ImportFile => 30.0,
        };

        // --- Capture material parameters ---
        let mat_choice  = self.materials_state.selected_material;
        let custom_n    = self.materials_state.custom_n;
        let custom_k    = self.materials_state.custom_k;

        // --- Capture simulation parameters ---
        let wl_start    = self.simulation_state.wavelength_start;
        let wl_end      = self.simulation_state.wavelength_end;
        let num_wl      = self.simulation_state.num_wavelengths;
        let env_n       = self.simulation_state.environment_n;
        let polarisation = self.simulation_state.polarisation;
        let compute_cd  = self.simulation_state.compute_cd;
        let use_gpu     = self.simulation_state.use_gpu;
        let debug       = self.simulation_state.debug_output;

        std::thread::spawn(move || {
            use std::sync::Arc;
            use lumina_core::fields::compute_circular_dichroism;
            use lumina_core::solver::cda::CdaSolver;
            use lumina_core::solver::{NearFieldPlane, OpticalSolver, SolverError};
            use lumina_core::types::{
                clausius_mossotti, radiative_correction, Dipole, IncidentField, SimulationParams,
            };
            use std::sync::atomic::{AtomicUsize, Ordering};
            use std::sync::Mutex;
            use rayon::prelude::*;
            use lumina_compute::ComputeBackend;
            use lumina_geometry::discretise::discretise_primitive;
            use lumina_geometry::primitives::*;
            use lumina_materials::johnson_christy::JohnsonChristyMaterial;
            use lumina_materials::palik::PalikMaterial;
            use lumina_materials::provider::MaterialProvider;
            use crate::panels::simulation::IncidentPolarisation;

            // Helper: send debug log message if enabled.
            let tx_debug = tx.clone();
            let dbg = move |msg: String| {
                if debug {
                    let _ = tx_debug.send(SimulationMessage::DebugLog(msg));
                }
            };

            // Create compute backend (GPU with CPU fallback).
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

            // Build material provider (needed by ImportFile arm before geometry is dispatched)
            let dyn_material: Option<Box<dyn MaterialProvider>> = match mat_choice {
                MaterialChoice::GoldJC    => Some(Box::new(JohnsonChristyMaterial::gold())),
                MaterialChoice::SilverJC  => Some(Box::new(JohnsonChristyMaterial::silver())),
                MaterialChoice::CopperJC  => Some(Box::new(JohnsonChristyMaterial::copper())),
                MaterialChoice::TiO2Palik => Some(Box::new(PalikMaterial::tio2())),
                MaterialChoice::SiO2Palik => Some(Box::new(PalikMaterial::sio2())),
                MaterialChoice::Custom    => None,
            };
            let custom_eps = Complex64::new(
                custom_n * custom_n - custom_k * custom_k,
                2.0 * custom_n * custom_k,
            );

            let params = SimulationParams {
                wavelength_range: [wl_start, wl_end],
                num_wavelengths: num_wl,
                environment_n: env_n,
                ..Default::default()
            };

            dbg(format!(
                "Wavelength range: {:.1}–{:.1} nm, {} points",
                wl_start, wl_end, num_wl
            ));
            dbg(format!("Environment n = {:.3}, ε_m = {:.3}", env_n, env_n * env_n));
            dbg(format!("Material: {:?}", mat_choice));
            dbg(format!("Polarisation: {:?}, CD: {}", polarisation, compute_cd));

            // Build geometry
            let primitive = match shape {
                ShapeType::Sphere => Primitive::Sphere(Sphere {
                    centre: [0.0, 0.0, 0.0],
                    radius: sphere_radius,
                }),
                ShapeType::Cylinder => Primitive::Cylinder(Cylinder {
                    base_centre: [0.0, 0.0, -cyl_length / 2.0],
                    axis: [0.0, 0.0, 1.0],
                    length: cyl_length,
                    radius: cyl_radius,
                }),
                ShapeType::Cuboid => Primitive::Cuboid(Cuboid {
                    centre: [0.0, 0.0, 0.0],
                    half_extents: [cub_hx, cub_hy, cub_hz],
                }),
                ShapeType::Ellipsoid => Primitive::Ellipsoid(Ellipsoid {
                    centre: [0.0, 0.0, 0.0],
                    semi_axes: [ell_a, ell_b, ell_c],
                }),
                ShapeType::Helix => {
                    let total_height = helix_pitch * helix_turns;
                    Primitive::Helix(Helix {
                        base_centre: [0.0, 0.0, -total_height / 2.0],
                        axis: [0.0, 0.0, 1.0],
                        radius: helix_r,
                        pitch: helix_pitch,
                        turns: helix_turns,
                        wire_radius: helix_wr,
                    })
                }
                ShapeType::ImportFile => {
                    match xyz_positions {
                        Some(ref positions_xyz) => {
                            // Skip discretisation — use raw parsed positions directly
                            let positions_copy: Vec<[f64; 3]> = positions_xyz.clone();
                            let lattice_len = positions_copy.len();

                            if lattice_len == 0 {
                                let _ = tx.send(SimulationMessage::Error(
                                    "Loaded .xyz file has no atoms.".to_string(),
                                ));
                                return;
                            }

                            // Use xyz positions directly as dipole positions
                            let half_extent = positions_copy
                                .iter()
                                .map(|p| (p[0]*p[0]+p[1]*p[1]+p[2]*p[2]).sqrt())
                                .fold(0.0_f64, f64::max);
                            let shape_half_extent_xyz = (half_extent * 2.0).max(5.0);

                            // Auto-detect nearest-neighbour distance from the positions.
                            // This is used as the FCD cell size and CM volume.
                            // Sample a few atoms to find the minimum nn distance robustly.
                            let nn_dist = {
                                let n_sample = lattice_len.min(20);
                                let mut min_d = f64::MAX;
                                for i in 0..n_sample {
                                    for j in 0..lattice_len {
                                        if i == j { continue; }
                                        let dx = positions_copy[i][0] - positions_copy[j][0];
                                        let dy = positions_copy[i][1] - positions_copy[j][1];
                                        let dz = positions_copy[i][2] - positions_copy[j][2];
                                        let d = (dx * dx + dy * dy + dz * dz).sqrt();
                                        if d > 1e-12 && d < min_d {
                                            min_d = d;
                                        }
                                    }
                                }
                                if min_d == f64::MAX { spacing } else { min_d }
                            };

                            let positions_for_thread = positions_copy;

                            // Determine epsilon_m and run simulation inline
                            let epsilon_m = env_n * env_n;
                            let mut solver_xyz = CdaSolver::with_fcd(1000, true, nn_dist);
                            solver_xyz.backend = Arc::clone(&backend);

                            let n = lattice_len;
                            let dim = 3 * n;
                            let solver_method = if n <= 1000 { "Direct (LU)" } else { "GMRES" };
                            dbg(format!("Import file: {} dipoles, matrix {}×{}", n, dim, dim));
                            dbg(format!(
                                "Auto-detected nn distance: {:.4} nm (FCD threshold: {:.4} nm)",
                                nn_dist, 2.0 * nn_dist
                            ));
                            dbg(format!("Solver: {}, FCD: on, cell_size: {:.4} nm", solver_method, nn_dist));
                            dbg(format!("Bounding radius: {:.3} nm", half_extent));

                            let wavelengths: Vec<f64> = (0..num_wl)
                                .map(|i| wl_start + (wl_end - wl_start) * i as f64 / (num_wl - 1).max(1) as f64)
                                .collect();

                            let sim_start = std::time::Instant::now();
                            let tx_par = Mutex::new(tx.clone());
                            let progress_counter = AtomicUsize::new(0);

                            dbg("Solving wavelengths in parallel (rayon)...".to_string());

                            let results: Result<Vec<(CrossSections, Vec<Dipole>)>, String> = wavelengths
                                .par_iter()
                                .map(|&wl| {
                                    let wl_start_time = std::time::Instant::now();
                                    let k = 2.0 * std::f64::consts::PI * env_n / wl;
                                    let epsilon = if let Some(ref mat) = dyn_material {
                                        mat.dielectric_function(wl).map_err(|e| format!(
                                            "Material out of range at {:.1} nm: {}", wl, e
                                        ))?
                                    } else { custom_eps };

                                    let dipoles: Vec<Dipole> = positions_for_thread.iter().map(|&pos| {
                                        let alpha_cm = clausius_mossotti(nn_dist.powi(3), epsilon, epsilon_m);
                                        let alpha = radiative_correction(alpha_cm, k);
                                        Dipole::isotropic(pos, alpha)
                                    }).collect();

                                    let cs = match solver_xyz.compute_cross_sections(&dipoles, wl, &params) {
                                        Ok(cs) => cs,
                                        Err(SolverError::ConvergenceFailure { max_iter, residual }) => {
                                            if let Ok(tx_lock) = tx_par.lock() {
                                                let _ = tx_lock.send(SimulationMessage::DebugLog(format!(
                                                    "WARNING: λ={:.1} nm did not converge after {} iterations \
                                                     (residual: {:.2e}). Consider increasing max_iterations, \
                                                     refining dipole spacing, or checking the material at this wavelength.",
                                                    wl, max_iter, residual
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
                                        Err(e) => return Err(format!("Solver error at {:.1} nm: {}", wl, e)),
                                    };

                                    let done = progress_counter.fetch_add(1, Ordering::Relaxed) + 1;
                                    if let Ok(tx_lock) = tx_par.lock() {
                                        let _ = tx_lock.send(SimulationMessage::Progress { done, total: num_wl });
                                        if debug {
                                            let elapsed = wl_start_time.elapsed();
                                            let _ = tx_lock.send(SimulationMessage::DebugLog(format!(
                                                "[{}/{}] λ={:.1} nm  ε=({:.2},{:.2})  C_ext={:.3e}  {:.1}ms",
                                                done, num_wl, wl, epsilon.re, epsilon.im,
                                                cs.extinction, elapsed.as_secs_f64() * 1000.0
                                            )));
                                        }
                                    }

                                    Ok((cs, dipoles))
                                })
                                .collect();

                            let spectra;
                            let peak_wl_xyz;
                            let peak_dipoles_xyz: Vec<Dipole>;

                            match results {
                                Ok(pairs) => {
                                    let mut peak_ext_xyz = 0.0_f64;
                                    let mut peak_idx = 0;
                                    let mut collected = Vec::with_capacity(pairs.len());
                                    for (i, (cs, _)) in pairs.iter().enumerate() {
                                        if cs.extinction > peak_ext_xyz {
                                            peak_ext_xyz = cs.extinction;
                                            peak_idx = i;
                                        }
                                        collected.push(cs.clone());
                                    }
                                    let n_failed = collected.iter()
                                        .filter(|cs| cs.extinction.is_nan()).count();
                                    spectra = collected;
                                    peak_wl_xyz = wavelengths[peak_idx];
                                    peak_dipoles_xyz = pairs.into_iter().nth(peak_idx)
                                        .map(|(_, d)| d).unwrap_or_default();
                                    if n_failed > 0 {
                                        let _ = tx.send(SimulationMessage::DebugLog(format!(
                                            "WARNING: {}/{} wavelengths did not converge. \
                                             These points appear as NaN in the spectra. \
                                             Try reducing dipole spacing or increasing solver iterations.",
                                            n_failed, num_wl
                                        )));
                                    }
                                }
                                Err(e) => {
                                    let _ = tx.send(SimulationMessage::Error(e));
                                    return;
                                }
                            }

                            dbg(format!(
                                "Spectra complete: {:.2}s total, peak C_ext at λ={:.1} nm",
                                sim_start.elapsed().as_secs_f64(),
                                peak_wl_xyz
                            ));

                            // Near/far-field at peak
                            let mut near_field = None;
                            let mut far_field = None;
                            if !peak_dipoles_xyz.is_empty() {
                                if let Ok(response) = solver_xyz.solve_dipoles(&peak_dipoles_xyz, peak_wl_xyz, &params) {
                                    let plane = NearFieldPlane {
                                        centre: [0.0, 0.0, 0.0],
                                        normal: [0.0, 0.0, 1.0],
                                        half_width: shape_half_extent_xyz,
                                        half_height: shape_half_extent_xyz,
                                        nx: 40, ny: 40,
                                    };
                                    if let Ok(nf_map) = solver_xyz.compute_near_field(&peak_dipoles_xyz, &response, &plane) {
                                        near_field = Some(nf_map);
                                    }
                                    far_field = Some(solver_xyz.compute_far_field(&peak_dipoles_xyz, &response, 36, 72));
                                }
                            }

                            let _ = tx.send(SimulationMessage::Complete { spectra, near_field, far_field });
                            return;
                        }
                        None => {
                            let _ = tx.send(SimulationMessage::Error(
                                "No .xyz file loaded. Use the Geometry panel to import one.".to_string(),
                            ));
                            return;
                        }
                    }
                }
            };

            let lattice = discretise_primitive(&primitive, spacing);
            let positions: Vec<[f64; 3]> = lattice.iter().map(|p| p.position).collect();

            let n = positions.len();
            let dim = 3 * n;
            let solver_method = if n <= 1000 { "Direct (LU)" } else { "GMRES" };
            dbg(format!("Geometry: {:?}, {} dipoles, matrix {}×{}", shape, n, dim, dim));
            dbg(format!("Solver: {}, FCD: on, spacing: {:.2} nm", solver_method, spacing));

            // Build incident fields
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

            // Select primary incident field
            let primary_incident = match polarisation {
                IncidentPolarisation::X | IncidentPolarisation::Circular => incident_x.clone(),
                IncidentPolarisation::Y => incident_y.clone(),
            };

            let mut solver = CdaSolver::with_incident(1000, true, spacing, primary_incident.clone());
            solver.backend = Arc::clone(&backend);
            let mut solver_y = CdaSolver::with_incident(1000, true, spacing, incident_y.clone());
            solver_y.backend = Arc::clone(&backend);

            let wavelengths: Vec<f64> = (0..num_wl)
                .map(|i| wl_start + (wl_end - wl_start) * i as f64 / (num_wl - 1).max(1) as f64)
                .collect();

            let epsilon_m = env_n * env_n;
            let passes = if compute_cd { 2 } else { 1 };
            let sim_start = std::time::Instant::now();
            let tx_par = Mutex::new(tx.clone());
            let progress_counter = AtomicUsize::new(0);

            dbg("Solving wavelengths in parallel (rayon)...".to_string());

            let results: Result<Vec<(CrossSections, Vec<Dipole>)>, String> = wavelengths
                .par_iter()
                .map(|&wl| {
                    let wl_start_time = std::time::Instant::now();
                    let k = 2.0 * std::f64::consts::PI * env_n / wl;

                    let epsilon = if let Some(ref mat) = dyn_material {
                        mat.dielectric_function(wl).map_err(|e| format!(
                            "Material out of range at {:.1} nm: {}", wl, e
                        ))?
                    } else {
                        custom_eps
                    };

                    let dipoles: Vec<Dipole> = positions
                        .iter()
                        .map(|&pos| {
                            let alpha_cm = clausius_mossotti(spacing.powi(3), epsilon, epsilon_m);
                            let alpha = radiative_correction(alpha_cm, k);
                            Dipole::isotropic(pos, alpha)
                        })
                        .collect();

                    let mut cs = match solver.compute_cross_sections(&dipoles, wl, &params) {
                        Ok(cs) => cs,
                        Err(SolverError::ConvergenceFailure { max_iter, residual }) => {
                            if let Ok(tx_lock) = tx_par.lock() {
                                let _ = tx_lock.send(SimulationMessage::DebugLog(format!(
                                    "WARNING: λ={:.1} nm did not converge after {} iterations \
                                     (residual: {:.2e}). Consider increasing max_iterations, \
                                     refining dipole spacing, or checking the material at this wavelength.",
                                    wl, max_iter, residual
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
                        Err(e) => return Err(format!("Solver error at {:.1} nm: {}", wl, e)),
                    };

                    // Optionally compute CD (skip if cross-sections failed to converge)
                    if compute_cd && !cs.extinction.is_nan() {
                        let resp_x = solver.solve_dipoles(&dipoles, wl, &params);
                        let resp_y = solver_y.solve_dipoles(&dipoles, wl, &params);
                        if let (Ok(rx), Ok(ry)) = (resp_x, resp_y) {
                            cs.circular_dichroism = Some(compute_circular_dichroism(&rx, &ry, k));
                        }
                    }

                    let done = progress_counter.fetch_add(1, Ordering::Relaxed) + 1;
                    if let Ok(tx_lock) = tx_par.lock() {
                        let _ = tx_lock.send(SimulationMessage::Progress {
                            done: done * passes,
                            total: num_wl * passes,
                        });
                        if debug {
                            let elapsed = wl_start_time.elapsed();
                            let _ = tx_lock.send(SimulationMessage::DebugLog(format!(
                                "[{}/{}] λ={:.1} nm  ε=({:.2},{:.2})  C_ext={:.3e}  C_abs={:.3e}  {:.1}ms",
                                done, num_wl, wl, epsilon.re, epsilon.im,
                                cs.extinction, cs.absorption,
                                elapsed.as_secs_f64() * 1000.0
                            )));
                        }
                    }

                    Ok((cs, dipoles))
                })
                .collect();

            let spectra;
            let peak_wl_idx;
            let all_dipoles_at_peak: Vec<Dipole>;

            match results {
                Ok(pairs) => {
                    let mut peak_ext = 0.0_f64;
                    let mut p_idx = 0;
                    let mut collected = Vec::with_capacity(pairs.len());
                    for (i, (cs, _)) in pairs.iter().enumerate() {
                        if cs.extinction > peak_ext {
                            peak_ext = cs.extinction;
                            p_idx = i;
                        }
                        collected.push(cs.clone());
                    }
                    let n_failed = collected.iter()
                        .filter(|cs| cs.extinction.is_nan()).count();
                    spectra = collected;
                    peak_wl_idx = p_idx;
                    all_dipoles_at_peak = pairs.into_iter().nth(p_idx)
                        .map(|(_, d)| d).unwrap_or_default();
                    if n_failed > 0 {
                        let _ = tx.send(SimulationMessage::DebugLog(format!(
                            "WARNING: {}/{} wavelengths did not converge. \
                             These points appear as NaN in the spectra. \
                             Try reducing dipole spacing or increasing solver iterations.",
                            n_failed, num_wl
                        )));
                    }
                }
                Err(e) => {
                    let _ = tx.send(SimulationMessage::Error(e));
                    return;
                }
            }

            dbg(format!(
                "Spectra complete: {:.2}s total, peak C_ext at λ={:.1} nm",
                sim_start.elapsed().as_secs_f64(),
                wavelengths[peak_wl_idx]
            ));

            // Compute near-field and far-field at the peak extinction wavelength
            let mut near_field = None;
            let mut far_field = None;

            if !all_dipoles_at_peak.is_empty() && !spectra.is_empty() {
                let peak_wl = wavelengths[peak_wl_idx];
                if let Ok(response) = solver.solve_dipoles(&all_dipoles_at_peak, peak_wl, &params) {
                    let plane = NearFieldPlane {
                        centre: [0.0, 0.0, 0.0],
                        normal: [0.0, 0.0, 1.0],
                        half_width: shape_half_extent * 2.0,
                        half_height: shape_half_extent * 2.0,
                        nx: 40,
                        ny: 40,
                    };
                    if let Ok(nf_map) = solver.compute_near_field(
                        &all_dipoles_at_peak,
                        &response,
                        &plane,
                    ) {
                        near_field = Some(nf_map);
                    }

                    let ff_map = solver.compute_far_field(
                        &all_dipoles_at_peak,
                        &response,
                        36,
                        72,
                    );
                    far_field = Some(ff_map);
                }
            }

            let _ = tx.send(SimulationMessage::Complete { spectra, near_field, far_field });
        });
    }
}

impl eframe::App for LuminaApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Poll for simulation results
        if let Some(rx) = &self.result_rx {
            while let Ok(msg) = rx.try_recv() {
                match msg {
                    SimulationMessage::Progress { done, total } => {
                        self.simulation_state.progress = done as f32 / total as f32;
                        ctx.request_repaint();
                    }
                    SimulationMessage::Complete { spectra, near_field, far_field } => {
                        self.results_state.spectra = Some(spectra);
                        self.results_state.near_field = near_field;
                        self.results_state.far_field = far_field;
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

        // Request repaint while simulation is running
        if self.simulation_state.is_running {
            ctx.request_repaint();
        }

        // Sidebar navigation
        egui::SidePanel::left("nav_panel")
            .resizable(false)
            .default_width(160.0)
            .show(ctx, |ui| {
                ui.heading("Lumina");
                ui.separator();

                ui.selectable_value(&mut self.active_panel, Panel::Geometry, "Geometry");
                ui.selectable_value(&mut self.active_panel, Panel::Materials, "Materials");
                ui.selectable_value(&mut self.active_panel, Panel::Simulation, "Simulation");
                ui.selectable_value(&mut self.active_panel, Panel::Results, "Results");
            });

        // Main content area
        let should_launch = std::cell::Cell::new(false);
        egui::CentralPanel::default().show(ctx, |ui| match self.active_panel {
            Panel::Geometry => self.geometry_state.ui(ui),
            Panel::Materials => self.materials_state.ui(ui),
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
