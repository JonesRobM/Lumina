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

        std::thread::spawn(move || {
            use std::sync::Arc;
            use lumina_core::fields::compute_circular_dichroism;
            use lumina_core::solver::cda::CdaSolver;
            use lumina_core::solver::{NearFieldPlane, OpticalSolver};
            use lumina_core::types::{
                clausius_mossotti, radiative_correction, Dipole, IncidentField, SimulationParams,
            };
            use lumina_compute::ComputeBackend;
            use lumina_geometry::discretise::discretise_primitive;
            use lumina_geometry::primitives::*;
            use lumina_materials::johnson_christy::JohnsonChristyMaterial;
            use lumina_materials::palik::PalikMaterial;
            use lumina_materials::provider::MaterialProvider;
            use crate::panels::simulation::IncidentPolarisation;

            // Create compute backend (GPU with CPU fallback).
            let backend: Arc<dyn ComputeBackend> = if use_gpu {
                #[cfg(feature = "gpu")]
                {
                    match lumina_compute::GpuBackend::new_blocking() {
                        Ok(gpu) => {
                            log::info!("GPU backend: {}", gpu.device_info().name);
                            Arc::new(gpu)
                        }
                        Err(e) => {
                            log::warn!("GPU unavailable ({}), falling back to CPU", e);
                            Arc::new(lumina_compute::CpuBackend::new())
                        }
                    }
                }
                #[cfg(not(feature = "gpu"))]
                {
                    Arc::new(lumina_compute::CpuBackend::new())
                }
            } else {
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

                            // Build the solver using a synthetic lattice from the .xyz coordinates
                            // We build dipoles inline below using these positions and current material
                            let positions_for_thread = positions_copy;

                            // Determine epsilon_m and run simulation inline
                            let epsilon_m = env_n * env_n;
                            let mut solver_xyz = CdaSolver::with_fcd(1000, true, spacing);
                            solver_xyz.backend = Arc::clone(&backend);

                            let wavelengths: Vec<f64> = (0..num_wl)
                                .map(|i| wl_start + (wl_end - wl_start) * i as f64 / (num_wl - 1).max(1) as f64)
                                .collect();

                            let mut spectra = Vec::new();
                            let mut peak_ext_xyz = 0.0_f64;
                            let mut peak_wl_xyz = wavelengths[0];
                            let mut peak_dipoles_xyz: Vec<Dipole> = Vec::new();

                            for (wi, &wl) in wavelengths.iter().enumerate() {
                                let k = 2.0 * std::f64::consts::PI * env_n / wl;
                                let epsilon = if let Some(ref mat) = dyn_material {
                                    match mat.dielectric_function(wl) {
                                        Ok(e) => e,
                                        Err(e) => {
                                            let _ = tx.send(SimulationMessage::Error(format!(
                                                "Material out of range at {:.1} nm: {}", wl, e)));
                                            return;
                                        }
                                    }
                                } else { custom_eps };

                                let dipoles: Vec<Dipole> = positions_for_thread.iter().map(|&pos| {
                                    let alpha_cm = clausius_mossotti(spacing.powi(3), epsilon, epsilon_m);
                                    let alpha = radiative_correction(alpha_cm, k);
                                    Dipole::isotropic(pos, alpha)
                                }).collect();

                                match solver_xyz.compute_cross_sections(&dipoles, wl, &params) {
                                    Ok(cs) => {
                                        if cs.extinction > peak_ext_xyz {
                                            peak_ext_xyz = cs.extinction;
                                            peak_wl_xyz = wl;
                                            peak_dipoles_xyz = dipoles.clone();
                                        }
                                        spectra.push(cs);
                                    }
                                    Err(e) => {
                                        let _ = tx.send(SimulationMessage::Error(format!(
                                            "Solver error at {:.1} nm: {}", wl, e)));
                                        return;
                                    }
                                }
                                let _ = tx.send(SimulationMessage::Progress { done: wi + 1, total: num_wl });
                            }

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

            let mut spectra = Vec::new();
            let epsilon_m = env_n * env_n;

            let mut peak_wl_idx = 0;
            let mut peak_ext = 0.0_f64;
            let mut all_dipoles_at_peak = Vec::new();

            // Total steps for progress: 1 or 2 solver passes
            let passes = if compute_cd { 2 } else { 1 };

            for (wi, &wl) in wavelengths.iter().enumerate() {
                let k = 2.0 * std::f64::consts::PI * env_n / wl;

                // Get dielectric function for this wavelength
                let epsilon = if let Some(ref mat) = dyn_material {
                    match mat.dielectric_function(wl) {
                        Ok(e) => e,
                        Err(e) => {
                            let _ = tx.send(SimulationMessage::Error(format!(
                                "Material out of range at {:.1} nm: {}",
                                wl, e
                            )));
                            return;
                        }
                    }
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

                // Solve primary (x-pol or selected)
                let cs_result = solver.compute_cross_sections(&dipoles, wl, &params);

                // Optionally compute CD
                let cd_value: Option<f64> = if compute_cd {
                    let resp_x = solver.solve_dipoles(&dipoles, wl, &params);
                    let resp_y = solver_y.solve_dipoles(&dipoles, wl, &params);
                    match (resp_x, resp_y) {
                        (Ok(rx), Ok(ry)) => {
                            Some(compute_circular_dichroism(&rx, &ry, k))
                        }
                        _ => None,
                    }
                } else {
                    None
                };

                match cs_result {
                    Ok(mut cs) => {
                        cs.circular_dichroism = cd_value;
                        if cs.extinction > peak_ext {
                            peak_ext = cs.extinction;
                            peak_wl_idx = wi;
                            all_dipoles_at_peak = dipoles.clone();
                        }
                        spectra.push(cs);
                    }
                    Err(e) => {
                        let _ = tx.send(SimulationMessage::Error(format!(
                            "Solver error at {:.1} nm: {}",
                            wl, e
                        )));
                        return;
                    }
                }

                let _ = tx.send(SimulationMessage::Progress {
                    done: (wi + 1) * passes,
                    total: num_wl * passes,
                });
            }

            // Compute near-field and far-field at the peak extinction wavelength
            let mut near_field = None;
            let mut far_field = None;

            if !all_dipoles_at_peak.is_empty() && !spectra.is_empty() {
                let peak_wl = wavelengths[peak_wl_idx];
                if let Ok(response) = solver.solve_dipoles(&all_dipoles_at_peak, peak_wl, &params) {
                    // Near-field (xy plane)
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

                    // Far-field pattern (36 theta x 72 phi)
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
