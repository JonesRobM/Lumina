//! Main application state and egui integration.

use std::sync::mpsc;
use eframe::egui;

use lumina_core::types::CrossSections;

use crate::panels;

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
    Complete(Vec<CrossSections>),
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
    /// Launch a simulation on a background thread.
    fn launch_simulation(&mut self) {
        let (tx, rx) = mpsc::channel();
        self.result_rx = Some(rx);
        self.simulation_state.is_running = true;
        self.simulation_state.progress = 0.0;

        let radius = self.geometry_state.sphere_radius;
        let spacing = self.geometry_state.dipole_spacing;
        let wl_start = self.simulation_state.wavelength_start;
        let wl_end = self.simulation_state.wavelength_end;
        let num_wl = self.simulation_state.num_wavelengths;
        let env_n = self.simulation_state.environment_n;

        std::thread::spawn(move || {
            use lumina_core::solver::cda::CdaSolver;
            use lumina_core::solver::OpticalSolver;
            use lumina_core::types::{
                clausius_mossotti, radiative_correction, Dipole, SimulationParams,
            };
            use lumina_geometry::discretise::discretise_primitive;
            use lumina_geometry::primitives::{Primitive, Sphere};
            use lumina_materials::johnson_christy::JohnsonChristyMaterial;
            use lumina_materials::provider::MaterialProvider;

            let sphere = Primitive::Sphere(Sphere {
                centre: [0.0, 0.0, 0.0],
                radius,
            });
            let lattice = discretise_primitive(&sphere, spacing);
            let positions: Vec<[f64; 3]> = lattice.iter().map(|p| p.position).collect();

            let gold = JohnsonChristyMaterial::gold();
            let solver = CdaSolver::default();
            let params = SimulationParams {
                wavelength_range: [wl_start, wl_end],
                num_wavelengths: num_wl,
                environment_n: env_n,
                ..Default::default()
            };

            let wavelengths: Vec<f64> = (0..num_wl)
                .map(|i| wl_start + (wl_end - wl_start) * i as f64 / (num_wl - 1).max(1) as f64)
                .collect();

            let mut spectra = Vec::new();
            let epsilon_m = env_n * env_n;

            for (wi, &wl) in wavelengths.iter().enumerate() {
                let k = 2.0 * std::f64::consts::PI * env_n / wl;
                let epsilon = match gold.dielectric_function(wl) {
                    Ok(e) => e,
                    Err(e) => {
                        let _ = tx.send(SimulationMessage::Error(format!(
                            "Material error at {:.1} nm: {}",
                            wl, e
                        )));
                        return;
                    }
                };

                let dipoles: Vec<Dipole> = positions
                    .iter()
                    .map(|&pos| {
                        let volume = spacing.powi(3);
                        let alpha_cm = clausius_mossotti(volume, epsilon, epsilon_m);
                        let alpha = radiative_correction(alpha_cm, k);
                        Dipole::isotropic(pos, alpha)
                    })
                    .collect();

                match solver.compute_cross_sections(&dipoles, wl, &params) {
                    Ok(cs) => spectra.push(cs),
                    Err(e) => {
                        let _ = tx.send(SimulationMessage::Error(format!(
                            "Solver error at {:.1} nm: {}",
                            wl, e
                        )));
                        return;
                    }
                }

                let _ = tx.send(SimulationMessage::Progress {
                    done: wi + 1,
                    total: num_wl,
                });
            }

            let _ = tx.send(SimulationMessage::Complete(spectra));
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
                    SimulationMessage::Complete(spectra) => {
                        self.results_state.spectra = Some(spectra);
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
