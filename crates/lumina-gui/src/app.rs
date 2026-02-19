//! Main application state and egui integration.

use eframe::egui;

use crate::panels;

/// The main LuminaCDA application.
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
        }
    }
}

impl eframe::App for LuminaApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Sidebar navigation
        egui::SidePanel::left("nav_panel")
            .resizable(false)
            .default_width(160.0)
            .show(ctx, |ui| {
                ui.heading("LuminaCDA");
                ui.separator();

                ui.selectable_value(&mut self.active_panel, Panel::Geometry, "Geometry");
                ui.selectable_value(&mut self.active_panel, Panel::Materials, "Materials");
                ui.selectable_value(&mut self.active_panel, Panel::Simulation, "Simulation");
                ui.selectable_value(&mut self.active_panel, Panel::Results, "Results");
            });

        // Main content area
        egui::CentralPanel::default().show(ctx, |ui| match self.active_panel {
            Panel::Geometry => self.geometry_state.ui(ui),
            Panel::Materials => self.materials_state.ui(ui),
            Panel::Simulation => self.simulation_state.ui(ui),
            Panel::Results => self.results_state.ui(ui),
        });
    }
}
