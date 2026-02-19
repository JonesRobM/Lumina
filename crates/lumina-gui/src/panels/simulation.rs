//! Simulation panel: wavelength range, solver parameters, and run controls.

use egui::Ui;

/// State for the simulation configuration panel.
#[derive(Debug)]
pub struct SimulationPanel {
    /// Start wavelength (nm).
    pub wavelength_start: f64,
    /// End wavelength (nm).
    pub wavelength_end: f64,
    /// Number of wavelength points.
    pub num_wavelengths: usize,
    /// Refractive index of the surrounding medium.
    pub environment_n: f64,
    /// Whether a simulation is currently running.
    pub is_running: bool,
    /// Progress (0.0 to 1.0).
    pub progress: f32,
}

impl Default for SimulationPanel {
    fn default() -> Self {
        Self {
            wavelength_start: 400.0,
            wavelength_end: 900.0,
            num_wavelengths: 100,
            environment_n: 1.0,
            is_running: false,
            progress: 0.0,
        }
    }
}

impl SimulationPanel {
    pub fn ui(&mut self, ui: &mut Ui) {
        ui.heading("Simulation");
        ui.separator();

        ui.add(
            egui::Slider::new(&mut self.wavelength_start, 200.0..=800.0)
                .text("Start wavelength (nm)"),
        );
        ui.add(
            egui::Slider::new(&mut self.wavelength_end, 400.0..=1500.0)
                .text("End wavelength (nm)"),
        );

        let mut num_wl = self.num_wavelengths as f64;
        ui.add(egui::Slider::new(&mut num_wl, 10.0..=500.0).text("Wavelength points"));
        self.num_wavelengths = num_wl as usize;

        ui.add_space(8.0);
        ui.add(
            egui::Slider::new(&mut self.environment_n, 1.0..=2.5)
                .text("Medium refractive index"),
        );

        ui.add_space(16.0);
        ui.separator();

        if self.is_running {
            ui.add(egui::ProgressBar::new(self.progress).text("Running..."));
            if ui.button("Cancel").clicked() {
                self.is_running = false;
            }
        } else if ui.button("Run Simulation").clicked() {
            // TODO: Launch simulation on a background thread.
            self.is_running = true;
            self.progress = 0.0;
        }
    }
}
