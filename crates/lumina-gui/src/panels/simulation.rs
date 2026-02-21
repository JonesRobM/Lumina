//! Simulation panel: wavelength range, solver parameters, and run controls.

use egui::Ui;

/// Incident field polarisation direction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IncidentPolarisation {
    /// x-polarised (default).
    X,
    /// y-polarised.
    Y,
    /// z-propagating, circular (used for CD — internally runs x + y separately).
    Circular,
}

impl Default for IncidentPolarisation {
    fn default() -> Self {
        Self::X
    }
}

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
    /// Incident field polarisation.
    pub polarisation: IncidentPolarisation,
    /// Whether to compute circular dichroism (requires two solver runs per wavelength).
    pub compute_cd: bool,
    /// Whether a simulation is currently running.
    pub is_running: bool,
    /// Progress (0.0 to 1.0).
    pub progress: f32,
    /// Set to true by UI when the user clicks "Run".
    pub launch_requested: bool,
    /// Error message from the last run, if any.
    pub error_message: Option<String>,
}

impl Default for SimulationPanel {
    fn default() -> Self {
        Self {
            wavelength_start: 400.0,
            wavelength_end: 800.0,
            num_wavelengths: 20,
            environment_n: 1.0,
            polarisation: IncidentPolarisation::X,
            compute_cd: false,
            is_running: false,
            progress: 0.0,
            launch_requested: false,
            error_message: None,
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
        ui.add(egui::Slider::new(&mut num_wl, 5.0..=200.0).text("Wavelength points"));
        self.num_wavelengths = num_wl as usize;

        ui.add_space(8.0);
        ui.add(
            egui::Slider::new(&mut self.environment_n, 1.0..=2.5)
                .text("Medium refractive index"),
        );

        ui.add_space(12.0);
        ui.separator();

        // Polarisation selector
        ui.label("Incident polarisation:");
        ui.horizontal(|ui| {
            ui.selectable_value(&mut self.polarisation, IncidentPolarisation::X, "x-pol");
            ui.selectable_value(&mut self.polarisation, IncidentPolarisation::Y, "y-pol");
            ui.selectable_value(
                &mut self.polarisation,
                IncidentPolarisation::Circular,
                "Circular",
            );
        });

        ui.add_space(4.0);

        // CD checkbox (only meaningful with circular or combined runs)
        ui.horizontal(|ui| {
            ui.checkbox(&mut self.compute_cd, "Compute circular dichroism (ΔC_ext)");
        });
        if self.compute_cd {
            ui.label(
                egui::RichText::new("  CD requires two solver calls per wavelength.")
                    .weak()
                    .small(),
            );
        }

        ui.add_space(16.0);
        ui.separator();

        if let Some(err) = &self.error_message {
            ui.colored_label(egui::Color32::RED, format!("Error: {}", err));
            ui.add_space(4.0);
        }

        if self.is_running {
            ui.add(egui::ProgressBar::new(self.progress).text(format!(
                "Running... {:.0}%",
                self.progress * 100.0
            )));
        } else if ui.button("Run Simulation").clicked() {
            self.error_message = None;
            self.launch_requested = true;
        }
    }
}
