//! Results panel: spectra display, near-field maps, and data export.

use egui::Ui;

/// State for the results display panel.
#[derive(Debug, Default)]
pub struct ResultsPanel {
    /// Whether results are available to display.
    pub has_results: bool,
    /// Which result view is active.
    pub active_view: ResultView,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ResultView {
    #[default]
    Spectra,
    NearField,
    DipoleData,
}

impl ResultsPanel {
    pub fn ui(&mut self, ui: &mut Ui) {
        ui.heading("Results");
        ui.separator();

        if !self.has_results {
            ui.label("No results yet. Run a simulation first.");
            return;
        }

        ui.horizontal(|ui| {
            ui.selectable_value(&mut self.active_view, ResultView::Spectra, "Spectra");
            ui.selectable_value(&mut self.active_view, ResultView::NearField, "Near-field");
            ui.selectable_value(&mut self.active_view, ResultView::DipoleData, "Dipole data");
        });

        ui.add_space(8.0);

        match self.active_view {
            ResultView::Spectra => {
                // TODO: Plot extinction/absorption/scattering vs wavelength.
                ui.group(|ui| {
                    ui.set_min_size(egui::vec2(500.0, 350.0));
                    ui.centered_and_justified(|ui| {
                        ui.label("Cross-section spectra plot (coming soon)");
                    });
                });
            }
            ResultView::NearField => {
                // TODO: Display |E|^2 heatmap.
                ui.group(|ui| {
                    ui.set_min_size(egui::vec2(500.0, 350.0));
                    ui.centered_and_justified(|ui| {
                        ui.label("Near-field intensity map (coming soon)");
                    });
                });
            }
            ResultView::DipoleData => {
                ui.label("Raw dipole moments and local fields.");
                if ui.button("Export to CSV").clicked() {
                    // TODO: Export dipole data.
                }
            }
        }
    }
}
