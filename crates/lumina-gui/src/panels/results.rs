//! Results panel: spectra display, near-field maps, and data export.

use egui::Ui;
use lumina_core::types::CrossSections;

/// State for the results display panel.
#[derive(Debug, Default)]
pub struct ResultsPanel {
    /// Whether results are available to display.
    pub has_results: bool,
    /// Which result view is active.
    pub active_view: ResultView,
    /// Computed spectra.
    pub spectra: Option<Vec<CrossSections>>,
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
            ui.label("No results yet. Configure geometry and run a simulation.");
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
                if let Some(spectra) = &self.spectra {
                    ui.label(format!("{} wavelength points computed", spectra.len()));
                    ui.add_space(4.0);

                    // Find peak extinction
                    if let Some(peak) = spectra
                        .iter()
                        .max_by(|a, b| a.extinction.partial_cmp(&b.extinction).unwrap())
                    {
                        ui.label(format!(
                            "Peak extinction: {:.2e} nm\u{00b2} at {:.1} nm",
                            peak.extinction, peak.wavelength_nm
                        ));
                    }

                    ui.add_space(8.0);

                    // Table of results
                    egui::ScrollArea::vertical()
                        .max_height(400.0)
                        .show(ui, |ui| {
                            egui::Grid::new("spectra_grid")
                                .striped(true)
                                .min_col_width(100.0)
                                .show(ui, |ui| {
                                    ui.strong("Wavelength (nm)");
                                    ui.strong("Extinction (nm\u{00b2})");
                                    ui.strong("Absorption (nm\u{00b2})");
                                    ui.strong("Scattering (nm\u{00b2})");
                                    ui.end_row();

                                    for cs in spectra {
                                        ui.label(format!("{:.1}", cs.wavelength_nm));
                                        ui.label(format!("{:.2e}", cs.extinction));
                                        ui.label(format!("{:.2e}", cs.absorption));
                                        ui.label(format!("{:.2e}", cs.scattering));
                                        ui.end_row();
                                    }
                                });
                        });

                    ui.add_space(8.0);
                    if ui.button("Export to CSV").clicked() {
                        // TODO: file dialog for save path
                        if let Err(e) = export_spectra_csv(spectra) {
                            log::error!("Failed to export: {}", e);
                        }
                    }
                }
            }
            ResultView::NearField => {
                ui.group(|ui| {
                    ui.set_min_size(egui::vec2(500.0, 350.0));
                    ui.centered_and_justified(|ui| {
                        ui.label("Near-field intensity map (coming soon)");
                    });
                });
            }
            ResultView::DipoleData => {
                ui.label("Raw dipole moments and local fields.");
                ui.label("Export functionality coming soon.");
            }
        }
    }
}

fn export_spectra_csv(spectra: &[CrossSections]) -> std::io::Result<()> {
    use std::io::Write;

    let path = "output/gui_spectra.csv";
    if let Some(parent) = std::path::Path::new(path).parent() {
        std::fs::create_dir_all(parent)?;
    }

    let mut file = std::fs::File::create(path)?;
    writeln!(file, "wavelength_nm,extinction_nm2,absorption_nm2,scattering_nm2")?;
    for cs in spectra {
        writeln!(
            file,
            "{:.2},{:.6e},{:.6e},{:.6e}",
            cs.wavelength_nm, cs.extinction, cs.absorption, cs.scattering
        )?;
    }

    Ok(())
}
