//! Results panel: spectra display, near-field maps, and data export.

use egui::Ui;
use lumina_core::types::{CrossSections, NearFieldMap};

/// State for the results display panel.
#[derive(Debug, Default)]
pub struct ResultsPanel {
    /// Whether results are available to display.
    pub has_results: bool,
    /// Which result view is active.
    pub active_view: ResultView,
    /// Computed spectra.
    pub spectra: Option<Vec<CrossSections>>,
    /// Near-field map (computed for one wavelength).
    pub near_field: Option<NearFieldMap>,
    /// Whether to show the data table beneath the spectra plot.
    pub show_table: bool,
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
            ResultView::Spectra => self.spectra_view(ui),
            ResultView::NearField => self.near_field_view(ui),
            ResultView::DipoleData => {
                ui.label("Raw dipole moments and local fields.");
                ui.label("Export functionality coming soon.");
            }
        }
    }

    fn spectra_view(&mut self, ui: &mut Ui) {
        if let Some(spectra) = &self.spectra {
            // Find peak extinction
            if let Some(peak) = spectra
                .iter()
                .max_by(|a, b| a.extinction.partial_cmp(&b.extinction).unwrap())
            {
                ui.label(format!(
                    "Peak extinction: {:.2e} nm\u{00b2} at {:.1} nm ({} points)",
                    peak.extinction, peak.wavelength_nm, spectra.len()
                ));
            }

            ui.add_space(4.0);

            // Interactive spectra plot
            let ext_points: egui_plot::PlotPoints = spectra
                .iter()
                .map(|cs| [cs.wavelength_nm, cs.extinction])
                .collect();
            let abs_points: egui_plot::PlotPoints = spectra
                .iter()
                .map(|cs| [cs.wavelength_nm, cs.absorption])
                .collect();
            let sca_points: egui_plot::PlotPoints = spectra
                .iter()
                .map(|cs| [cs.wavelength_nm, cs.scattering])
                .collect();

            let ext_line = egui_plot::Line::new(ext_points)
                .name("Extinction")
                .color(egui::Color32::from_rgb(220, 50, 50))
                .width(2.0);
            let abs_line = egui_plot::Line::new(abs_points)
                .name("Absorption")
                .color(egui::Color32::from_rgb(50, 120, 220))
                .width(2.0);
            let sca_line = egui_plot::Line::new(sca_points)
                .name("Scattering")
                .color(egui::Color32::from_rgb(50, 180, 80))
                .width(2.0);

            egui_plot::Plot::new("spectra_plot")
                .height(300.0)
                .x_axis_label("Wavelength (nm)")
                .y_axis_label("Cross-section (nm\u{00b2})")
                .legend(egui_plot::Legend::default())
                .show(ui, |plot_ui| {
                    plot_ui.line(ext_line);
                    plot_ui.line(abs_line);
                    plot_ui.line(sca_line);
                });

            ui.add_space(8.0);

            ui.horizontal(|ui| {
                ui.checkbox(&mut self.show_table, "Show data table");

                if ui.button("Export to CSV").clicked() {
                    if let Err(e) = export_spectra_csv(spectra) {
                        log::error!("Failed to export: {}", e);
                    }
                }
            });

            if self.show_table {
                ui.add_space(4.0);
                egui::ScrollArea::vertical()
                    .max_height(250.0)
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
            }
        }
    }

    fn near_field_view(&mut self, ui: &mut Ui) {
        if let Some(nf) = &self.near_field {
            ui.label(format!(
                "Near-field map: {}x{} grid, {} points",
                nf.nx, nf.ny, nf.field_intensity.len()
            ));
            ui.add_space(4.0);

            // Render as a heatmap using egui_plot
            let nx = nf.nx;
            let ny = nf.ny;

            if nx > 0 && ny > 0 && nf.field_intensity.len() == nx * ny {
                // Find intensity range for colour mapping
                let i_min = nf.field_intensity.iter().cloned().fold(f64::INFINITY, f64::min);
                let i_max = nf.field_intensity.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                let i_range = (i_max - i_min).max(1e-30);

                // Build coloured points for the heatmap
                // Bin into 16 colour levels
                let num_bins = 16usize;
                let mut bins: Vec<Vec<[f64; 2]>> = vec![Vec::new(); num_bins];

                for iy in 0..ny {
                    for ix in 0..nx {
                        let idx = iy * nx + ix;
                        let intensity = nf.field_intensity[idx];
                        let frac = ((intensity - i_min) / i_range).clamp(0.0, 1.0);
                        let bin = ((frac * (num_bins - 1) as f64).round() as usize).min(num_bins - 1);

                        // Use extent for positioning
                        let x = nf.extent[0]
                            + (nf.extent[1] - nf.extent[0]) * ix as f64 / (nx - 1).max(1) as f64;
                        let y = nf.extent[2]
                            + (nf.extent[3] - nf.extent[2]) * iy as f64 / (ny - 1).max(1) as f64;
                        bins[bin].push([x, y]);
                    }
                }

                egui_plot::Plot::new("nearfield_plot")
                    .height(350.0)
                    .data_aspect(1.0)
                    .x_axis_label("Position (nm)")
                    .y_axis_label("Position (nm)")
                    .show(ui, |plot_ui| {
                        for (bin_idx, bin_points) in bins.iter().enumerate() {
                            if bin_points.is_empty() {
                                continue;
                            }
                            let frac = bin_idx as f32 / (num_bins - 1).max(1) as f32;
                            // Hot colourmap: black → red → yellow → white
                            let r = (frac * 3.0).min(1.0);
                            let g = ((frac - 0.33) * 3.0).clamp(0.0, 1.0);
                            let b = ((frac - 0.67) * 3.0).clamp(0.0, 1.0);
                            let colour = egui::Color32::from_rgb(
                                (r * 255.0) as u8,
                                (g * 255.0) as u8,
                                (b * 255.0) as u8,
                            );

                            let plot_points: egui_plot::PlotPoints =
                                bin_points.iter().copied().collect();
                            let points = egui_plot::Points::new(plot_points)
                                .radius(4.0)
                                .color(colour);
                            plot_ui.points(points);
                        }
                    });

                ui.add_space(4.0);
                ui.label(format!(
                    "|E|² range: [{:.2e}, {:.2e}]",
                    i_min, i_max
                ));
            }
        } else {
            ui.group(|ui| {
                ui.set_min_size(egui::vec2(500.0, 350.0));
                ui.centered_and_justified(|ui| {
                    ui.label("Near-field map not computed. Run a simulation first.");
                });
            });
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
