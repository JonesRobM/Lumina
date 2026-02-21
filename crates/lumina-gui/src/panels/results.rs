//! Results panel: spectra display, near-field maps, far-field patterns, and data export.

use egui::Ui;
use lumina_core::types::{CrossSections, FarFieldMap, NearFieldMap};

/// State for the results display panel.
#[derive(Debug, Default)]
pub struct ResultsPanel {
    /// Whether results are available to display.
    pub has_results: bool,
    /// Which result view is active.
    pub active_view: ResultView,
    /// Computed spectra.
    pub spectra: Option<Vec<CrossSections>>,
    /// Near-field map (computed at peak extinction wavelength).
    pub near_field: Option<NearFieldMap>,
    /// Far-field radiation pattern (computed at peak extinction wavelength).
    pub far_field: Option<FarFieldMap>,
    /// Whether to show the data table beneath the spectra plot.
    pub show_table: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ResultView {
    #[default]
    Spectra,
    NearField,
    FarField,
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
            ui.selectable_value(&mut self.active_view, ResultView::FarField, "Far-field");
        });

        ui.add_space(8.0);

        match self.active_view {
            ResultView::Spectra => self.spectra_view(ui),
            ResultView::NearField => self.near_field_view(ui),
            ResultView::FarField => self.far_field_view(ui),
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

            let ext_points: egui_plot::PlotPoints =
                spectra.iter().map(|cs| [cs.wavelength_nm, cs.extinction]).collect();
            let abs_points: egui_plot::PlotPoints =
                spectra.iter().map(|cs| [cs.wavelength_nm, cs.absorption]).collect();
            let sca_points: egui_plot::PlotPoints =
                spectra.iter().map(|cs| [cs.wavelength_nm, cs.scattering]).collect();

            // Check if CD data is present
            let has_cd = spectra.iter().any(|cs| cs.circular_dichroism.is_some());

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

                    if has_cd {
                        let cd_points: egui_plot::PlotPoints = spectra
                            .iter()
                            .map(|cs| [cs.wavelength_nm, cs.circular_dichroism.unwrap_or(0.0)])
                            .collect();
                        let cd_line = egui_plot::Line::new(cd_points)
                            .name("\u{0394}C_ext (CD)")
                            .color(egui::Color32::from_rgb(180, 80, 200))
                            .width(2.0)
                            .style(egui_plot::LineStyle::dashed_loose());
                        plot_ui.line(cd_line);
                    }
                });

            ui.add_space(8.0);

            ui.horizontal(|ui| {
                ui.checkbox(&mut self.show_table, "Show data table");

                if ui.button("Export CSV…").clicked() {
                    export_spectra_csv_dialog(spectra);
                }

                if ui.button("Export JSON…").clicked() {
                    export_spectra_json_dialog(spectra);
                }
            });

            if self.show_table {
                ui.add_space(4.0);
                egui::ScrollArea::vertical()
                    .max_height(200.0)
                    .show(ui, |ui| {
                        egui::Grid::new("spectra_grid")
                            .striped(true)
                            .min_col_width(100.0)
                            .show(ui, |ui| {
                                ui.strong("Wavelength (nm)");
                                ui.strong("Extinction (nm\u{00b2})");
                                ui.strong("Absorption (nm\u{00b2})");
                                ui.strong("Scattering (nm\u{00b2})");
                                if has_cd {
                                    ui.strong("\u{0394}C_ext (nm\u{00b2})");
                                }
                                ui.end_row();

                                for cs in spectra {
                                    ui.label(format!("{:.1}", cs.wavelength_nm));
                                    ui.label(format!("{:.2e}", cs.extinction));
                                    ui.label(format!("{:.2e}", cs.absorption));
                                    ui.label(format!("{:.2e}", cs.scattering));
                                    if has_cd {
                                        ui.label(format!(
                                            "{:.2e}",
                                            cs.circular_dichroism.unwrap_or(0.0)
                                        ));
                                    }
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
                "Near-field map: {}×{} grid, {} points",
                nf.nx, nf.ny, nf.field_intensity.len()
            ));
            ui.add_space(4.0);

            let nx = nf.nx;
            let ny = nf.ny;

            if nx > 0 && ny > 0 && nf.field_intensity.len() == nx * ny {
                let i_min = nf.field_intensity.iter().cloned().fold(f64::INFINITY, f64::min);
                let i_max = nf.field_intensity.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                let i_range = (i_max - i_min).max(1e-30);

                let num_bins = 16usize;
                let mut bins: Vec<Vec<[f64; 2]>> = vec![Vec::new(); num_bins];

                for iy in 0..ny {
                    for ix in 0..nx {
                        let idx = iy * nx + ix;
                        let intensity = nf.field_intensity[idx];
                        let frac = ((intensity - i_min) / i_range).clamp(0.0, 1.0);
                        let bin = ((frac * (num_bins - 1) as f64).round() as usize).min(num_bins - 1);
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
                            let points = egui_plot::Points::new(plot_points).radius(4.0).color(colour);
                            plot_ui.points(points);
                        }
                    });

                ui.add_space(4.0);
                ui.label(format!("|E|\u{00b2} range: [{:.2e}, {:.2e}]", i_min, i_max));
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

    fn far_field_view(&mut self, ui: &mut Ui) {
        if let Some(ff) = &self.far_field {
            ui.label(format!(
                "Far-field pattern at {:.1} nm  ({}\u{00d7}{} = {} points)",
                ff.wavelength_nm, ff.n_theta, ff.n_phi, ff.theta.len()
            ));
            ui.add_space(4.0);

            let n_theta = ff.n_theta;
            let n_phi = ff.n_phi;

            if n_theta < 2 || n_phi < 1 {
                ui.label("Insufficient sampling for plot.");
                return;
            }

            // Normalise intensity
            let i_max = ff.intensity.iter().cloned().fold(0.0_f64, f64::max).max(1e-30);

            // Extract E-plane (φ=0, ip=0) and H-plane (φ=π/2, ip≈n_phi/4)
            let e_plane_phi_idx = 0usize;
            let h_plane_phi_idx = (n_phi / 4).min(n_phi - 1);

            let e_plane: egui_plot::PlotPoints = (0..n_theta)
                .map(|it| {
                    let idx = it * n_phi + e_plane_phi_idx;
                    let r = (ff.intensity[idx] / i_max).sqrt();
                    let theta = ff.theta[idx];
                    // Polar to Cartesian: x = r·sin(θ), y = r·cos(θ)
                    [r * theta.sin(), r * theta.cos()]
                })
                .collect();

            let h_plane: egui_plot::PlotPoints = (0..n_theta)
                .map(|it| {
                    let idx = it * n_phi + h_plane_phi_idx;
                    let r = (ff.intensity[idx] / i_max).sqrt();
                    let theta = ff.theta[idx];
                    [r * theta.sin(), r * theta.cos()]
                })
                .collect();

            // Mirror the E-plane (negative x side, φ=π, ip≈n_phi/2)
            let e_plane_back_phi_idx = (n_phi / 2).min(n_phi - 1);
            let e_plane_back: egui_plot::PlotPoints = (0..n_theta)
                .map(|it| {
                    let idx = it * n_phi + e_plane_back_phi_idx;
                    let r = (ff.intensity[idx] / i_max).sqrt();
                    let theta = ff.theta[idx];
                    [-(r * theta.sin()), r * theta.cos()]
                })
                .collect();

            let e_line = egui_plot::Line::new(e_plane)
                .name("E-plane (\u{03c6}=0\u{00b0})")
                .color(egui::Color32::from_rgb(220, 50, 50))
                .width(2.0);
            let e_line_back = egui_plot::Line::new(e_plane_back)
                .color(egui::Color32::from_rgb(220, 50, 50))
                .width(2.0);
            let h_line = egui_plot::Line::new(h_plane)
                .name("H-plane (\u{03c6}=90\u{00b0})")
                .color(egui::Color32::from_rgb(50, 120, 220))
                .width(2.0);

            egui_plot::Plot::new("farfield_plot")
                .height(350.0)
                .data_aspect(1.0)
                .x_axis_label("sin(\u{03b8})")
                .y_axis_label("cos(\u{03b8})")
                .legend(egui_plot::Legend::default())
                .show(ui, |plot_ui| {
                    plot_ui.line(e_line);
                    plot_ui.line(e_line_back);
                    plot_ui.line(h_line);
                });

            ui.add_space(4.0);
            ui.label("Polar amplitude pattern (normalised \u{221a}(I/I_max)). E-plane = xz, H-plane = yz.");
        } else {
            ui.group(|ui| {
                ui.set_min_size(egui::vec2(500.0, 350.0));
                ui.centered_and_justified(|ui| {
                    ui.label("Far-field pattern not computed. Run a simulation first.");
                });
            });
        }
    }
}

/// Open a save-file dialog and write spectra to CSV.
fn export_spectra_csv_dialog(spectra: &[CrossSections]) {
    use std::io::Write;

    let path = rfd::FileDialog::new()
        .set_title("Export spectra as CSV")
        .set_file_name("spectra.csv")
        .add_filter("CSV", &["csv"])
        .save_file();

    if let Some(path) = path {
        let result = (|| -> std::io::Result<()> {
            if let Some(parent) = path.parent() {
                std::fs::create_dir_all(parent)?;
            }
            let mut file = std::fs::File::create(&path)?;
            let has_cd = spectra.iter().any(|cs| cs.circular_dichroism.is_some());
            if has_cd {
                writeln!(file, "wavelength_nm,extinction_nm2,absorption_nm2,scattering_nm2,circular_dichroism_nm2")?;
            } else {
                writeln!(file, "wavelength_nm,extinction_nm2,absorption_nm2,scattering_nm2")?;
            }
            for cs in spectra {
                if has_cd {
                    writeln!(
                        file,
                        "{:.2},{:.6e},{:.6e},{:.6e},{:.6e}",
                        cs.wavelength_nm, cs.extinction, cs.absorption, cs.scattering,
                        cs.circular_dichroism.unwrap_or(0.0),
                    )?;
                } else {
                    writeln!(
                        file,
                        "{:.2},{:.6e},{:.6e},{:.6e}",
                        cs.wavelength_nm, cs.extinction, cs.absorption, cs.scattering,
                    )?;
                }
            }
            Ok(())
        })();
        if let Err(e) = result {
            log::error!("Failed to export CSV: {}", e);
        }
    }
}

/// Open a save-file dialog and write spectra to JSON.
fn export_spectra_json_dialog(spectra: &[CrossSections]) {
    let path = rfd::FileDialog::new()
        .set_title("Export spectra as JSON")
        .set_file_name("spectra.json")
        .add_filter("JSON", &["json"])
        .save_file();

    if let Some(path) = path {
        if let Some(parent) = path.parent() {
            let _ = std::fs::create_dir_all(parent);
        }
        match serde_json::to_string_pretty(spectra) {
            Ok(json) => {
                if let Err(e) = std::fs::write(&path, json) {
                    log::error!("Failed to write JSON: {}", e);
                }
            }
            Err(e) => {
                log::error!("Failed to serialise JSON: {}", e);
            }
        }
    }
}
