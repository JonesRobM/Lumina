//! Materials panel: material selection and dielectric function visualisation.

use egui::Ui;

/// State for the materials configuration panel.
#[derive(Debug)]
pub struct MaterialsPanel {
    /// Currently selected material.
    pub selected_material: MaterialChoice,
    /// Custom refractive index (real part).
    pub custom_n: f64,
    /// Custom extinction coefficient.
    pub custom_k: f64,
    /// Cached dielectric function data for plotting.
    epsilon_cache: Option<(MaterialChoice, Vec<[f64; 3]>)>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MaterialChoice {
    GoldJC,
    SilverJC,
    CopperJC,
    Custom,
}

impl Default for MaterialsPanel {
    fn default() -> Self {
        Self {
            selected_material: MaterialChoice::GoldJC,
            custom_n: 1.5,
            custom_k: 0.0,
            epsilon_cache: None,
        }
    }
}

impl MaterialsPanel {
    pub fn ui(&mut self, ui: &mut Ui) {
        ui.heading("Materials");
        ui.separator();

        let prev_material = self.selected_material;

        ui.horizontal(|ui| {
            ui.label("Material:");
            ui.selectable_value(&mut self.selected_material, MaterialChoice::GoldJC, "Au (J&C)");
            ui.selectable_value(
                &mut self.selected_material,
                MaterialChoice::SilverJC,
                "Ag (J&C)",
            );
            ui.selectable_value(
                &mut self.selected_material,
                MaterialChoice::CopperJC,
                "Cu (J&C)",
            );
            ui.selectable_value(&mut self.selected_material, MaterialChoice::Custom, "Custom");
        });

        if self.selected_material == MaterialChoice::Custom {
            ui.add_space(8.0);
            ui.add(egui::Slider::new(&mut self.custom_n, 0.01..=10.0).text("n (refractive index)"));
            ui.add(
                egui::Slider::new(&mut self.custom_k, 0.0..=20.0).text("k (extinction coeff.)"),
            );
        }

        // Invalidate cache if material changed
        if prev_material != self.selected_material {
            self.epsilon_cache = None;
        }

        // Compute and cache ε(λ) data
        if self.epsilon_cache.is_none() && self.selected_material != MaterialChoice::Custom {
            self.epsilon_cache = Some((self.selected_material, compute_epsilon_data(self.selected_material)));
        }

        // Plot ε(λ)
        ui.add_space(16.0);

        if self.selected_material == MaterialChoice::Custom {
            let eps1 = self.custom_n * self.custom_n - self.custom_k * self.custom_k;
            let eps2 = 2.0 * self.custom_n * self.custom_k;
            ui.label(format!("ε = {:.2} + {:.2}i (constant)", eps1, eps2));
        } else if let Some((_, data)) = &self.epsilon_cache {
            let eps1_points: egui_plot::PlotPoints = data
                .iter()
                .map(|d| [d[0], d[1]])
                .collect();
            let eps2_points: egui_plot::PlotPoints = data
                .iter()
                .map(|d| [d[0], d[2]])
                .collect();

            let eps1_line = egui_plot::Line::new(eps1_points)
                .name("ε₁ (real)")
                .color(egui::Color32::from_rgb(50, 120, 220))
                .width(2.0);
            let eps2_line = egui_plot::Line::new(eps2_points)
                .name("ε₂ (imag)")
                .color(egui::Color32::from_rgb(220, 100, 50))
                .width(2.0);

            egui_plot::Plot::new("epsilon_plot")
                .height(250.0)
                .x_axis_label("Wavelength (nm)")
                .y_axis_label("ε")
                .legend(egui_plot::Legend::default())
                .show(ui, |plot_ui| {
                    plot_ui.line(eps1_line);
                    plot_ui.line(eps2_line);
                });
        }
    }
}

/// Compute ε₁(λ) and ε₂(λ) data for a Johnson & Christy material.
fn compute_epsilon_data(material: MaterialChoice) -> Vec<[f64; 3]> {
    use lumina_materials::johnson_christy::JohnsonChristyMaterial;
    use lumina_materials::provider::MaterialProvider;

    let mat = match material {
        MaterialChoice::GoldJC => JohnsonChristyMaterial::gold(),
        MaterialChoice::SilverJC | MaterialChoice::CopperJC => {
            // Silver and copper use gold data for now (placeholder)
            JohnsonChristyMaterial::gold()
        }
        MaterialChoice::Custom => return Vec::new(),
    };

    let mut data = Vec::new();
    // Sample at 2nm intervals across J&C range
    let wl_start = 200.0;
    let wl_end = 890.0;
    let num_points = 346;

    for i in 0..num_points {
        let wl = wl_start + (wl_end - wl_start) * i as f64 / (num_points - 1) as f64;
        if let Ok(eps) = mat.dielectric_function(wl) {
            data.push([wl, eps.re, eps.im]);
        }
    }

    data
}
