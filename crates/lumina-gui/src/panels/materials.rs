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
        }
    }
}

impl MaterialsPanel {
    pub fn ui(&mut self, ui: &mut Ui) {
        ui.heading("Materials");
        ui.separator();

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

        // TODO: Plot epsilon(lambda) using egui_plot or similar.
        ui.add_space(16.0);
        ui.group(|ui| {
            ui.set_min_size(egui::vec2(400.0, 250.0));
            ui.centered_and_justified(|ui| {
                ui.label("Dielectric function plot (coming soon)");
            });
        });
    }
}
