//! Geometry panel: shape builder, file import, and dipole lattice preview.

use egui::Ui;

/// State for the geometry configuration panel.
#[derive(Debug)]
pub struct GeometryPanel {
    /// Currently selected primitive type.
    pub selected_shape: ShapeType,
    /// Sphere radius (nm).
    pub sphere_radius: f64,
    /// Dipole spacing (nm).
    pub dipole_spacing: f64,
    /// Number of dipoles (computed from shape + spacing).
    pub dipole_count: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShapeType {
    Sphere,
    Cylinder,
    Cuboid,
    Ellipsoid,
    ImportFile,
}

impl Default for GeometryPanel {
    fn default() -> Self {
        Self {
            selected_shape: ShapeType::Sphere,
            sphere_radius: 20.0,
            dipole_spacing: 2.0,
            dipole_count: 0,
        }
    }
}

impl GeometryPanel {
    pub fn ui(&mut self, ui: &mut Ui) {
        ui.heading("Geometry");
        ui.separator();

        ui.horizontal(|ui| {
            ui.label("Shape:");
            ui.selectable_value(&mut self.selected_shape, ShapeType::Sphere, "Sphere");
            ui.selectable_value(&mut self.selected_shape, ShapeType::Cylinder, "Cylinder");
            ui.selectable_value(&mut self.selected_shape, ShapeType::Cuboid, "Cuboid");
            ui.selectable_value(&mut self.selected_shape, ShapeType::Ellipsoid, "Ellipsoid");
            ui.selectable_value(&mut self.selected_shape, ShapeType::ImportFile, "Import File");
        });

        ui.add_space(8.0);

        match self.selected_shape {
            ShapeType::Sphere => {
                ui.add(egui::Slider::new(&mut self.sphere_radius, 1.0..=200.0).text("Radius (nm)"));
            }
            _ => {
                ui.label("Shape parameters not yet implemented.");
            }
        }

        ui.add_space(8.0);
        ui.add(
            egui::Slider::new(&mut self.dipole_spacing, 0.5..=10.0).text("Dipole spacing (nm)"),
        );

        ui.add_space(8.0);
        ui.label(format!("Estimated dipoles: {}", self.dipole_count));

        // TODO: 3D/2D viewport for dipole lattice preview.
        ui.add_space(16.0);
        ui.group(|ui| {
            ui.set_min_size(egui::vec2(400.0, 300.0));
            ui.centered_and_justified(|ui| {
                ui.label("3D Viewport (coming soon)");
            });
        });
    }
}
