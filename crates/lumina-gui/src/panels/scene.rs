//! Scene panel: unified multi-object assembly editor.
//!
//! Replaces the separate Geometry + Materials panels with a single panel that
//! manages a [`SceneSpec`] — a list of [`ObjectSpec`]s each with their own
//! shape, material, dipole spacing, and world-space transform.

use egui::Ui;

use lumina_geometry::scene::{ObjectSpec, SceneSpec, ShapeSpec};

/// Colors used to distinguish objects in the scatter preview.
const OBJECT_COLORS: &[egui::Color32] = &[
    egui::Color32::from_rgb(80, 140, 220),
    egui::Color32::from_rgb(220, 100, 50),
    egui::Color32::from_rgb(80, 200, 120),
    egui::Color32::from_rgb(200, 80, 200),
    egui::Color32::from_rgb(200, 200, 50),
    egui::Color32::from_rgb(80, 200, 200),
    egui::Color32::from_rgb(200, 50, 80),
    egui::Color32::from_rgb(150, 150, 150),
];

/// Projection plane for the dipole scatter preview.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProjectionPlane {
    XY,
    XZ,
    YZ,
}

/// State for the unified scene editor panel.
pub struct ScenePanel {
    /// Source-of-truth scene specification.
    pub scene: SceneSpec,
    /// Index of the currently selected object.
    pub selected: usize,
    /// Per-object cached world-space positions (None if not yet computed or error).
    cached: Vec<Option<Vec<[f64; 3]>>>,
    /// Whether any cache entry needs refreshing.
    cache_dirty: bool,
    /// Projection plane for the scatter preview.
    pub projection: ProjectionPlane,
}

impl Default for ScenePanel {
    fn default() -> Self {
        Self {
            scene: SceneSpec::default(),
            selected: 0,
            cached: vec![None],
            cache_dirty: true,
            projection: ProjectionPlane::XY,
        }
    }
}

impl ScenePanel {
    pub fn ui(&mut self, ui: &mut Ui) {
        ui.heading("Scene");
        ui.separator();

        // ── Object list row ──────────────────────────────────────────────────
        ui.horizontal(|ui| {
            ui.label("Objects:");
            for i in 0..self.scene.objects.len() {
                let color = OBJECT_COLORS[i % OBJECT_COLORS.len()];
                let label = &self.scene.objects[i].name;
                // Tiny color swatch
                let (rect, _) = ui.allocate_exact_size(
                    egui::vec2(8.0, 14.0),
                    egui::Sense::hover(),
                );
                ui.painter().rect_filled(rect, 2.0, color);
                if ui.selectable_label(i == self.selected, label.as_str()).clicked() {
                    self.selected = i;
                }
            }
        });

        // ── List management buttons ──────────────────────────────────────────
        ui.horizontal(|ui| {
            if ui.button("+ Add object").clicked() {
                let idx = self.scene.objects.len() + 1;
                let mut obj = ObjectSpec::default_sphere();
                obj.name = format!("Object {idx}");
                self.scene.objects.push(obj);
                self.selected = self.scene.objects.len() - 1;
                self.cached.push(None);
                self.cache_dirty = true;
            }

            let can_remove = self.scene.objects.len() > 1;
            ui.add_enabled_ui(can_remove, |ui| {
                if ui.button("Remove").clicked() {
                    self.scene.objects.remove(self.selected);
                    self.cached.remove(self.selected);
                    if self.selected >= self.scene.objects.len() {
                        self.selected = self.scene.objects.len() - 1;
                    }
                    self.cache_dirty = true;
                }
            });

            let can_up = self.selected > 0;
            let can_down = self.selected + 1 < self.scene.objects.len();
            ui.add_enabled_ui(can_up, |ui| {
                if ui.small_button("↑").clicked() {
                    self.scene.objects.swap(self.selected, self.selected - 1);
                    self.cached.swap(self.selected, self.selected - 1);
                    self.selected -= 1;
                }
            });
            ui.add_enabled_ui(can_down, |ui| {
                if ui.small_button("↓").clicked() {
                    self.scene.objects.swap(self.selected, self.selected + 1);
                    self.cached.swap(self.selected, self.selected + 1);
                    self.selected += 1;
                }
            });
        });

        ui.separator();

        // ── Selected object editor ───────────────────────────────────────────
        if let Some(obj) = self.scene.objects.get_mut(self.selected) {
            // Name
            ui.horizontal(|ui| {
                ui.label("Name:");
                ui.text_edit_singleline(&mut obj.name);
            });

            ui.add_space(4.0);

            // Shape type selector
            let prev_tag = shape_tag(&obj.shape);
            ui.horizontal_wrapped(|ui| {
                ui.label("Shape:");
                let mut new_tag = prev_tag;
                if ui.selectable_label(prev_tag == 0, "Sphere").clicked()    { new_tag = 0; }
                if ui.selectable_label(prev_tag == 1, "Cylinder").clicked()  { new_tag = 1; }
                if ui.selectable_label(prev_tag == 2, "Cuboid").clicked()    { new_tag = 2; }
                if ui.selectable_label(prev_tag == 3, "Ellipsoid").clicked() { new_tag = 3; }
                if ui.selectable_label(prev_tag == 4, "Helix").clicked()     { new_tag = 4; }
                if ui.selectable_label(prev_tag == 5, "File").clicked()      { new_tag = 5; }
                if new_tag != prev_tag {
                    obj.shape = default_shape(new_tag);
                    self.cache_dirty = true;
                }
            });

            ui.add_space(4.0);

            // Shape parameters
            if shape_params_ui(ui, &mut obj.shape) {
                self.cache_dirty = true;
            }

            ui.add_space(6.0);

            // Material
            ui.horizontal(|ui| {
                ui.label("Material:");
                egui::ComboBox::from_id_salt(("mat", self.selected))
                    .selected_text(&obj.material)
                    .show_ui(ui, |ui| {
                        for &mat in &["Au_JC", "Ag_JC", "Cu_JC", "TiO2_Palik", "SiO2_Palik"] {
                            ui.selectable_value(&mut obj.material, mat.to_string(), mat);
                        }
                    });
            });

            // Dipole spacing
            let prev_spacing = obj.dipole_spacing;
            ui.add(
                egui::Slider::new(&mut obj.dipole_spacing, 0.5..=10.0)
                    .text("Dipole spacing (nm)"),
            );
            if (obj.dipole_spacing - prev_spacing).abs() > 1e-9 {
                self.cache_dirty = true;
            }

            ui.add_space(6.0);
            ui.separator();

            // Transform
            ui.label("Transform:");
            let t = &mut obj.transform;
            let mut tf_dirty = false;

            ui.horizontal(|ui| {
                ui.label("Position (nm):");
                tf_dirty |= ui.add(egui::DragValue::new(&mut t.position[0]).speed(0.5).prefix("x: ")).changed();
                tf_dirty |= ui.add(egui::DragValue::new(&mut t.position[1]).speed(0.5).prefix("y: ")).changed();
                tf_dirty |= ui.add(egui::DragValue::new(&mut t.position[2]).speed(0.5).prefix("z: ")).changed();
            });
            ui.horizontal(|ui| {
                ui.label("Rotation (°):");
                tf_dirty |= ui.add(egui::DragValue::new(&mut t.rotation_deg[0]).speed(1.0).prefix("Rx: ")).changed();
                tf_dirty |= ui.add(egui::DragValue::new(&mut t.rotation_deg[1]).speed(1.0).prefix("Ry: ")).changed();
                tf_dirty |= ui.add(egui::DragValue::new(&mut t.rotation_deg[2]).speed(1.0).prefix("Rz: ")).changed();
            });
            ui.horizontal(|ui| {
                ui.label("Scale:");
                tf_dirty |= ui.add(
                    egui::DragValue::new(&mut t.scale)
                        .speed(0.01)
                        .range(0.01..=10.0_f64),
                ).changed();
            });
            if tf_dirty {
                self.cache_dirty = true;
            }

            // Per-object dipole count
            let count = self.cached.get(self.selected)
                .and_then(|c| c.as_ref())
                .map(|v| v.len())
                .unwrap_or(0);
            ui.add_space(4.0);
            ui.label(format!("Dipoles (this object): {count}"));
        }

        // ── Refresh cache ───────────────────────────────────────────────────
        if self.cache_dirty {
            self.recompute_cache();
        }

        ui.separator();

        // Total dipole count
        let total: usize = self.cached.iter()
            .filter_map(|c| c.as_ref())
            .map(|v| v.len())
            .sum();
        ui.label(format!("Total dipoles: {total}"));

        // Projection selector
        ui.horizontal(|ui| {
            ui.label("Projection:");
            ui.selectable_value(&mut self.projection, ProjectionPlane::XY, "XY");
            ui.selectable_value(&mut self.projection, ProjectionPlane::XZ, "XZ");
            ui.selectable_value(&mut self.projection, ProjectionPlane::YZ, "YZ");
        });

        // Scatter preview — one series per object, color-coded
        let (ax_a, ax_b, lbl_a, lbl_b) = match self.projection {
            ProjectionPlane::XY => (0, 1, "x (nm)", "y (nm)"),
            ProjectionPlane::XZ => (0, 2, "x (nm)", "z (nm)"),
            ProjectionPlane::YZ => (1, 2, "y (nm)", "z (nm)"),
        };

        egui_plot::Plot::new("scene_scatter")
            .height(300.0)
            .data_aspect(1.0)
            .x_axis_label(lbl_a)
            .y_axis_label(lbl_b)
            .show(ui, |plot_ui| {
                for (i, cache_entry) in self.cached.iter().enumerate() {
                    let Some(positions) = cache_entry else { continue };
                    let color = OBJECT_COLORS[i % OBJECT_COLORS.len()];
                    let name = self.scene.objects.get(i)
                        .map(|o| o.name.as_str())
                        .unwrap_or("?");
                    let pts_data: egui_plot::PlotPoints = positions
                        .iter()
                        .map(|p| [p[ax_a], p[ax_b]])
                        .collect();
                    plot_ui.points(
                        egui_plot::Points::new(pts_data)
                            .radius(2.0)
                            .color(color)
                            .name(name),
                    );
                }
            });
    }

    /// Recompute world-space position cache for each object.
    fn recompute_cache(&mut self) {
        self.cached.resize(self.scene.objects.len(), None);
        for (i, obj) in self.scene.objects.iter().enumerate() {
            self.cached[i] = obj.world_positions().ok();
        }
        self.cache_dirty = false;
    }
}

// ── Helpers ──────────────────────────────────────────────────────────────────

/// Return an integer tag identifying the shape variant (for comparisons).
fn shape_tag(shape: &ShapeSpec) -> u8 {
    match shape {
        ShapeSpec::Sphere { .. }    => 0,
        ShapeSpec::Cylinder { .. }  => 1,
        ShapeSpec::Cuboid { .. }    => 2,
        ShapeSpec::Ellipsoid { .. } => 3,
        ShapeSpec::Helix { .. }     => 4,
        ShapeSpec::File { .. }      => 5,
    }
}

/// Return the default shape for a given tag (used when switching shape type).
fn default_shape(tag: u8) -> ShapeSpec {
    match tag {
        0 => ShapeSpec::Sphere { radius: 20.0 },
        1 => ShapeSpec::Cylinder { radius: 10.0, length: 30.0, axis: [0.0, 0.0, 1.0] },
        2 => ShapeSpec::Cuboid { half_extents: [10.0, 10.0, 10.0] },
        3 => ShapeSpec::Ellipsoid { semi_axes: [15.0, 10.0, 8.0] },
        4 => ShapeSpec::Helix {
            radius: 15.0, pitch: 20.0, turns: 3.0,
            wire_radius: 4.0, axis: [0.0, 0.0, 1.0],
        },
        _ => ShapeSpec::File { path: String::new() },
    }
}

/// Render shape-specific parameter widgets; returns true if any parameter changed.
fn shape_params_ui(ui: &mut Ui, shape: &mut ShapeSpec) -> bool {
    match shape {
        ShapeSpec::Sphere { radius } => {
            slider_changed(ui, radius, 1.0..=200.0, "Radius (nm)")
        }
        ShapeSpec::Cylinder { radius, length, .. } => {
            let mut d = false;
            d |= slider_changed(ui, radius, 1.0..=200.0, "Radius (nm)");
            d |= slider_changed(ui, length, 1.0..=500.0, "Length (nm)");
            d
        }
        ShapeSpec::Cuboid { half_extents } => {
            let mut d = false;
            d |= slider_changed(ui, &mut half_extents[0], 1.0..=200.0, "Half-x (nm)");
            d |= slider_changed(ui, &mut half_extents[1], 1.0..=200.0, "Half-y (nm)");
            d |= slider_changed(ui, &mut half_extents[2], 1.0..=200.0, "Half-z (nm)");
            d
        }
        ShapeSpec::Ellipsoid { semi_axes } => {
            let mut d = false;
            d |= slider_changed(ui, &mut semi_axes[0], 1.0..=200.0, "Semi-axis a (nm)");
            d |= slider_changed(ui, &mut semi_axes[1], 1.0..=200.0, "Semi-axis b (nm)");
            d |= slider_changed(ui, &mut semi_axes[2], 1.0..=200.0, "Semi-axis c (nm)");
            d
        }
        ShapeSpec::Helix { radius, pitch, turns, wire_radius, .. } => {
            let mut d = false;
            d |= slider_changed(ui, radius, 1.0..=200.0, "Coil radius (nm)");
            d |= slider_changed(ui, wire_radius, 0.5..=50.0, "Wire radius (nm)");
            d |= slider_changed(ui, pitch, 1.0..=200.0, "Pitch (nm/turn)");
            d |= slider_changed(ui, turns, 0.5..=20.0, "Turns");
            d
        }
        ShapeSpec::File { path } => {
            let mut changed = false;
            ui.horizontal(|ui| {
                ui.label("Path:");
                if ui.text_edit_singleline(path).changed() {
                    changed = true;
                }
                if ui.button("Browse…").clicked() {
                    if let Some(picked) = rfd::FileDialog::new()
                        .set_title("Open structure file")
                        .add_filter("Structure files", &["xyz", "obj"])
                        .pick_file()
                    {
                        *path = picked.display().to_string();
                        changed = true;
                    }
                }
            });
            ui.label(
                egui::RichText::new(
                    ".xyz: atom positions in Å (converted to nm). \
                     .obj: mesh volume filled with dipoles. \
                     Set dipole spacing to the inter-atomic distance for .xyz files.",
                )
                .weak()
                .small(),
            );
            changed
        }
    }
}

/// Add a slider and return true if the value changed.
fn slider_changed(
    ui: &mut Ui,
    value: &mut f64,
    range: std::ops::RangeInclusive<f64>,
    text: &str,
) -> bool {
    let prev = *value;
    ui.add(egui::Slider::new(value, range).text(text));
    (*value - prev).abs() > 1e-9
}

