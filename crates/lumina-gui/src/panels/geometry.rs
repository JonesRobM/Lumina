//! Geometry panel: shape builder, file import, and dipole lattice preview.

use egui::Ui;
use rfd;

/// An entry in the directory file listing.
#[derive(Debug, Clone)]
pub struct ImportFileEntry {
    /// Full path to the file.
    pub path: std::path::PathBuf,
    /// Filename without the directory prefix.
    pub filename: String,
    /// File extension ("xyz" or "obj").
    pub extension: String,
}

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

    // Cylinder parameters
    /// Cylinder radius (nm).
    pub cylinder_radius: f64,
    /// Cylinder length (nm).
    pub cylinder_length: f64,

    // Cuboid parameters
    /// Cuboid half-extent in x (nm).
    pub cuboid_hx: f64,
    /// Cuboid half-extent in y (nm).
    pub cuboid_hy: f64,
    /// Cuboid half-extent in z (nm).
    pub cuboid_hz: f64,

    // Ellipsoid parameters
    /// Ellipsoid semi-axis a (nm).
    pub ellipsoid_a: f64,
    /// Ellipsoid semi-axis b (nm).
    pub ellipsoid_b: f64,
    /// Ellipsoid semi-axis c (nm).
    pub ellipsoid_c: f64,

    // Helix parameters
    /// Helix coil radius (nm).
    pub helix_radius: f64,
    /// Wire cross-section radius (nm).
    pub helix_wire_radius: f64,
    /// Rise per full turn (nm).
    pub helix_pitch: f64,
    /// Number of turns.
    pub helix_turns: f64,

    /// Which projection plane to show in the 2D scatter preview.
    pub projection_plane: ProjectionPlane,

    /// Cached dipole positions for the scatter preview.
    pub cached_positions: Option<Vec<[f64; 3]>>,
    /// Whether the cache needs refreshing.
    cache_dirty: bool,

    /// Path to the loaded .xyz file (for display).
    pub xyz_path_display: Option<String>,
    /// Error message from last .xyz parse attempt.
    pub xyz_error: Option<String>,

    /// Loaded OBJ mesh (stored so it can be re-discretised when spacing changes).
    pub obj_mesh: Option<lumina_geometry::parsers::obj::ObjMesh>,
    /// Path to the loaded .obj file (for display).
    pub obj_path_display: Option<String>,
    /// Error message from last .obj parse attempt.
    pub obj_error: Option<String>,

    /// Directory path for browsing structure files.
    pub import_directory: Option<std::path::PathBuf>,
    /// Display string for the selected directory.
    pub import_dir_display: Option<String>,
    /// Cached list of .xyz and .obj files in the import directory.
    pub available_files: Vec<ImportFileEntry>,
    /// Currently selected filename from the dropdown.
    pub selected_import_file: Option<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShapeType {
    Sphere,
    Cylinder,
    Cuboid,
    Ellipsoid,
    Helix,
    ImportFile,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProjectionPlane {
    XY,
    XZ,
    YZ,
}

impl Default for GeometryPanel {
    fn default() -> Self {
        Self {
            selected_shape: ShapeType::Sphere,
            sphere_radius: 20.0,
            dipole_spacing: 2.0,
            dipole_count: 0,
            cylinder_radius: 10.0,
            cylinder_length: 30.0,
            cuboid_hx: 10.0,
            cuboid_hy: 10.0,
            cuboid_hz: 10.0,
            ellipsoid_a: 15.0,
            ellipsoid_b: 10.0,
            ellipsoid_c: 8.0,
            helix_radius: 15.0,
            helix_wire_radius: 4.0,
            helix_pitch: 20.0,
            helix_turns: 3.0,
            projection_plane: ProjectionPlane::XY,
            cached_positions: None,
            cache_dirty: true,
            xyz_path_display: None,
            xyz_error: None,
            obj_mesh: None,
            obj_path_display: None,
            obj_error: None,
            import_directory: None,
            import_dir_display: None,
            available_files: Vec::new(),
            selected_import_file: None,
        }
    }
}

impl GeometryPanel {
    pub fn ui(&mut self, ui: &mut Ui) {
        ui.heading("Geometry");
        ui.separator();

        let prev_shape = self.selected_shape;

        ui.horizontal(|ui| {
            ui.label("Shape:");
            ui.selectable_value(&mut self.selected_shape, ShapeType::Sphere, "Sphere");
            ui.selectable_value(&mut self.selected_shape, ShapeType::Cylinder, "Cylinder");
            ui.selectable_value(&mut self.selected_shape, ShapeType::Cuboid, "Cuboid");
            ui.selectable_value(&mut self.selected_shape, ShapeType::Ellipsoid, "Ellipsoid");
            ui.selectable_value(&mut self.selected_shape, ShapeType::Helix, "Helix");
            ui.selectable_value(&mut self.selected_shape, ShapeType::ImportFile, "Import File");
        });

        if prev_shape != self.selected_shape {
            self.cache_dirty = true;
        }

        ui.add_space(8.0);

        match self.selected_shape {
            ShapeType::Sphere => {
                if slider_changed(ui, &mut self.sphere_radius, 1.0..=200.0, "Radius (nm)") {
                    self.cache_dirty = true;
                }
            }
            ShapeType::Cylinder => {
                if slider_changed(ui, &mut self.cylinder_radius, 1.0..=100.0, "Radius (nm)") {
                    self.cache_dirty = true;
                }
                if slider_changed(ui, &mut self.cylinder_length, 1.0..=200.0, "Length (nm)") {
                    self.cache_dirty = true;
                }
            }
            ShapeType::Cuboid => {
                if slider_changed(ui, &mut self.cuboid_hx, 1.0..=100.0, "Half-extent X (nm)") {
                    self.cache_dirty = true;
                }
                if slider_changed(ui, &mut self.cuboid_hy, 1.0..=100.0, "Half-extent Y (nm)") {
                    self.cache_dirty = true;
                }
                if slider_changed(ui, &mut self.cuboid_hz, 1.0..=100.0, "Half-extent Z (nm)") {
                    self.cache_dirty = true;
                }
            }
            ShapeType::Ellipsoid => {
                if slider_changed(ui, &mut self.ellipsoid_a, 1.0..=100.0, "Semi-axis a (nm)") {
                    self.cache_dirty = true;
                }
                if slider_changed(ui, &mut self.ellipsoid_b, 1.0..=100.0, "Semi-axis b (nm)") {
                    self.cache_dirty = true;
                }
                if slider_changed(ui, &mut self.ellipsoid_c, 1.0..=100.0, "Semi-axis c (nm)") {
                    self.cache_dirty = true;
                }
            }
            ShapeType::Helix => {
                if slider_changed(ui, &mut self.helix_radius, 1.0..=100.0, "Coil radius (nm)") {
                    self.cache_dirty = true;
                }
                if slider_changed(ui, &mut self.helix_wire_radius, 0.5..=20.0, "Wire radius (nm)") {
                    self.cache_dirty = true;
                }
                if slider_changed(ui, &mut self.helix_pitch, 1.0..=100.0, "Pitch (nm/turn)") {
                    self.cache_dirty = true;
                }
                if slider_changed(ui, &mut self.helix_turns, 0.5..=10.0, "Turns") {
                    self.cache_dirty = true;
                }
            }
            ShapeType::ImportFile => {
                // ── Directory browser ──
                ui.horizontal(|ui| {
                    if ui.button("Set folder…").clicked() {
                        self.pick_import_directory();
                    }
                    if let Some(dir) = &self.import_dir_display {
                        ui.label(dir.as_str());
                    } else {
                        ui.label(
                            egui::RichText::new("No folder selected").weak(),
                        );
                    }
                });

                if self.import_directory.is_some() {
                    ui.horizontal(|ui| {
                        if ui.button("Refresh").clicked() {
                            self.scan_directory();
                        }
                        ui.label(
                            egui::RichText::new(format!(
                                "{} file(s)",
                                self.available_files.len()
                            ))
                            .weak()
                            .small(),
                        );
                    });

                    if !self.available_files.is_empty() {
                        let selected_text = self
                            .selected_import_file
                            .clone()
                            .unwrap_or_else(|| "Select a file…".into());

                        let mut new_selection: Option<ImportFileEntry> = None;

                        egui::ComboBox::from_label("File")
                            .selected_text(&selected_text)
                            .width(250.0)
                            .show_ui(ui, |ui| {
                                for entry in &self.available_files {
                                    let is_selected = self
                                        .selected_import_file
                                        .as_deref()
                                        == Some(entry.filename.as_str());
                                    if ui
                                        .selectable_label(is_selected, &entry.filename)
                                        .clicked()
                                    {
                                        new_selection = Some(entry.clone());
                                    }
                                }
                            });

                        if let Some(entry) = new_selection {
                            self.selected_import_file = Some(entry.filename.clone());
                            self.load_file_from_entry(&entry);
                        }
                    }
                }

                ui.add_space(4.0);

                // Display current loaded file info
                if let Some(path) = &self.obj_path_display {
                    ui.label(format!("Loaded OBJ: {}", path));
                    if let Some(mesh) = &self.obj_mesh {
                        ui.label(format!(
                            "  {} vertices, {} triangles",
                            mesh.vertices.len(),
                            mesh.faces.len()
                        ));
                    }
                } else if let Some(path) = &self.xyz_path_display {
                    ui.label(format!("Loaded XYZ: {}", path));
                    if let Some(positions) = &self.cached_positions {
                        ui.label(format!("  {} atoms", positions.len()));
                    }
                }

                if let Some(err) = &self.xyz_error {
                    ui.colored_label(egui::Color32::RED, format!("Parse error: {}", err));
                }
                if let Some(err) = &self.obj_error {
                    ui.colored_label(egui::Color32::RED, format!("Parse error: {}", err));
                }

                ui.add_space(4.0);
                ui.separator();
                ui.label(
                    egui::RichText::new("Or open individual file:")
                        .weak()
                        .small(),
                );

                ui.horizontal(|ui| {
                    if ui.button("Open .xyz file…").clicked() {
                        self.load_xyz_file();
                    }
                    if ui.button("Open .obj file…").clicked() {
                        self.load_obj_file();
                    }
                });

                ui.label(
                    egui::RichText::new(
                        ".xyz: Ångströms → nm.  .obj: coordinates in nm, mesh volume filled with dipoles.",
                    )
                    .weak()
                    .small(),
                );
            }
        }

        ui.add_space(8.0);
        let prev_spacing = self.dipole_spacing;
        ui.add(
            egui::Slider::new(&mut self.dipole_spacing, 0.5..=10.0).text("Dipole spacing (nm)"),
        );
        if (self.dipole_spacing - prev_spacing).abs() > 1e-6 {
            self.cache_dirty = true;
        }

        // Recompute lattice if dirty
        if self.cache_dirty {
            self.recompute_lattice();
        }

        ui.add_space(8.0);
        ui.label(format!("Dipoles: {}", self.dipole_count));

        // 2D scatter preview
        ui.add_space(8.0);
        ui.horizontal(|ui| {
            ui.label("Projection:");
            ui.selectable_value(&mut self.projection_plane, ProjectionPlane::XY, "XY");
            ui.selectable_value(&mut self.projection_plane, ProjectionPlane::XZ, "XZ");
            ui.selectable_value(&mut self.projection_plane, ProjectionPlane::YZ, "YZ");
        });

        if let Some(positions) = &self.cached_positions {
            let (ax_a, ax_b, label_a, label_b) = match self.projection_plane {
                ProjectionPlane::XY => (0, 1, "x (nm)", "y (nm)"),
                ProjectionPlane::XZ => (0, 2, "x (nm)", "z (nm)"),
                ProjectionPlane::YZ => (1, 2, "y (nm)", "z (nm)"),
            };

            let plot_points: egui_plot::PlotPoints = positions
                .iter()
                .map(|p| [p[ax_a], p[ax_b]])
                .collect();

            let points = egui_plot::Points::new(plot_points)
                .radius(2.5)
                .color(egui::Color32::from_rgb(80, 140, 220));

            egui_plot::Plot::new("dipole_scatter")
                .height(300.0)
                .data_aspect(1.0)
                .x_axis_label(label_a)
                .y_axis_label(label_b)
                .show(ui, |plot_ui| {
                    plot_ui.points(points);
                });
        }
    }

    /// Open a file dialog, load and parse a .xyz file, and store positions.
    fn load_xyz_file(&mut self) {
        use lumina_geometry::parsers::xyz::parse_xyz;

        let path = rfd::FileDialog::new()
            .set_title("Open XYZ file")
            .add_filter("XYZ", &["xyz"])
            .pick_file();

        if let Some(path) = path {
            self.xyz_error = None;
            match std::fs::read_to_string(&path) {
                Ok(content) => match parse_xyz(&content) {
                    Ok(points) => {
                        self.xyz_path_display =
                            Some(path.file_name().unwrap_or_default().to_string_lossy().into_owned());
                        self.cached_positions = Some(points.iter().map(|p| p.position).collect());
                        self.dipole_count = points.len();
                        self.cache_dirty = false;
                    }
                    Err(e) => {
                        self.xyz_error = Some(format!("{:?}", e));
                    }
                },
                Err(e) => {
                    self.xyz_error = Some(format!("Could not read file: {}", e));
                }
            }
        }
    }

    /// Open a file dialog, load and parse an .obj file, discretise, and store positions.
    fn load_obj_file(&mut self) {
        use lumina_geometry::parsers::obj::parse_obj;
        use lumina_geometry::discretise::discretise_mesh;

        let path = rfd::FileDialog::new()
            .set_title("Open OBJ file")
            .add_filter("OBJ", &["obj"])
            .pick_file();

        if let Some(path) = path {
            self.obj_error = None;
            self.xyz_error = None;
            match std::fs::read_to_string(&path) {
                Ok(content) => match parse_obj(&content) {
                    Ok(mesh) => {
                        self.obj_path_display =
                            Some(path.file_name().unwrap_or_default().to_string_lossy().into_owned());
                        self.xyz_path_display = None; // Clear xyz if loading obj
                        let lattice = discretise_mesh(&mesh, self.dipole_spacing);
                        self.dipole_count = lattice.len();
                        self.cached_positions = Some(lattice.iter().map(|p| p.position).collect());
                        self.obj_mesh = Some(mesh);
                        self.cache_dirty = false;
                    }
                    Err(e) => {
                        self.obj_error = Some(format!("{:?}", e));
                    }
                },
                Err(e) => {
                    self.obj_error = Some(format!("Could not read file: {}", e));
                }
            }
        }
    }

    /// Open a folder picker and scan the selected directory for structure files.
    fn pick_import_directory(&mut self) {
        let path = rfd::FileDialog::new()
            .set_title("Select structure folder")
            .pick_folder();

        if let Some(path) = path {
            self.import_dir_display = Some(path.display().to_string());
            self.import_directory = Some(path);
            self.scan_directory();
        }
    }

    /// Scan the import directory for .xyz and .obj files.
    fn scan_directory(&mut self) {
        self.available_files.clear();
        self.selected_import_file = None;

        let Some(dir) = &self.import_directory else {
            return;
        };

        let Ok(entries) = std::fs::read_dir(dir) else {
            return;
        };

        for entry in entries.flatten() {
            let path = entry.path();
            if !path.is_file() {
                continue;
            }
            let ext = path
                .extension()
                .and_then(|e| e.to_str())
                .map(|e| e.to_lowercase());
            match ext.as_deref() {
                Some("xyz") | Some("obj") => {
                    let filename = path
                        .file_name()
                        .unwrap_or_default()
                        .to_string_lossy()
                        .into_owned();
                    self.available_files.push(ImportFileEntry {
                        path,
                        filename,
                        extension: ext.unwrap(),
                    });
                }
                _ => {}
            }
        }

        self.available_files.sort_by(|a, b| a.filename.cmp(&b.filename));
    }

    /// Load and parse a file selected from the directory browser.
    fn load_file_from_entry(&mut self, entry: &ImportFileEntry) {
        self.xyz_error = None;
        self.obj_error = None;

        let content = match std::fs::read_to_string(&entry.path) {
            Ok(c) => c,
            Err(e) => {
                let msg = format!("Could not read file: {}", e);
                if entry.extension == "xyz" {
                    self.xyz_error = Some(msg);
                } else {
                    self.obj_error = Some(msg);
                }
                return;
            }
        };

        match entry.extension.as_str() {
            "xyz" => {
                use lumina_geometry::parsers::xyz::parse_xyz;
                match parse_xyz(&content) {
                    Ok(points) => {
                        self.xyz_path_display = Some(entry.filename.clone());
                        self.obj_path_display = None;
                        self.obj_mesh = None;
                        self.cached_positions =
                            Some(points.iter().map(|p| p.position).collect());
                        self.dipole_count = points.len();
                        self.cache_dirty = false;
                    }
                    Err(e) => {
                        self.xyz_error = Some(format!("{:?}", e));
                    }
                }
            }
            "obj" => {
                use lumina_geometry::parsers::obj::parse_obj;
                use lumina_geometry::discretise::discretise_mesh;
                match parse_obj(&content) {
                    Ok(mesh) => {
                        self.obj_path_display = Some(entry.filename.clone());
                        self.xyz_path_display = None;
                        let lattice = discretise_mesh(&mesh, self.dipole_spacing);
                        self.dipole_count = lattice.len();
                        self.cached_positions =
                            Some(lattice.iter().map(|p| p.position).collect());
                        self.obj_mesh = Some(mesh);
                        self.cache_dirty = false;
                    }
                    Err(e) => {
                        self.obj_error = Some(format!("{:?}", e));
                    }
                }
            }
            _ => {}
        }
    }

    fn recompute_lattice(&mut self) {
        use lumina_geometry::discretise::discretise_primitive;
        use lumina_geometry::primitives::*;

        let primitive = match self.selected_shape {
            ShapeType::Sphere => Primitive::Sphere(Sphere {
                centre: [0.0, 0.0, 0.0],
                radius: self.sphere_radius,
            }),
            ShapeType::Cylinder => Primitive::Cylinder(Cylinder {
                base_centre: [0.0, 0.0, -self.cylinder_length / 2.0],
                axis: [0.0, 0.0, 1.0],
                length: self.cylinder_length,
                radius: self.cylinder_radius,
            }),
            ShapeType::Cuboid => Primitive::Cuboid(Cuboid {
                centre: [0.0, 0.0, 0.0],
                half_extents: [self.cuboid_hx, self.cuboid_hy, self.cuboid_hz],
            }),
            ShapeType::Ellipsoid => Primitive::Ellipsoid(Ellipsoid {
                centre: [0.0, 0.0, 0.0],
                semi_axes: [self.ellipsoid_a, self.ellipsoid_b, self.ellipsoid_c],
            }),
            ShapeType::Helix => {
                use lumina_geometry::primitives::Helix;
                let total_height = self.helix_pitch * self.helix_turns;
                Primitive::Helix(Helix {
                    base_centre: [0.0, 0.0, -total_height / 2.0],
                    axis: [0.0, 0.0, 1.0],
                    radius: self.helix_radius,
                    pitch: self.helix_pitch,
                    turns: self.helix_turns,
                    wire_radius: self.helix_wire_radius,
                })
            }
            ShapeType::ImportFile => {
                if let Some(mesh) = &self.obj_mesh {
                    // Re-discretise the stored OBJ mesh at the current spacing
                    use lumina_geometry::discretise::discretise_mesh;
                    let lattice = discretise_mesh(mesh, self.dipole_spacing);
                    self.dipole_count = lattice.len();
                    self.cached_positions = Some(lattice.iter().map(|p| p.position).collect());
                } else {
                    // XYZ positions are fixed — no re-discretisation
                    self.dipole_count = self.cached_positions.as_ref().map(|v| v.len()).unwrap_or(0);
                }
                self.cache_dirty = false;
                return;
            }
        };

        let lattice = discretise_primitive(&primitive, self.dipole_spacing);
        self.dipole_count = lattice.len();
        self.cached_positions = Some(lattice.iter().map(|p| p.position).collect());
        self.cache_dirty = false;
    }
}

/// Helper: add a slider and return true if the value changed.
fn slider_changed(
    ui: &mut Ui,
    value: &mut f64,
    range: std::ops::RangeInclusive<f64>,
    text: &str,
) -> bool {
    let prev = *value;
    ui.add(egui::Slider::new(value, range).text(text));
    (*value - prev).abs() > 1e-6
}
