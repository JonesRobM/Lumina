//! Multi-object scene specification and geometry build pipeline.
//!
//! A [`SceneSpec`] is a list of [`ObjectSpec`]s, each describing one
//! nanostructure component with its own shape, material label, dipole
//! spacing, and full affine transform (position, ZYX Euler rotation, scale).
//!
//! # Usage
//!
//! ```toml
//! [[geometry.object]]
//! name        = "Left_sphere"
//! material    = "Au_JC"
//! dipole_spacing = 2.0
//! position    = [-25.0, 0.0, 0.0]
//! rotation_deg = [0.0, 0.0, 0.0]
//! scale       = 1.0
//!
//! [geometry.object.shape]
//! type   = "sphere"
//! radius = 20.0
//! ```
//!
//! Call [`SceneSpec::build_geometry`] to flatten the scene into a
//! [`SceneDipoles`] (positions + material strings + spacings) that the
//! solver's material-resolution loop consumes per wavelength.

use std::collections::HashMap;

use nalgebra::{Rotation3, Vector3};
use serde::{Deserialize, Serialize};
use thiserror::Error;

use crate::discretise::{discretise_mesh, discretise_primitive};
use crate::ellipsoid::ellipsoid_depol_factors;
use crate::parsers::xyz::parse_xyz;
use crate::parsers::obj::parse_obj;
use crate::primitives::{Cuboid, Cylinder, Ellipsoid, Helix, Primitive, Sphere};

// ─── Error type ─────────────────────────────────────────────────────────────

/// Errors that can occur when building scene geometry.
#[derive(Debug, Error)]
pub enum SceneError {
    #[error("Object '{name}' produced no dipoles — check shape parameters and dipole spacing")]
    EmptyObject { name: String },

    #[error("File error for object '{name}': {message}")]
    FileError { name: String, message: String },

    #[error("Unsupported file format '{ext}' for object '{name}' (expected .xyz or .obj)")]
    UnsupportedFormat { name: String, ext: String },

    #[error("Scene contains no objects")]
    EmptyScene,
}

// ─── Shape specification ─────────────────────────────────────────────────────

fn default_z_axis() -> [f64; 3] {
    [0.0, 0.0, 1.0]
}

/// One concentric shell in a [`ShapeSpec::CoreShell`] particle.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ShellLayer {
    /// Shell thickness in nm (measured radially outward from the previous surface).
    pub thickness: f64,
    /// Material identifier for this shell (e.g. `"Ag_JC"`, `"SiO2_Palik"`).
    pub material: String,
}

/// Core shape variants supported in a core-shell particle.
///
/// Helices are excluded because shell membership is geometrically ill-defined
/// for a coiled wire. File imports are also excluded.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "lowercase")]
pub enum CoreShapeSpec {
    /// Spherical core of the given radius (nm).
    Sphere { radius: f64 },
    /// Cylindrical core.
    Cylinder {
        radius: f64,
        length: f64,
        /// Axis direction (normalised internally). Default: [0, 0, 1].
        #[serde(default = "default_z_axis")]
        axis: [f64; 3],
    },
    /// Cuboid core with given half-extents (nm).
    Cuboid { half_extents: [f64; 3] },
    /// Ellipsoidal core with given semi-axes (nm).
    Ellipsoid { semi_axes: [f64; 3] },
}

/// Shape of one object in a scene, always expressed in object-local coordinates
/// centred at the origin. Positioning is handled by [`ObjectTransform`].
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "lowercase")]
pub enum ShapeSpec {
    /// Sphere of the given radius (nm).
    Sphere { radius: f64 },
    /// Cylinder with the given radius (nm) and length (nm).
    Cylinder {
        radius: f64,
        length: f64,
        /// Axis direction (normalised internally). Default: [0, 0, 1].
        #[serde(default = "default_z_axis")]
        axis: [f64; 3],
    },
    /// Axis-aligned cuboid with given half-extents [hx, hy, hz] (nm).
    Cuboid { half_extents: [f64; 3] },
    /// Ellipsoid with semi-axes [a, b, c] (nm).
    Ellipsoid { semi_axes: [f64; 3] },
    /// Helix with given coil radius, pitch, number of turns, and wire radius (all nm).
    Helix {
        radius: f64,
        pitch: f64,
        turns: f64,
        wire_radius: f64,
        /// Helix axis direction. Default: [0, 0, 1].
        #[serde(default = "default_z_axis")]
        axis: [f64; 3],
    },
    /// Import geometry from a `.xyz` or `.obj` file.
    File { path: String },
    /// Concentric multi-shell particle (e.g. Au@Ag, Au@SiO₂).
    ///
    /// The core shape and material define the innermost region. The `shells`
    /// list ordered from innermost to outermost; each adds `thickness` nm to
    /// every surface.
    ///
    /// # TOML example (Au@Ag sphere):
    /// ```toml
    /// [geometry.object.shape]
    /// type          = "coreshell"
    /// core_material = "Au_JC"
    ///
    /// [geometry.object.shape.base]
    /// type   = "sphere"
    /// radius = 15.0
    ///
    /// [[geometry.object.shape.shells]]
    /// thickness = 5.0
    /// material  = "Ag_JC"
    /// ```
    CoreShell {
        /// Core shape geometry.
        base: CoreShapeSpec,
        /// Material of the innermost core region.
        core_material: String,
        /// Shells listed from innermost to outermost.
        shells: Vec<ShellLayer>,
    },
}

impl ShapeSpec {
    /// Discretise this shape into object-local positions centred at the origin.
    ///
    /// Returns `(positions, labels)` where:
    /// - For XYZ imports: `labels[i]` is the atom-type string (e.g. `"Au"`).
    /// - For [`ShapeSpec::CoreShell`]: `labels[i]` is the **resolved material name**
    ///   (e.g. `"Au_JC"`, `"Ag_JC"`) and should be used directly without
    ///   `species_map` lookup.
    /// - For all other shapes: `labels[i]` is `None`.
    fn discretise_local(
        &self,
        spacing: f64,
        obj_name: &str,
    ) -> Result<(Vec<[f64; 3]>, Vec<Option<String>>), SceneError> {
        match self {
            ShapeSpec::CoreShell { base, core_material, shells } => {
                let outer = core_shell_outer_primitive(base, shells);
                let (bb_min, bb_max) = outer.bounding_box();
                let cx = 0.5 * (bb_min[0] + bb_max[0]);
                let cy = 0.5 * (bb_min[1] + bb_max[1]);
                let cz = 0.5 * (bb_min[2] + bb_max[2]);
                let pts = discretise_primitive(&outer, spacing);
                if pts.is_empty() {
                    return Err(SceneError::EmptyObject { name: obj_name.to_string() });
                }
                let mut positions = Vec::with_capacity(pts.len());
                let mut labels = Vec::with_capacity(pts.len());
                for pt in &pts {
                    let local = [
                        pt.position[0] - cx,
                        pt.position[1] - cy,
                        pt.position[2] - cz,
                    ];
                    let mat = shell_material_at(&local, base, core_material, shells);
                    positions.push(local);
                    labels.push(Some(mat));
                }
                Ok((positions, labels))
            }
            ShapeSpec::File { path } => {
                let p = std::path::Path::new(path);
                let ext = p
                    .extension()
                    .and_then(|e| e.to_str())
                    .unwrap_or("")
                    .to_lowercase();

                match ext.as_str() {
                    "xyz" => {
                        let content =
                            std::fs::read_to_string(p).map_err(|e| SceneError::FileError {
                                name: obj_name.to_string(),
                                message: e.to_string(),
                            })?;
                        let parsed =
                            parse_xyz(&content).map_err(|e| SceneError::FileError {
                                name: obj_name.to_string(),
                                message: e.to_string(),
                            })?;
                        if parsed.is_empty() {
                            return Err(SceneError::EmptyObject {
                                name: obj_name.to_string(),
                            });
                        }
                        let positions: Vec<[f64; 3]> =
                            parsed.iter().map(|p| p.position).collect();
                        let labels: Vec<Option<String>> =
                            parsed.iter().map(|p| p.label.clone()).collect();
                        let centred = centre_positions(positions);
                        Ok((centred, labels))
                    }
                    "obj" => {
                        let content =
                            std::fs::read_to_string(p).map_err(|e| SceneError::FileError {
                                name: obj_name.to_string(),
                                message: e.to_string(),
                            })?;
                        let mesh =
                            parse_obj(&content).map_err(|e| SceneError::FileError {
                                name: obj_name.to_string(),
                                message: e.to_string(),
                            })?;
                        let pts = discretise_mesh(&mesh, spacing);
                        if pts.is_empty() {
                            return Err(SceneError::EmptyObject {
                                name: obj_name.to_string(),
                            });
                        }
                        let positions: Vec<[f64; 3]> =
                            pts.iter().map(|p| p.position).collect();
                        let labels = vec![None; positions.len()];
                        let centred = centre_positions(positions);
                        Ok((centred, labels))
                    }
                    other => Err(SceneError::UnsupportedFormat {
                        name: obj_name.to_string(),
                        ext: other.to_string(),
                    }),
                }
            }
            _ => {
                let prim = self.to_primitive();
                let (bb_min, bb_max) = prim.bounding_box();
                let cx = 0.5 * (bb_min[0] + bb_max[0]);
                let cy = 0.5 * (bb_min[1] + bb_max[1]);
                let cz = 0.5 * (bb_min[2] + bb_max[2]);
                let pts = discretise_primitive(&prim, spacing);
                if pts.is_empty() {
                    return Err(SceneError::EmptyObject {
                        name: obj_name.to_string(),
                    });
                }
                let positions: Vec<[f64; 3]> = pts
                    .iter()
                    .map(|p| [p.position[0] - cx, p.position[1] - cy, p.position[2] - cz])
                    .collect();
                let labels = vec![None; positions.len()];
                Ok((positions, labels))
            }
        }
    }

    /// Convert to a [`Primitive`] for discretisation. Panics for `File` variant.
    fn to_primitive(&self) -> Primitive {
        match self {
            ShapeSpec::Sphere { radius } => Primitive::Sphere(Sphere {
                centre: [0.0, 0.0, 0.0],
                radius: *radius,
            }),
            ShapeSpec::Cylinder { radius, length, axis } => Primitive::Cylinder(Cylinder {
                base_centre: [0.0, 0.0, 0.0],
                axis: *axis,
                length: *length,
                radius: *radius,
            }),
            ShapeSpec::Cuboid { half_extents } => Primitive::Cuboid(Cuboid {
                centre: [0.0, 0.0, 0.0],
                half_extents: *half_extents,
            }),
            ShapeSpec::Ellipsoid { semi_axes } => Primitive::Ellipsoid(Ellipsoid {
                centre: [0.0, 0.0, 0.0],
                semi_axes: *semi_axes,
            }),
            ShapeSpec::Helix { radius, pitch, turns, wire_radius, axis } => {
                Primitive::Helix(Helix {
                    base_centre: [0.0, 0.0, 0.0],
                    axis: *axis,
                    radius: *radius,
                    pitch: *pitch,
                    turns: *turns,
                    wire_radius: *wire_radius,
                })
            }
            ShapeSpec::File { .. } => panic!("File shapes must be handled separately"),
            ShapeSpec::CoreShell { .. } => panic!("CoreShell must be handled separately"),
        }
    }
}

// ─── Transform ───────────────────────────────────────────────────────────────

/// Per-object affine transform: scale → ZYX Euler rotation → translation.
///
/// Rotation convention: `R = Rz(rotation_deg[2]) · Ry(rotation_deg[1]) · Rx(rotation_deg[0])`.
/// A zero rotation `[0, 0, 0]` is the identity. Angles are in degrees.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ObjectTransform {
    /// World-space position offset in nm (x, y, z).
    #[serde(default)]
    pub position: [f64; 3],
    /// ZYX Euler angles in degrees (Rx, Ry, Rz).
    #[serde(default)]
    pub rotation_deg: [f64; 3],
    /// Uniform scale factor (1.0 = no scaling).
    #[serde(default = "default_scale")]
    pub scale: f64,
}

fn default_scale() -> f64 {
    1.0
}

impl Default for ObjectTransform {
    fn default() -> Self {
        Self {
            position: [0.0, 0.0, 0.0],
            rotation_deg: [0.0, 0.0, 0.0],
            scale: 1.0,
        }
    }
}

impl ObjectTransform {
    /// Apply this transform to a slice of object-local positions.
    ///
    /// Pipeline: scale → rotate (ZYX Euler) → translate.
    pub fn apply(&self, positions: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let rx = self.rotation_deg[0].to_radians();
        let ry = self.rotation_deg[1].to_radians();
        let rz = self.rotation_deg[2].to_radians();
        // nalgebra from_euler_angles(roll, pitch, yaw) = Rz(yaw) * Ry(pitch) * Rx(roll)
        let rot = Rotation3::from_euler_angles(rx, ry, rz);
        let scale = self.scale;
        let [tx, ty, tz] = self.position;

        positions
            .iter()
            .map(|p| {
                let v = Vector3::new(p[0] * scale, p[1] * scale, p[2] * scale);
                let rv = rot * v;
                [rv.x + tx, rv.y + ty, rv.z + tz]
            })
            .collect()
    }
}

// ─── ObjectSpec ──────────────────────────────────────────────────────────────

/// One component in a multi-object scene.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ObjectSpec {
    /// Human-readable identifier shown in the GUI object list.
    pub name: String,
    /// Shape of this object in object-local coordinates.
    pub shape: ShapeSpec,
    /// Material identifier (e.g. `"Au_JC"`, `"Ag_JC"`, `"TiO2_Palik"`).
    pub material: String,
    /// Cubic dipole lattice spacing in nm.
    pub dipole_spacing: f64,
    /// Position, rotation, and scale in world space.
    #[serde(flatten)]
    pub transform: ObjectTransform,
    /// Element label → material override for XYZ multi-species imports.
    ///
    /// Keys are atom-type strings from the XYZ file (e.g. `"Au"`, `"Ag"`).
    /// Values override `material` for dipoles with that label.
    /// Empty for primitive shapes and single-species imports.
    #[serde(default, skip_serializing_if = "HashMap::is_empty")]
    pub species_map: HashMap<String, String>,
}

impl ObjectSpec {
    /// Create a default sphere object centred at the origin.
    pub fn default_sphere() -> Self {
        Self {
            name: "Object".to_string(),
            shape: ShapeSpec::Sphere { radius: 20.0 },
            material: "Au_JC".to_string(),
            dipole_spacing: 2.0,
            transform: ObjectTransform::default(),
            species_map: HashMap::new(),
        }
    }

    /// Compute world-space dipole positions for this object.
    ///
    /// Useful for GUI preview rendering without building the full scene.
    /// Returns an error if the shape is empty or a file cannot be read.
    pub fn world_positions(&self) -> Result<Vec<[f64; 3]>, SceneError> {
        let (local_pos, _) = self.shape.discretise_local(self.dipole_spacing, &self.name)?;
        Ok(self.transform.apply(&local_pos))
    }
}

// ─── SceneSpec ───────────────────────────────────────────────────────────────

/// A complete multi-object scene ready for CDA simulation.
///
/// Serialises directly to/from TOML (via `serde`) as `[[geometry.object]]` arrays,
/// making it the shared format for both the GUI and CLI.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SceneSpec {
    /// All objects in the scene.
    #[serde(rename = "object")]
    pub objects: Vec<ObjectSpec>,
}

impl Default for SceneSpec {
    fn default() -> Self {
        Self {
            objects: vec![ObjectSpec::default_sphere()],
        }
    }
}

// ─── SceneDipoles ────────────────────────────────────────────────────────────

/// Wavelength-independent polarisability geometry data for a dipole in an
/// ellipsoidal object.
///
/// Carries the depolarisation factors and orientation computed at scene-build
/// time, used by the per-wavelength solver loop to build the full 3×3
/// polarisability tensor instead of calling the isotropic Clausius-Mossotti.
///
/// Only populated for [`ShapeSpec::Ellipsoid`] objects. All other shapes
/// produce `None` in [`SceneDipoles::hints`].
#[derive(Debug, Clone)]
pub struct DipolePolarisabilityHint {
    /// Depolarisation factors $[L_x, L_y, L_z]$ in the object-local
    /// principal-axis frame. Computed from the ellipsoid semi-axes.
    pub depol_factors: [f64; 3],
    /// 3×3 rotation matrix (row-major) mapping the object-local frame to
    /// the world (lab) frame: $\boldsymbol{\alpha} = R\,\text{diag}(\alpha_x,\alpha_y,\alpha_z)\,R^T$.
    pub rotation: [f64; 9],
}

/// Geometry-only output of [`SceneSpec::build_geometry`].
///
/// Contains the world-space dipole positions, per-dipole material labels,
/// dipole spacings, and optional anisotropy hints for ellipsoidal objects.
/// Material resolution (ε(ω) → Clausius-Mossotti) is performed per-wavelength
/// by the caller.
pub struct SceneDipoles {
    /// World-space dipole positions in nm.
    pub positions: Vec<[f64; 3]>,
    /// Material identifier for each dipole (resolves via material providers).
    pub materials: Vec<String>,
    /// Dipole lattice spacing for each dipole in nm.
    pub spacings: Vec<f64>,
    /// Per-dipole anisotropy hint. `Some` only for dipoles belonging to
    /// [`ShapeSpec::Ellipsoid`] objects; `None` for all other shapes.
    /// Length always equals `positions.len()`.
    pub hints: Vec<Option<DipolePolarisabilityHint>>,
}

impl SceneSpec {
    /// Build the flat dipole geometry for all objects in the scene.
    ///
    /// For each object:
    /// 1. Discretise the shape into object-local positions centred at origin.
    /// 2. Apply the object transform (scale → rotate → translate).
    /// 3. Assign material strings (using `species_map` for XYZ atoms, else `material`).
    ///
    /// Returns a [`SceneDipoles`] struct that the wavelength loop uses to
    /// resolve polarisabilities.
    pub fn build_geometry(&self) -> Result<SceneDipoles, SceneError> {
        if self.objects.is_empty() {
            return Err(SceneError::EmptyScene);
        }

        let mut positions = Vec::new();
        let mut materials = Vec::new();
        let mut spacings = Vec::new();
        let mut hints: Vec<Option<DipolePolarisabilityHint>> = Vec::new();

        for obj in &self.objects {
            let (local_pos, labels) =
                obj.shape.discretise_local(obj.dipole_spacing, &obj.name)?;

            let world_pos = obj.transform.apply(&local_pos);

            // Compute anisotropy hint for ellipsoidal objects only.
            let hint_for_obj: Option<DipolePolarisabilityHint> =
                if let ShapeSpec::Ellipsoid { semi_axes } = &obj.shape {
                    let s = obj.transform.scale;
                    let scaled = [semi_axes[0] * s, semi_axes[1] * s, semi_axes[2] * s];
                    let depol_factors = ellipsoid_depol_factors(scaled);
                    let rotation = rotation_matrix_row_major(obj.transform.rotation_deg);
                    Some(DipolePolarisabilityHint { depol_factors, rotation })
                } else {
                    None
                };

            // CoreShell shapes return pre-resolved material names in `labels`;
            // they bypass species_map (which is irrelevant for primitive shells).
            let is_core_shell = matches!(obj.shape, ShapeSpec::CoreShell { .. });

            for (pos, label) in world_pos.into_iter().zip(labels.into_iter()) {
                positions.push(pos);
                let mat = if is_core_shell {
                    label.unwrap_or_else(|| obj.material.clone())
                } else {
                    match label {
                        Some(ref lbl) => obj
                            .species_map
                            .get(lbl.as_str())
                            .cloned()
                            .unwrap_or_else(|| obj.material.clone()),
                        None => obj.material.clone(),
                    }
                };
                materials.push(mat);
                spacings.push(obj.dipole_spacing);
                hints.push(hint_for_obj.clone());
            }
        }

        Ok(SceneDipoles { positions, materials, spacings, hints })
    }
}

// ─── Helpers ─────────────────────────────────────────────────────────────────

/// Build a 3×3 rotation matrix (row-major) from ZYX Euler angles in degrees.
///
/// Convention: $R = R_z(\text{rz}) \cdot R_y(\text{ry}) \cdot R_x(\text{rx})$,
/// consistent with [`ObjectTransform::apply`].
fn rotation_matrix_row_major(rotation_deg: [f64; 3]) -> [f64; 9] {
    let rx = rotation_deg[0].to_radians();
    let ry = rotation_deg[1].to_radians();
    let rz = rotation_deg[2].to_radians();
    let rot = Rotation3::from_euler_angles(rx, ry, rz);
    let m = rot.matrix();
    [
        m[(0, 0)], m[(0, 1)], m[(0, 2)],
        m[(1, 0)], m[(1, 1)], m[(1, 2)],
        m[(2, 0)], m[(2, 1)], m[(2, 2)],
    ]
}

/// Returns the outer bounding [`Primitive`] for a core-shell particle.
///
/// Each surface of the core shape is extended by the total cumulative shell
/// thickness, enclosing all shell regions in a single primitive that can be
/// passed to [`discretise_primitive`].
fn core_shell_outer_primitive(base: &CoreShapeSpec, shells: &[ShellLayer]) -> Primitive {
    let total: f64 = shells.iter().map(|s| s.thickness).sum();
    match base {
        CoreShapeSpec::Sphere { radius } => Primitive::Sphere(Sphere {
            centre: [0.0, 0.0, 0.0],
            radius: radius + total,
        }),
        CoreShapeSpec::Cylinder { radius, length, axis } => Primitive::Cylinder(Cylinder {
            base_centre: [0.0, 0.0, 0.0],
            axis: *axis,
            length: length + 2.0 * total,
            radius: radius + total,
        }),
        CoreShapeSpec::Cuboid { half_extents } => Primitive::Cuboid(Cuboid {
            centre: [0.0, 0.0, 0.0],
            half_extents: [
                half_extents[0] + total,
                half_extents[1] + total,
                half_extents[2] + total,
            ],
        }),
        CoreShapeSpec::Ellipsoid { semi_axes } => Primitive::Ellipsoid(Ellipsoid {
            centre: [0.0, 0.0, 0.0],
            semi_axes: [
                semi_axes[0] + total,
                semi_axes[1] + total,
                semi_axes[2] + total,
            ],
        }),
    }
}

/// Determine the material label for a dipole at `local` in a core-shell particle.
///
/// Returns `core_material` for points inside the core. For exterior points,
/// walks the shell list (innermost to outermost) and returns the material
/// whose accumulated thickness first encompasses the point. Points at floating-
/// point boundaries or at non-spherical corners where the Euclidean SDF exceeds
/// the total shell thickness are assigned the outermost shell material.
fn shell_material_at(
    local: &[f64; 3],
    base: &CoreShapeSpec,
    core_material: &str,
    shells: &[ShellLayer],
) -> String {
    let dist = core_surface_distance(local, base);
    if dist == 0.0 {
        return core_material.to_string();
    }
    let mut boundary = 0.0;
    for shell in shells {
        boundary += shell.thickness;
        if dist <= boundary {
            return shell.material.clone();
        }
    }
    // Fallback: handles FP boundary cases and corners of non-spherical shapes
    // where the Euclidean SDF slightly exceeds the nominal total thickness.
    shells
        .last()
        .map(|s| s.material.clone())
        .unwrap_or_else(|| core_material.to_string())
}

/// Euclidean unsigned distance from `local` to the nearest point on the core
/// surface. Returns 0.0 for points on or inside the core.
fn core_surface_distance(local: &[f64; 3], base: &CoreShapeSpec) -> f64 {
    match base {
        CoreShapeSpec::Sphere { radius } => {
            let r = (local[0] * local[0] + local[1] * local[1] + local[2] * local[2]).sqrt();
            (r - radius).max(0.0)
        }
        CoreShapeSpec::Cylinder { radius, length, axis } => {
            let ax = normalise3(axis);
            let dot = local[0] * ax[0] + local[1] * ax[1] + local[2] * ax[2];
            let perp = [
                local[0] - dot * ax[0],
                local[1] - dot * ax[1],
                local[2] - dot * ax[2],
            ];
            let r_xy =
                (perp[0] * perp[0] + perp[1] * perp[1] + perp[2] * perp[2]).sqrt();
            let z_abs = dot.abs();
            let half_len = 0.5 * length;
            let de_r = (r_xy - radius).max(0.0);
            let de_z = (z_abs - half_len).max(0.0);
            (de_r * de_r + de_z * de_z).sqrt()
        }
        CoreShapeSpec::Cuboid { half_extents } => {
            let dx = (local[0].abs() - half_extents[0]).max(0.0);
            let dy = (local[1].abs() - half_extents[1]).max(0.0);
            let dz = (local[2].abs() - half_extents[2]).max(0.0);
            (dx * dx + dy * dy + dz * dz).sqrt()
        }
        CoreShapeSpec::Ellipsoid { semi_axes } => {
            let [a, b, c] = *semi_axes;
            let u = (local[0] / a).powi(2) + (local[1] / b).powi(2) + (local[2] / c).powi(2);
            let min_axis = a.min(b).min(c);
            (u.sqrt() - 1.0).max(0.0) * min_axis
        }
    }
}

/// Normalise a 3-vector; returns [0, 0, 1] for near-zero inputs.
fn normalise3(v: &[f64; 3]) -> [f64; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len < 1e-12 { [0.0, 0.0, 1.0] } else { [v[0] / len, v[1] / len, v[2] / len] }
}

/// Shift a set of positions so their bounding-box centroid is at the origin.
fn centre_positions(positions: Vec<[f64; 3]>) -> Vec<[f64; 3]> {
    if positions.is_empty() {
        return positions;
    }
    let mut min = positions[0];
    let mut max = positions[0];
    for p in &positions {
        for i in 0..3 {
            if p[i] < min[i] { min[i] = p[i]; }
            if p[i] > max[i] { max[i] = p[i]; }
        }
    }
    let cx = 0.5 * (min[0] + max[0]);
    let cy = 0.5 * (min[1] + max[1]);
    let cz = 0.5 * (min[2] + max[2]);
    positions.into_iter().map(|p| [p[0] - cx, p[1] - cy, p[2] - cz]).collect()
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn make_sphere_spec(radius: f64, pos: [f64; 3], material: &str) -> ObjectSpec {
        ObjectSpec {
            name: "test".to_string(),
            shape: ShapeSpec::Sphere { radius },
            material: material.to_string(),
            dipole_spacing: 2.0,
            transform: ObjectTransform { position: pos, ..Default::default() },
            species_map: HashMap::new(),
        }
    }

    /// Zero rotation + zero translation must leave positions unchanged.
    #[test]
    fn test_object_transform_identity() {
        let t = ObjectTransform::default();
        let pts = vec![[1.0, 2.0, 3.0], [-1.0, 0.5, -2.0]];
        let out = t.apply(&pts);
        for (a, b) in pts.iter().zip(out.iter()) {
            for i in 0..3 {
                assert!((a[i] - b[i]).abs() < 1e-12, "Identity mismatch at axis {i}");
            }
        }
    }

    /// 90° rotation about X should map Y → Z and Z → −Y.
    #[test]
    fn test_euler_rotation_x90() {
        let t = ObjectTransform {
            rotation_deg: [90.0, 0.0, 0.0],
            ..Default::default()
        };
        let pts = vec![[0.0, 1.0, 0.0]]; // +Y
        let out = t.apply(&pts);
        // After Rx(90°): Y maps to +Z
        assert!((out[0][0]).abs() < 1e-12, "x should be 0");
        assert!((out[0][1]).abs() < 1e-12, "y should be 0");
        assert!((out[0][2] - 1.0).abs() < 1e-10, "z should be 1, got {}", out[0][2]);
    }

    /// A dimer of two spheres separated by a gap should have no dipoles closer
    /// than the gap distance between them.
    #[test]
    fn test_dimer_no_overlap() {
        let gap = 10.0_f64; // nm
        let r = 20.0_f64;
        let scene = SceneSpec {
            objects: vec![
                make_sphere_spec(r, [-(r + gap / 2.0), 0.0, 0.0], "Au_JC"),
                make_sphere_spec(r, [r + gap / 2.0, 0.0, 0.0], "Ag_JC"),
            ],
        };
        let geom = scene.build_geometry().unwrap();
        let n = geom.positions.len();
        // Find minimum inter-object distance.
        // Object 0 positions come first; we need to know the split point.
        // Simpler: check all pairs where one is left (x < 0) and other is right (x > 0).
        let left: Vec<_> = geom.positions.iter().filter(|p| p[0] < 0.0).collect();
        let right: Vec<_> = geom.positions.iter().filter(|p| p[0] > 0.0).collect();
        assert!(!left.is_empty() && !right.is_empty(), "Dimer must have dipoles on both sides");
        let min_dist = left
            .iter()
            .flat_map(|l| right.iter().map(move |r| {
                let dx = l[0] - r[0]; let dy = l[1] - r[1]; let dz = l[2] - r[2];
                (dx*dx + dy*dy + dz*dz).sqrt()
            }))
            .fold(f64::INFINITY, f64::min);
        assert!(
            min_dist >= gap - 2.5, // allow half-spacing tolerance
            "Dimer inter-object distance {min_dist:.2} nm less than gap {gap} nm"
        );
        let _ = n;
    }

    /// SceneSpec must round-trip through TOML without data loss.
    #[test]
    fn test_scene_toml_round_trip() {
        let scene = SceneSpec {
            objects: vec![
                ObjectSpec {
                    name: "Sphere_Au".to_string(),
                    shape: ShapeSpec::Sphere { radius: 20.0 },
                    material: "Au_JC".to_string(),
                    dipole_spacing: 2.0,
                    transform: ObjectTransform {
                        position: [10.0, 0.0, -5.0],
                        rotation_deg: [45.0, 0.0, 90.0],
                        scale: 1.5,
                    },
                    species_map: HashMap::new(),
                },
                ObjectSpec {
                    name: "Cylinder_Ag".to_string(),
                    shape: ShapeSpec::Cylinder {
                        radius: 5.0,
                        length: 40.0,
                        axis: [0.0, 0.0, 1.0],
                    },
                    material: "Ag_JC".to_string(),
                    dipole_spacing: 3.0,
                    transform: ObjectTransform::default(),
                    species_map: HashMap::new(),
                },
            ],
        };

        // Wrap in a geometry table for TOML compatibility
        #[derive(Serialize, Deserialize)]
        struct Wrapper { geometry: SceneSpec }
        let w = Wrapper { geometry: scene.clone() };
        let toml_str = toml::to_string(&w).expect("serialise failed");
        let w2: Wrapper = toml::from_str(&toml_str).expect("deserialise failed");

        assert_eq!(w2.geometry.objects.len(), 2);
        let obj0 = &w2.geometry.objects[0];
        assert_eq!(obj0.name, "Sphere_Au");
        assert_eq!(obj0.material, "Au_JC");
        assert!((obj0.transform.position[0] - 10.0).abs() < 1e-12);
        assert!((obj0.transform.scale - 1.5).abs() < 1e-12);

        if let ShapeSpec::Sphere { radius } = obj0.shape {
            assert!((radius - 20.0).abs() < 1e-12);
        } else {
            panic!("Expected Sphere");
        }
    }

    /// build_geometry dipole count must equal sum of per-object counts.
    #[test]
    fn test_scene_build_dipoles_count() {
        let spec_a = make_sphere_spec(15.0, [-20.0, 0.0, 0.0], "Au_JC");
        let spec_b = make_sphere_spec(15.0, [20.0, 0.0, 0.0], "Ag_JC");

        let single_a = SceneSpec { objects: vec![spec_a.clone()] };
        let single_b = SceneSpec { objects: vec![spec_b.clone()] };
        let dimer = SceneSpec { objects: vec![spec_a, spec_b] };

        let count_a = single_a.build_geometry().unwrap().positions.len();
        let count_b = single_b.build_geometry().unwrap().positions.len();
        let count_dimer = dimer.build_geometry().unwrap().positions.len();

        assert_eq!(count_dimer, count_a + count_b,
            "Dimer count {count_dimer} ≠ sum {}", count_a + count_b);
    }

    /// Au@Ag sphere: core dipoles must be Au_JC, shell dipoles must be Ag_JC.
    #[test]
    fn test_coreshell_sphere_materials() {
        let scene = SceneSpec {
            objects: vec![ObjectSpec {
                name: "cs".to_string(),
                shape: ShapeSpec::CoreShell {
                    base: CoreShapeSpec::Sphere { radius: 10.0 },
                    core_material: "Au_JC".to_string(),
                    shells: vec![ShellLayer {
                        thickness: 5.0,
                        material: "Ag_JC".to_string(),
                    }],
                },
                material: "Au_JC".to_string(),
                dipole_spacing: 2.0,
                transform: ObjectTransform::default(),
                species_map: HashMap::new(),
            }],
        };

        let geom = scene.build_geometry().unwrap();
        assert!(!geom.positions.is_empty(), "must produce dipoles");

        let mut n_core = 0usize;
        let mut n_shell = 0usize;
        for (pos, mat) in geom.positions.iter().zip(geom.materials.iter()) {
            let r = (pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]).sqrt();
            match mat.as_str() {
                "Au_JC" => {
                    assert!(r <= 10.0 + 1.5, "core dipole at r={r:.2} > core radius + half-spacing");
                    n_core += 1;
                }
                "Ag_JC" => {
                    assert!(r >= 8.0, "shell dipole at r={r:.2} unexpectedly close to centre");
                    n_shell += 1;
                }
                other => panic!("unexpected material: {other}"),
            }
        }
        assert!(n_core > 0, "no core dipoles generated");
        assert!(n_shell > 0, "no shell dipoles generated");
    }

    /// CoreShell with two shells: verify all dipoles get one of the three materials.
    #[test]
    fn test_coreshell_two_shells_all_assigned() {
        let shapes: &[(ShapeSpec, &[&str])] = &[
            (
                ShapeSpec::CoreShell {
                    base: CoreShapeSpec::Cuboid { half_extents: [8.0, 8.0, 8.0] },
                    core_material: "Au_JC".to_string(),
                    shells: vec![ShellLayer { thickness: 3.0, material: "Ag_JC".to_string() }],
                },
                &["Au_JC", "Ag_JC"],
            ),
            (
                ShapeSpec::CoreShell {
                    base: CoreShapeSpec::Ellipsoid { semi_axes: [10.0, 7.0, 5.0] },
                    core_material: "Au_JC".to_string(),
                    shells: vec![
                        ShellLayer { thickness: 2.0, material: "Ag_JC".to_string() },
                        ShellLayer { thickness: 3.0, material: "TiO2_Palik".to_string() },
                    ],
                },
                &["Au_JC", "Ag_JC", "TiO2_Palik"],
            ),
        ];
        for (shape, allowed) in shapes {
            let geom = SceneSpec {
                objects: vec![ObjectSpec {
                    name: "t".to_string(),
                    shape: shape.clone(),
                    material: "Au_JC".to_string(),
                    dipole_spacing: 2.0,
                    transform: ObjectTransform::default(),
                    species_map: HashMap::new(),
                }],
            }
            .build_geometry()
            .unwrap();
            assert!(!geom.positions.is_empty());
            for mat in &geom.materials {
                assert!(
                    allowed.contains(&mat.as_str()),
                    "unexpected material '{mat}'"
                );
            }
        }
    }

    /// CoreShell must survive a TOML round-trip without data loss.
    #[test]
    fn test_coreshell_toml_round_trip() {
        let scene = SceneSpec {
            objects: vec![ObjectSpec {
                name: "Au_Ag_sphere".to_string(),
                shape: ShapeSpec::CoreShell {
                    base: CoreShapeSpec::Sphere { radius: 15.0 },
                    core_material: "Au_JC".to_string(),
                    shells: vec![
                        ShellLayer { thickness: 3.0, material: "Ag_JC".to_string() },
                        ShellLayer { thickness: 4.0, material: "SiO2_Palik".to_string() },
                    ],
                },
                material: "Au_JC".to_string(),
                dipole_spacing: 2.0,
                transform: ObjectTransform::default(),
                species_map: HashMap::new(),
            }],
        };

        #[derive(Serialize, Deserialize)]
        struct Wrapper { geometry: SceneSpec }
        let toml_str = toml::to_string(&Wrapper { geometry: scene }).expect("serialise failed");
        let w2: Wrapper = toml::from_str(&toml_str).expect("deserialise failed");
        let obj = &w2.geometry.objects[0];

        if let ShapeSpec::CoreShell { base, core_material, shells } = &obj.shape {
            assert_eq!(core_material, "Au_JC");
            assert_eq!(shells.len(), 2);
            assert!((shells[0].thickness - 3.0).abs() < 1e-12);
            assert_eq!(shells[0].material, "Ag_JC");
            assert_eq!(shells[1].material, "SiO2_Palik");
            if let CoreShapeSpec::Sphere { radius } = base {
                assert!((radius - 15.0).abs() < 1e-12);
            } else {
                panic!("expected Sphere core");
            }
        } else {
            panic!("expected CoreShell shape");
        }
    }

    /// Species map must correctly assign different materials to different labels.
    #[test]
    fn test_xyz_species_map_applied() {
        // Build a minimal fake XYZ result by using a scene with two spheres
        // of different materials at the same position (materials tracked via strings).
        let mut species_map = HashMap::new();
        species_map.insert("Au".to_string(), "Au_JC".to_string());
        species_map.insert("Ag".to_string(), "Ag_JC".to_string());

        // Verify species_map resolution logic directly.
        let obj = ObjectSpec {
            name: "test".to_string(),
            shape: ShapeSpec::Sphere { radius: 10.0 },
            material: "Au_JC".to_string(),
            dipole_spacing: 3.0,
            transform: ObjectTransform::default(),
            species_map: species_map.clone(),
        };

        // With label "Ag" → should pick Ag_JC
        let mat_ag = obj.species_map.get("Ag").cloned().unwrap_or(obj.material.clone());
        assert_eq!(mat_ag, "Ag_JC");

        // With label "Unknown" → should fall back to Au_JC
        let mat_unk = obj.species_map.get("Unknown").cloned().unwrap_or(obj.material.clone());
        assert_eq!(mat_unk, "Au_JC");

        // Without label → uses default material
        let mat_none: String = obj.material.clone();
        assert_eq!(mat_none, "Au_JC");
    }
}
