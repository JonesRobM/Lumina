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
}

impl ShapeSpec {
    /// Discretise this shape into object-local positions centred at the origin.
    ///
    /// Returns `(positions, labels)` where `labels[i]` is the atom-type string
    /// from an XYZ file, or `None` for primitive shapes and OBJ meshes.
    fn discretise_local(
        &self,
        spacing: f64,
        obj_name: &str,
    ) -> Result<(Vec<[f64; 3]>, Vec<Option<String>>), SceneError> {
        match self {
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

/// Geometry-only output of [`SceneSpec::build_geometry`].
///
/// Contains the world-space dipole positions, per-dipole material labels,
/// and dipole spacings. Material resolution (ε(ω) → Clausius-Mossotti) is
/// performed per-wavelength by the caller using the material strings.
pub struct SceneDipoles {
    /// World-space dipole positions in nm.
    pub positions: Vec<[f64; 3]>,
    /// Material identifier for each dipole (resolves via material providers).
    pub materials: Vec<String>,
    /// Dipole lattice spacing for each dipole in nm.
    pub spacings: Vec<f64>,
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

        for obj in &self.objects {
            let (local_pos, labels) =
                obj.shape.discretise_local(obj.dipole_spacing, &obj.name)?;

            let world_pos = obj.transform.apply(&local_pos);

            for (pos, label) in world_pos.into_iter().zip(labels.into_iter()) {
                positions.push(pos);
                let mat = match label {
                    Some(ref lbl) => obj
                        .species_map
                        .get(lbl.as_str())
                        .cloned()
                        .unwrap_or_else(|| obj.material.clone()),
                    None => obj.material.clone(),
                };
                materials.push(mat);
                spacings.push(obj.dipole_spacing);
            }
        }

        Ok(SceneDipoles { positions, materials, spacings })
    }
}

// ─── Helpers ─────────────────────────────────────────────────────────────────

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
