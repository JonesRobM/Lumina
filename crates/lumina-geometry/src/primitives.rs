//! Parametric geometric primitives.
//!
//! Each primitive defines a closed volume in 3D space that can be discretised
//! into a dipole lattice by the [`discretise`](crate::discretise) module.
//! Primitives are fully described by their TOML parameters and can be
//! rescaled interactively in the GUI.

use serde::{Deserialize, Serialize};

/// A geometric shape that can be discretised into dipoles.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum Primitive {
    Sphere(Sphere),
    Cylinder(Cylinder),
    Cuboid(Cuboid),
    Helix(Helix),
    Ellipsoid(Ellipsoid),
}

/// A sphere defined by its centre and radius.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Sphere {
    /// Centre position (nm).
    pub centre: [f64; 3],
    /// Radius (nm).
    pub radius: f64,
}

/// A cylinder defined by its two end-cap centres and radius.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Cylinder {
    /// Centre of the bottom end-cap (nm).
    pub base_centre: [f64; 3],
    /// Axis direction (normalised internally).
    pub axis: [f64; 3],
    /// Length along the axis (nm).
    pub length: f64,
    /// Radius (nm).
    pub radius: f64,
}

/// An axis-aligned cuboid.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Cuboid {
    /// Centre position (nm).
    pub centre: [f64; 3],
    /// Half-extents along x, y, z (nm).
    pub half_extents: [f64; 3],
}

/// A helix defined by its axis, radius, pitch, and number of turns.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Helix {
    /// Centre of the helix base (nm).
    pub base_centre: [f64; 3],
    /// Axis direction.
    pub axis: [f64; 3],
    /// Helix radius (nm).
    pub radius: f64,
    /// Pitch: distance between turns along the axis (nm).
    pub pitch: f64,
    /// Number of turns.
    pub turns: f64,
    /// Wire radius (nm) — cross-section of the helical wire.
    pub wire_radius: f64,
}

/// An ellipsoid defined by its centre and semi-axis lengths.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Ellipsoid {
    /// Centre position (nm).
    pub centre: [f64; 3],
    /// Semi-axis lengths along x, y, z (nm).
    pub semi_axes: [f64; 3],
}

impl Primitive {
    /// Check whether a point lies inside this primitive.
    pub fn contains(&self, point: &[f64; 3]) -> bool {
        match self {
            Primitive::Sphere(s) => {
                let dx = point[0] - s.centre[0];
                let dy = point[1] - s.centre[1];
                let dz = point[2] - s.centre[2];
                dx * dx + dy * dy + dz * dz <= s.radius * s.radius
            }
            Primitive::Cuboid(c) => {
                (point[0] - c.centre[0]).abs() <= c.half_extents[0]
                    && (point[1] - c.centre[1]).abs() <= c.half_extents[1]
                    && (point[2] - c.centre[2]).abs() <= c.half_extents[2]
            }
            Primitive::Ellipsoid(e) => {
                let dx = (point[0] - e.centre[0]) / e.semi_axes[0];
                let dy = (point[1] - e.centre[1]) / e.semi_axes[1];
                let dz = (point[2] - e.centre[2]) / e.semi_axes[2];
                dx * dx + dy * dy + dz * dz <= 1.0
            }
            // TODO: Implement containment for Cylinder and Helix.
            _ => todo!("Containment test for this primitive"),
        }
    }

    /// Axis-aligned bounding box: returns (min_corner, max_corner).
    pub fn bounding_box(&self) -> ([f64; 3], [f64; 3]) {
        match self {
            Primitive::Sphere(s) => (
                [
                    s.centre[0] - s.radius,
                    s.centre[1] - s.radius,
                    s.centre[2] - s.radius,
                ],
                [
                    s.centre[0] + s.radius,
                    s.centre[1] + s.radius,
                    s.centre[2] + s.radius,
                ],
            ),
            Primitive::Cuboid(c) => (
                [
                    c.centre[0] - c.half_extents[0],
                    c.centre[1] - c.half_extents[1],
                    c.centre[2] - c.half_extents[2],
                ],
                [
                    c.centre[0] + c.half_extents[0],
                    c.centre[1] + c.half_extents[1],
                    c.centre[2] + c.half_extents[2],
                ],
            ),
            Primitive::Ellipsoid(e) => (
                [
                    e.centre[0] - e.semi_axes[0],
                    e.centre[1] - e.semi_axes[1],
                    e.centre[2] - e.semi_axes[2],
                ],
                [
                    e.centre[0] + e.semi_axes[0],
                    e.centre[1] + e.semi_axes[1],
                    e.centre[2] + e.semi_axes[2],
                ],
            ),
            _ => todo!("Bounding box for this primitive"),
        }
    }
}
