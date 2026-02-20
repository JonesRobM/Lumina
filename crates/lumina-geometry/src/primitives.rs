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
            Primitive::Cylinder(cyl) => {
                // Normalise axis
                let ax_len = (cyl.axis[0] * cyl.axis[0]
                    + cyl.axis[1] * cyl.axis[1]
                    + cyl.axis[2] * cyl.axis[2])
                    .sqrt();
                if ax_len < 1e-15 {
                    return false;
                }
                let ax = [cyl.axis[0] / ax_len, cyl.axis[1] / ax_len, cyl.axis[2] / ax_len];

                // Vector from base centre to point
                let dp = [
                    point[0] - cyl.base_centre[0],
                    point[1] - cyl.base_centre[1],
                    point[2] - cyl.base_centre[2],
                ];

                // Projection onto axis
                let proj = dp[0] * ax[0] + dp[1] * ax[1] + dp[2] * ax[2];
                if proj < 0.0 || proj > cyl.length {
                    return false;
                }

                // Perpendicular distance squared
                let perp_sq = dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2] - proj * proj;
                perp_sq <= cyl.radius * cyl.radius
            }
            Primitive::Helix(h) => {
                // Normalise axis
                let ax_len = (h.axis[0] * h.axis[0]
                    + h.axis[1] * h.axis[1]
                    + h.axis[2] * h.axis[2])
                    .sqrt();
                if ax_len < 1e-15 {
                    return false;
                }
                let ax = [h.axis[0] / ax_len, h.axis[1] / ax_len, h.axis[2] / ax_len];

                // Build local coordinate system: ax, u, v
                let seed = if ax[0].abs() < 0.9 {
                    [1.0, 0.0, 0.0]
                } else {
                    [0.0, 1.0, 0.0]
                };
                let dot_sa = seed[0] * ax[0] + seed[1] * ax[1] + seed[2] * ax[2];
                let u_raw = [
                    seed[0] - dot_sa * ax[0],
                    seed[1] - dot_sa * ax[1],
                    seed[2] - dot_sa * ax[2],
                ];
                let u_len = (u_raw[0] * u_raw[0] + u_raw[1] * u_raw[1] + u_raw[2] * u_raw[2]).sqrt();
                let u = [u_raw[0] / u_len, u_raw[1] / u_len, u_raw[2] / u_len];
                let v = [
                    ax[1] * u[2] - ax[2] * u[1],
                    ax[2] * u[0] - ax[0] * u[2],
                    ax[0] * u[1] - ax[1] * u[0],
                ];

                // Vector from base centre to point
                let dp = [
                    point[0] - h.base_centre[0],
                    point[1] - h.base_centre[1],
                    point[2] - h.base_centre[2],
                ];

                // Project onto local coordinates
                let z = dp[0] * ax[0] + dp[1] * ax[1] + dp[2] * ax[2];
                let total_height = h.pitch * h.turns;
                if z < -h.wire_radius || z > total_height + h.wire_radius {
                    return false;
                }

                let pu = dp[0] * u[0] + dp[1] * u[1] + dp[2] * u[2];
                let pv = dp[0] * v[0] + dp[1] * v[1] + dp[2] * v[2];

                // Find nearest point on helix centreline by scanning parameter t
                // The centreline is: C(t) = base + t*pitch/(2π)*ax + R*(cos(t)*u + sin(t)*v)
                // where t ranges from 0 to 2π*turns.
                let t_max = 2.0 * std::f64::consts::PI * h.turns;
                // Initial estimate from axial position
                let t_guess = if h.pitch > 1e-15 {
                    z * 2.0 * std::f64::consts::PI / h.pitch
                } else {
                    0.0
                };
                let t_guess = t_guess.clamp(0.0, t_max);

                // Search in a window around the guess
                let num_samples = 32;
                let t_window = 2.0 * std::f64::consts::PI; // one full turn
                let t_lo = (t_guess - t_window).max(0.0);
                let t_hi = (t_guess + t_window).min(t_max);

                let mut min_dist_sq = f64::INFINITY;
                for si in 0..=num_samples {
                    let t = t_lo + (t_hi - t_lo) * si as f64 / num_samples as f64;
                    let cz = t * h.pitch / (2.0 * std::f64::consts::PI);
                    let cu = h.radius * t.cos();
                    let cv = h.radius * t.sin();
                    let dz = z - cz;
                    let du = pu - cu;
                    let dv = pv - cv;
                    let dist_sq = dz * dz + du * du + dv * dv;
                    if dist_sq < min_dist_sq {
                        min_dist_sq = dist_sq;
                    }
                }

                min_dist_sq <= h.wire_radius * h.wire_radius
            }
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
            Primitive::Cylinder(cyl) => {
                // Normalise axis
                let ax_len = (cyl.axis[0] * cyl.axis[0]
                    + cyl.axis[1] * cyl.axis[1]
                    + cyl.axis[2] * cyl.axis[2])
                    .sqrt();
                let ax = if ax_len > 1e-15 {
                    [cyl.axis[0] / ax_len, cyl.axis[1] / ax_len, cyl.axis[2] / ax_len]
                } else {
                    [0.0, 0.0, 1.0]
                };

                // Top centre
                let top = [
                    cyl.base_centre[0] + ax[0] * cyl.length,
                    cyl.base_centre[1] + ax[1] * cyl.length,
                    cyl.base_centre[2] + ax[2] * cyl.length,
                ];

                // For each axis i, the radius contribution perpendicular to the
                // cylinder axis is r * sqrt(1 - ax[i]^2).
                let mut bb_min = [0.0f64; 3];
                let mut bb_max = [0.0f64; 3];
                for i in 0..3 {
                    let r_extent = cyl.radius * (1.0 - ax[i] * ax[i]).max(0.0).sqrt();
                    let lo = cyl.base_centre[i].min(top[i]) - r_extent;
                    let hi = cyl.base_centre[i].max(top[i]) + r_extent;
                    bb_min[i] = lo;
                    bb_max[i] = hi;
                }
                (bb_min, bb_max)
            }
            Primitive::Helix(h) => {
                let ax_len = (h.axis[0] * h.axis[0]
                    + h.axis[1] * h.axis[1]
                    + h.axis[2] * h.axis[2])
                    .sqrt();
                let ax = if ax_len > 1e-15 {
                    [h.axis[0] / ax_len, h.axis[1] / ax_len, h.axis[2] / ax_len]
                } else {
                    [0.0, 0.0, 1.0]
                };

                // The helix sweeps a tube of radius wire_radius around a helical
                // centreline at distance h.radius from the axis. The bounding box
                // extends (h.radius + h.wire_radius) perpendicular to the axis,
                // and from 0 to pitch*turns along the axis (plus wire_radius).
                let total_height = h.pitch * h.turns;
                let outer_r = h.radius + h.wire_radius;

                let mut bb_min = [0.0f64; 3];
                let mut bb_max = [0.0f64; 3];
                for i in 0..3 {
                    let r_extent = outer_r * (1.0 - ax[i] * ax[i]).max(0.0).sqrt();
                    let base = h.base_centre[i];
                    let top = h.base_centre[i] + total_height * ax[i];
                    bb_min[i] = base.min(top) - r_extent - h.wire_radius;
                    bb_max[i] = base.max(top) + r_extent + h.wire_radius;
                }
                (bb_min, bb_max)
            }
        }
    }
}
