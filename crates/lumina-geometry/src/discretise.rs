//! Discretisation of geometric shapes into dipole lattices.
//!
//! Given a [`Primitive`](crate::primitives::Primitive), this module generates
//! a cubic lattice of points within the shape's volume, spaced by the
//! user-specified dipole spacing $d$. The spacing must satisfy the CDA
//! validity criterion $d \ll \lambda / (2\pi |m|)$ where $m$ is the complex
//! refractive index.

use crate::primitives::Primitive;

/// A point in a dipole lattice, produced by discretisation.
#[derive(Debug, Clone)]
pub struct LatticePoint {
    /// Position in 3D space (nm).
    pub position: [f64; 3],
}

/// Discretise a primitive shape into a cubic lattice of dipole positions.
///
/// # Arguments
/// * `primitive` - The shape to discretise.
/// * `spacing` - Lattice spacing in nanometres.
///
/// # Returns
/// A vector of lattice points lying within the shape.
pub fn discretise_primitive(primitive: &Primitive, spacing: f64) -> Vec<LatticePoint> {
    assert!(spacing > 0.0, "Dipole spacing must be positive");

    let (min, max) = primitive.bounding_box();
    let mut points = Vec::new();

    let mut x = min[0];
    while x <= max[0] {
        let mut y = min[1];
        while y <= max[1] {
            let mut z = min[2];
            while z <= max[2] {
                let p = [x, y, z];
                if primitive.contains(&p) {
                    points.push(LatticePoint { position: p });
                }
                z += spacing;
            }
            y += spacing;
        }
        x += spacing;
    }

    points
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::Sphere;

    #[test]
    fn test_sphere_discretisation_count() {
        let sphere = Primitive::Sphere(Sphere {
            centre: [0.0, 0.0, 0.0],
            radius: 10.0,
        });
        let points = discretise_primitive(&sphere, 2.0);

        // A sphere of radius 10 nm with 2 nm spacing should produce roughly
        // (4/3 * pi * 10^3) / 2^3 ~ 524 dipoles. Allow generous tolerance.
        assert!(points.len() > 300, "Too few dipoles: {}", points.len());
        assert!(points.len() < 800, "Too many dipoles: {}", points.len());
    }

    #[test]
    fn test_all_points_inside_sphere() {
        let sphere = Primitive::Sphere(Sphere {
            centre: [5.0, 5.0, 5.0],
            radius: 8.0,
        });
        let points = discretise_primitive(&sphere, 1.5);

        for p in &points {
            let dx = p.position[0] - 5.0;
            let dy = p.position[1] - 5.0;
            let dz = p.position[2] - 5.0;
            let r_sq = dx * dx + dy * dy + dz * dz;
            assert!(r_sq <= 8.0 * 8.0 + 1e-10);
        }
    }
}
