//! Discretisation of geometric shapes into dipole lattices.
//!
//! Given a [`Primitive`](crate::primitives::Primitive), this module generates
//! a cubic lattice of points within the shape's volume, spaced by the
//! user-specified dipole spacing $d$. The spacing must satisfy the CDA
//! validity criterion $d \ll \lambda / (2\pi |m|)$ where $m$ is the complex
//! refractive index.

use crate::parsers::obj::ObjMesh;
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

    // Centre the lattice on the midpoint of the bounding box so that the
    // grid is symmetric about the shape centre.
    let cx = 0.5 * (min[0] + max[0]);
    let cy = 0.5 * (min[1] + max[1]);
    let cz = 0.5 * (min[2] + max[2]);

    let nx = ((max[0] - cx) / spacing).floor() as i32;
    let ny = ((max[1] - cy) / spacing).floor() as i32;
    let nz = ((max[2] - cz) / spacing).floor() as i32;

    for ix in -nx..=nx {
        let x = cx + ix as f64 * spacing;
        for iy in -ny..=ny {
            let y = cy + iy as f64 * spacing;
            for iz in -nz..=nz {
                let z = cz + iz as f64 * spacing;
                let p = [x, y, z];
                if primitive.contains(&p) {
                    points.push(LatticePoint { position: p });
                }
            }
        }
    }

    points
}

/// Discretise a closed triangle mesh into a cubic lattice of dipole positions.
///
/// Uses ray-casting (Möller–Trumbore along +z) to test whether each grid point
/// lies inside the mesh. The lattice is centred on the bounding-box midpoint,
/// following the same convention as [`discretise_primitive`].
///
/// # Arguments
/// * `mesh` - A closed triangle mesh parsed from an OBJ file.
/// * `spacing` - Lattice spacing in nanometres.
///
/// # Returns
/// A vector of lattice points lying within the mesh volume.
pub fn discretise_mesh(mesh: &ObjMesh, spacing: f64) -> Vec<LatticePoint> {
    assert!(spacing > 0.0, "Dipole spacing must be positive");

    let (min, max) = mesh.bounding_box();
    let mut points = Vec::new();

    // Centre the lattice on the midpoint of the bounding box.
    let cx = 0.5 * (min[0] + max[0]);
    let cy = 0.5 * (min[1] + max[1]);
    let cz = 0.5 * (min[2] + max[2]);

    let nx = ((max[0] - cx) / spacing).floor() as i32;
    let ny = ((max[1] - cy) / spacing).floor() as i32;
    let nz = ((max[2] - cz) / spacing).floor() as i32;

    for ix in -nx..=nx {
        let x = cx + ix as f64 * spacing;
        for iy in -ny..=ny {
            let y = cy + iy as f64 * spacing;
            for iz in -nz..=nz {
                let z = cz + iz as f64 * spacing;
                let p = [x, y, z];
                if mesh_contains(&mesh.vertices, &mesh.faces, &p) {
                    points.push(LatticePoint { position: p });
                }
            }
        }
    }

    points
}

/// Test whether a point is inside a closed triangle mesh using ray casting.
///
/// Casts a ray from `point` along +z and counts intersections with the mesh
/// triangles. An odd intersection count means the point is inside.
///
/// A tiny perturbation is applied to the xy coordinates to avoid the degenerate
/// case where the ray passes exactly through a shared edge or vertex of the
/// mesh, which would cause double-counting.
fn mesh_contains(vertices: &[[f64; 3]], faces: &[[usize; 3]], point: &[f64; 3]) -> bool {
    // Perturb xy to avoid hitting shared edges/vertices exactly.
    let perturbed = [point[0] + 1.23e-10, point[1] + 2.34e-10, point[2]];
    let mut count = 0u32;
    for face in faces {
        if ray_triangle_z(
            &perturbed,
            &vertices[face[0]],
            &vertices[face[1]],
            &vertices[face[2]],
        )
        .is_some()
        {
            count += 1;
        }
    }
    count % 2 == 1
}

/// Möller–Trumbore ray-triangle intersection for a ray along +z.
///
/// Returns `Some(t)` if the ray from `origin` in direction [0, 0, 1] intersects
/// the triangle `(v0, v1, v2)` at parameter `t > 0`. Returns `None` otherwise.
fn ray_triangle_z(origin: &[f64; 3], v0: &[f64; 3], v1: &[f64; 3], v2: &[f64; 3]) -> Option<f64> {
    // Edge vectors
    let e1 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
    let e2 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];

    // P = D × E2 where D = [0, 0, 1]
    // P = [0*e2z - 1*e2y, 1*e2x - 0*e2z, 0*e2y - 0*e2x] = [-e2y, e2x, 0]
    let px = -e2[1];
    let py = e2[0];

    // det = E1 · P = e1x * (-e2y) + e1y * e2x
    let det = e1[0] * px + e1[1] * py;
    if det.abs() < 1e-12 {
        return None; // Ray parallel to triangle
    }
    let inv_det = 1.0 / det;

    // T = origin - V0
    let tx = origin[0] - v0[0];
    let ty = origin[1] - v0[1];
    let tz = origin[2] - v0[2];

    // u = (T · P) * inv_det
    let u = (tx * px + ty * py) * inv_det;
    if !(0.0..=1.0).contains(&u) {
        return None;
    }

    // Q = T × E1
    let qx = ty * e1[2] - tz * e1[1];
    let qy = tz * e1[0] - tx * e1[2];
    let qz = tx * e1[1] - ty * e1[0];

    // v = (D · Q) * inv_det, D = [0, 0, 1] so D · Q = qz
    let v = qz * inv_det;
    if v < 0.0 || u + v > 1.0 {
        return None;
    }

    // t = (E2 · Q) * inv_det
    let t = (e2[0] * qx + e2[1] * qy + e2[2] * qz) * inv_det;
    if t > 0.0 {
        Some(t)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parsers::obj::parse_obj;
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

    /// Build a cube OBJ string for testing. Cube spans ±h on each axis (6 quad faces).
    fn cube_obj(h: f64) -> String {
        format!(
            "v {h} {h} -{h}\nv {h} -{h} -{h}\nv -{h} -{h} -{h}\nv -{h} {h} -{h}\n\
             v {h} {h} {h}\nv {h} -{h} {h}\nv -{h} -{h} {h}\nv -{h} {h} {h}\n\
             f 1 2 3 4\nf 5 8 7 6\nf 1 5 6 2\nf 3 7 8 4\nf 1 4 8 5\nf 2 6 7 3\n"
        )
    }

    #[test]
    fn test_ray_triangle_z_hit() {
        // Triangle in z=5 plane: (0,0,5), (10,0,5), (0,10,5)
        let v0 = [0.0, 0.0, 5.0];
        let v1 = [10.0, 0.0, 5.0];
        let v2 = [0.0, 10.0, 5.0];
        // Ray from (1,1,0) along +z should hit at t=5
        let origin = [1.0, 1.0, 0.0];
        let t = ray_triangle_z(&origin, &v0, &v1, &v2);
        assert!(t.is_some());
        assert!((t.unwrap() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_ray_triangle_z_miss() {
        let v0 = [0.0, 0.0, 5.0];
        let v1 = [10.0, 0.0, 5.0];
        let v2 = [0.0, 10.0, 5.0];
        // Ray from (20,20,0) misses the triangle entirely
        let origin = [20.0, 20.0, 0.0];
        assert!(ray_triangle_z(&origin, &v0, &v1, &v2).is_none());
    }

    #[test]
    fn test_ray_triangle_z_behind() {
        let v0 = [0.0, 0.0, -5.0];
        let v1 = [10.0, 0.0, -5.0];
        let v2 = [0.0, 10.0, -5.0];
        // Triangle at z=-5 is behind origin at z=0 → t < 0, should return None
        let origin = [1.0, 1.0, 0.0];
        assert!(ray_triangle_z(&origin, &v0, &v1, &v2).is_none());
    }

    #[test]
    fn test_mesh_contains_cube() {
        let mesh = parse_obj(&cube_obj(5.0)).unwrap();
        // Origin should be inside
        assert!(mesh_contains(&mesh.vertices, &mesh.faces, &[0.0, 0.0, 0.0]));
        // Points well inside
        assert!(mesh_contains(&mesh.vertices, &mesh.faces, &[3.0, 3.0, 3.0]));
        assert!(mesh_contains(&mesh.vertices, &mesh.faces, &[-4.0, -4.0, -4.0]));
        // Points outside
        assert!(!mesh_contains(&mesh.vertices, &mesh.faces, &[6.0, 0.0, 0.0]));
        assert!(!mesh_contains(&mesh.vertices, &mesh.faces, &[0.0, 6.0, 0.0]));
        assert!(!mesh_contains(&mesh.vertices, &mesh.faces, &[0.0, 0.0, 6.0]));
        assert!(!mesh_contains(&mesh.vertices, &mesh.faces, &[10.0, 10.0, 10.0]));
    }

    #[test]
    fn test_discretise_cube_mesh() {
        // Cube ±5nm at 2nm spacing: same as parametric cuboid test
        // Centred grid: -4,-2,0,2,4 per axis (5 points) → 5³ = 125 interior points
        // But boundary points at ±5 land exactly ON the surface; ray-casting may
        // include or exclude them. With the full ±5 grid: points from -4 to +4
        // in steps of 2 → 5 per axis → 125 minimum.
        let mesh = parse_obj(&cube_obj(5.0)).unwrap();
        let points = discretise_mesh(&mesh, 2.0);
        // Should get at least 125 points (interior) and at most ~729 (if boundary included)
        assert!(
            points.len() >= 100 && points.len() <= 800,
            "Cube mesh dipole count = {}, expected ~125-343",
            points.len()
        );
        // All points should be within the cube ±5
        for p in &points {
            assert!(p.position[0].abs() <= 5.0 + 1e-10);
            assert!(p.position[1].abs() <= 5.0 + 1e-10);
            assert!(p.position[2].abs() <= 5.0 + 1e-10);
        }
    }
}
