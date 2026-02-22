//! Parser for Wavefront `.obj` mesh files.
//!
//! OBJ files define 3D meshes as collections of vertices (`v`) and faces (`f`).
//! The mesh surface is used to define the boundary of a nanostructure, which
//! is then filled with dipoles by the discretisation module.
//!
//! Coordinates are treated as nanometres (no unit conversion).
//! Faces are triangulated on parse (quads and n-gons use fan triangulation).

use super::ParseError;

/// A parsed triangle mesh from an OBJ file.
#[derive(Debug, Clone)]
pub struct ObjMesh {
    /// Vertex positions in nm (0-indexed).
    pub vertices: Vec<[f64; 3]>,
    /// Triangulated face indices (0-indexed into `vertices`).
    pub faces: Vec<[usize; 3]>,
}

/// Parse an OBJ file, extracting vertices and triangulated faces.
///
/// Handles the following OBJ elements:
/// - `v x y z` — vertex position
/// - `f v1 v2 v3 ...` — face (triangulated via fan from first vertex)
/// - `f v1/vt1 ...` or `f v1/vt1/vn1 ...` or `f v1//vn1 ...` — takes vertex index only
///
/// All other lines (`vt`, `vn`, `g`, `o`, `mtllib`, `usemtl`, `s`, `#`, blank) are ignored.
pub fn parse_obj(content: &str) -> Result<ObjMesh, ParseError> {
    let mut vertices = Vec::new();
    let mut faces = Vec::new();

    for (line_idx, line) in content.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let mut parts = line.split_whitespace();
        let keyword = match parts.next() {
            Some(k) => k,
            None => continue,
        };

        match keyword {
            "v" => {
                let coords: Vec<&str> = parts.collect();
                if coords.len() < 3 {
                    return Err(ParseError::FormatError {
                        line: line_idx + 1,
                        message: format!("Vertex needs at least 3 coordinates, got {}", coords.len()),
                    });
                }
                let x: f64 = coords[0].parse().map_err(|_| ParseError::FormatError {
                    line: line_idx + 1,
                    message: format!("Invalid x coordinate: {}", coords[0]),
                })?;
                let y: f64 = coords[1].parse().map_err(|_| ParseError::FormatError {
                    line: line_idx + 1,
                    message: format!("Invalid y coordinate: {}", coords[1]),
                })?;
                let z: f64 = coords[2].parse().map_err(|_| ParseError::FormatError {
                    line: line_idx + 1,
                    message: format!("Invalid z coordinate: {}", coords[2]),
                })?;
                vertices.push([x, y, z]);
            }
            "f" => {
                let indices: Result<Vec<usize>, ParseError> = parts
                    .enumerate()
                    .map(|(i, token)| {
                        // Handle v, v/vt, v/vt/vn, v//vn formats — take first element
                        let idx_str = token.split('/').next().unwrap_or(token);
                        let idx: usize = idx_str.parse().map_err(|_| ParseError::FormatError {
                            line: line_idx + 1,
                            message: format!("Invalid face index at position {}: {}", i + 1, token),
                        })?;
                        if idx == 0 {
                            return Err(ParseError::FormatError {
                                line: line_idx + 1,
                                message: "Face index 0 is invalid (OBJ indices are 1-based)".into(),
                            });
                        }
                        Ok(idx - 1) // Convert 1-based to 0-based
                    })
                    .collect();
                let indices = indices?;

                if indices.len() < 3 {
                    return Err(ParseError::FormatError {
                        line: line_idx + 1,
                        message: format!("Face needs at least 3 vertices, got {}", indices.len()),
                    });
                }

                // Fan-triangulate from first vertex: (v0,v1,v2), (v0,v2,v3), ...
                for i in 1..indices.len() - 1 {
                    faces.push([indices[0], indices[i], indices[i + 1]]);
                }
            }
            // Ignore everything else (vt, vn, g, o, s, mtllib, usemtl, etc.)
            _ => {}
        }
    }

    // Validate
    if vertices.len() < 4 {
        return Err(ParseError::FormatError {
            line: 0,
            message: format!(
                "Mesh needs at least 4 vertices to form a closed volume, got {}",
                vertices.len()
            ),
        });
    }
    if faces.is_empty() {
        return Err(ParseError::FormatError {
            line: 0,
            message: "No faces found in OBJ file".into(),
        });
    }

    // Check all face indices are in bounds
    let n = vertices.len();
    for (fi, face) in faces.iter().enumerate() {
        for &idx in face {
            if idx >= n {
                return Err(ParseError::FormatError {
                    line: 0,
                    message: format!(
                        "Face {} references vertex index {} but only {} vertices exist",
                        fi + 1,
                        idx + 1,
                        n
                    ),
                });
            }
        }
    }

    Ok(ObjMesh { vertices, faces })
}

impl ObjMesh {
    /// Compute the axis-aligned bounding box as `(min_corner, max_corner)`.
    pub fn bounding_box(&self) -> ([f64; 3], [f64; 3]) {
        let mut min = [f64::INFINITY; 3];
        let mut max = [f64::NEG_INFINITY; 3];
        for v in &self.vertices {
            for i in 0..3 {
                if v[i] < min[i] {
                    min[i] = v[i];
                }
                if v[i] > max[i] {
                    max[i] = v[i];
                }
            }
        }
        (min, max)
    }
}

/// Helper: build a cube OBJ string for testing. Cube spans ±half_size on each axis.
#[cfg(test)]
fn cube_obj(half_size: f64) -> String {
    let h = half_size;
    format!(
        "# Unit cube\n\
         v {h} {h} -{h}\n\
         v {h} -{h} -{h}\n\
         v -{h} -{h} -{h}\n\
         v -{h} {h} -{h}\n\
         v {h} {h} {h}\n\
         v {h} -{h} {h}\n\
         v -{h} -{h} {h}\n\
         v -{h} {h} {h}\n\
         f 1 2 3 4\n\
         f 5 8 7 6\n\
         f 1 5 6 2\n\
         f 3 7 8 4\n\
         f 1 4 8 5\n\
         f 2 6 7 3\n"
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_cube_obj() {
        let obj = cube_obj(5.0);
        let mesh = parse_obj(&obj).unwrap();
        assert_eq!(mesh.vertices.len(), 8);
        // 6 quad faces → 12 triangles
        assert_eq!(mesh.faces.len(), 12);
        // All indices in [0, 7]
        for face in &mesh.faces {
            for &idx in face {
                assert!(idx < 8);
            }
        }
    }

    #[test]
    fn test_parse_triangle_faces() {
        let obj = "\
            v 0 0 0\nv 1 0 0\nv 0 1 0\nv 0 0 1\n\
            f 1 2 3\nf 1 2 4\nf 1 3 4\nf 2 3 4\n";
        let mesh = parse_obj(obj).unwrap();
        assert_eq!(mesh.vertices.len(), 4);
        assert_eq!(mesh.faces.len(), 4); // Already triangles, no fan-split
    }

    #[test]
    fn test_parse_quad_faces_triangulated() {
        let obj = "\
            v 0 0 0\nv 10 0 0\nv 10 10 0\nv 0 10 0\n\
            v 0 0 10\nv 10 0 10\nv 10 10 10\nv 0 10 10\n\
            f 1 2 3 4\nf 5 6 7 8\n";
        let mesh = parse_obj(obj).unwrap();
        // 2 quads → 4 triangles
        assert_eq!(mesh.faces.len(), 4);
    }

    #[test]
    fn test_parse_with_normals_and_texcoords() {
        let obj = "\
            v 0 0 0\nv 1 0 0\nv 0 1 0\nv 0 0 1\n\
            vt 0 0\nvt 1 0\nvt 0 1\n\
            vn 0 0 1\nvn 0 1 0\nvn 1 0 0\n\
            f 1/1/1 2/2/2 3/3/3\n\
            f 1//1 2//2 4//3\n";
        let mesh = parse_obj(obj).unwrap();
        assert_eq!(mesh.faces.len(), 2);
        assert_eq!(mesh.faces[0], [0, 1, 2]);
        assert_eq!(mesh.faces[1], [0, 1, 3]);
    }

    #[test]
    fn test_parse_empty_returns_error() {
        let result = parse_obj("");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_no_faces_returns_error() {
        let obj = "v 0 0 0\nv 1 0 0\nv 0 1 0\nv 0 0 1\n";
        let result = parse_obj(obj);
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_out_of_bounds_face_index() {
        let obj = "v 0 0 0\nv 1 0 0\nv 0 1 0\nv 0 0 1\nf 1 2 5\n";
        let result = parse_obj(obj);
        assert!(result.is_err());
    }

    #[test]
    fn test_bounding_box() {
        let obj = cube_obj(5.0);
        let mesh = parse_obj(&obj).unwrap();
        let (min, max) = mesh.bounding_box();
        for i in 0..3 {
            assert!((min[i] - (-5.0)).abs() < 1e-10);
            assert!((max[i] - 5.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_comments_and_whitespace_ignored() {
        let obj = "\
            # This is a comment\n\
            \n\
            v 0 0 0\n\
            v 1 0 0\n\
            # Another comment\n\
            v 0 1 0\n\
            v 0 0 1\n\
            \n\
            f 1 2 3\n\
            f 1 2 4\n\
            f 1 3 4\n\
            f 2 3 4\n";
        let mesh = parse_obj(obj).unwrap();
        assert_eq!(mesh.vertices.len(), 4);
        assert_eq!(mesh.faces.len(), 4);
    }
}
