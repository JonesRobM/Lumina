//! Parser for `.xyz` molecular coordinate files.
//!
//! The XYZ format is a simple plain-text format:
//! ```text
//! <num_atoms>
//! <comment line>
//! <element> <x> <y> <z>
//! <element> <x> <y> <z>
//! ...
//! ```
//!
//! Coordinates are in angstroms. The parser converts to nanometres
//! (1 angstrom = 0.1 nm) for internal use.

use super::{ParseError, ParsedPoint};

const ANGSTROM_TO_NM: f64 = 0.1;

/// Parse an XYZ file from a string.
pub fn parse_xyz(content: &str) -> Result<Vec<ParsedPoint>, ParseError> {
    let lines: Vec<&str> = content.lines().collect();

    if lines.len() < 3 {
        return Err(ParseError::FormatError {
            line: 1,
            message: "XYZ file must have at least 3 lines".into(),
        });
    }

    let num_atoms: usize = lines[0].trim().parse().map_err(|_| ParseError::FormatError {
        line: 1,
        message: "First line must be the number of atoms".into(),
    })?;

    // Line 2 is the comment line (ignored)

    let mut points = Vec::with_capacity(num_atoms);
    for (idx, line) in lines[2..].iter().enumerate() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 4 {
            return Err(ParseError::FormatError {
                line: idx + 3,
                message: format!("Expected 'element x y z', got '{}'", line),
            });
        }

        let x: f64 = parts[1].parse().map_err(|_| ParseError::FormatError {
            line: idx + 3,
            message: format!("Invalid x coordinate: {}", parts[1]),
        })?;
        let y: f64 = parts[2].parse().map_err(|_| ParseError::FormatError {
            line: idx + 3,
            message: format!("Invalid y coordinate: {}", parts[2]),
        })?;
        let z: f64 = parts[3].parse().map_err(|_| ParseError::FormatError {
            line: idx + 3,
            message: format!("Invalid z coordinate: {}", parts[3]),
        })?;

        points.push(ParsedPoint {
            position: [x * ANGSTROM_TO_NM, y * ANGSTROM_TO_NM, z * ANGSTROM_TO_NM],
            label: Some(parts[0].to_string()),
        });
    }

    if points.len() != num_atoms {
        return Err(ParseError::FormatError {
            line: 1,
            message: format!("Header says {} atoms but found {}", num_atoms, points.len()),
        });
    }

    Ok(points)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_simple_xyz() {
        let content = "3\nTest molecule\nAu 0.0 0.0 0.0\nAu 2.88 0.0 0.0\nAu 1.44 2.49 0.0\n";
        let points = parse_xyz(content).unwrap();
        assert_eq!(points.len(), 3);
        assert_eq!(points[0].label.as_deref(), Some("Au"));
        // 2.88 angstrom = 0.288 nm
        assert!((points[1].position[0] - 0.288).abs() < 1e-10);
    }

    // ── Mackay icosahedron generator ────────────────────────────────────

    /// Generate Mackay icosahedron coordinates for shell number `k`.
    ///
    /// Uses the 12 vertices of a regular icosahedron as basis directions.
    /// Atoms are placed at all unique positions on vertices, edges, and
    /// triangular faces of each shell, following the standard close-packed
    /// Mackay construction.
    ///
    /// Returns positions in angstroms (Au nearest-neighbour = 2.884 Å).
    fn mackay_icosahedron(k: usize) -> Vec<[f64; 3]> {
        // Au lattice constant 4.078 Å → nn distance = a/√2
        let nn: f64 = 4.078 / std::f64::consts::SQRT_2; // ≈ 2.884 Å

        // 12 vertices of a regular icosahedron (unit vectors)
        let phi: f64 = (1.0 + 5.0_f64.sqrt()) / 2.0; // golden ratio
        let norm = (1.0 + phi * phi).sqrt();
        let a = 1.0 / norm;
        let b = phi / norm;

        let vertices: [[f64; 3]; 12] = [
            [ 0.0,  a,  b], [ 0.0,  a, -b], [ 0.0, -a,  b], [ 0.0, -a, -b],
            [ a,  b, 0.0], [ a, -b, 0.0], [-a,  b, 0.0], [-a, -b, 0.0],
            [ b, 0.0,  a], [-b, 0.0,  a], [ b, 0.0, -a], [-b, 0.0, -a],
        ];

        // 30 edges of the icosahedron (pairs of vertex indices)
        let edges: Vec<[usize; 2]> = find_icosahedron_edges(&vertices);
        assert_eq!(edges.len(), 30, "Icosahedron must have 30 edges");

        // 20 triangular faces
        let faces: Vec<[usize; 3]> = find_icosahedron_faces(&vertices, &edges);
        assert_eq!(faces.len(), 20, "Icosahedron must have 20 faces");

        let mut positions: Vec<[f64; 3]> = Vec::new();

        // Shell 0: just the origin
        positions.push([0.0, 0.0, 0.0]);

        for shell in 1..=k {
            let s = shell as f64;
            let r = s * nn;

            // 12 vertex atoms: at r * vertex_unit_vector
            for v in &vertices {
                positions.push([r * v[0], r * v[1], r * v[2]]);
            }

            // Edge atoms: (shell - 1) atoms per edge, interpolated
            if shell > 1 {
                for edge in &edges {
                    let v0 = vertices[edge[0]];
                    let v1 = vertices[edge[1]];
                    for j in 1..shell {
                        let t = j as f64 / s;
                        let pos = [
                            r * (v0[0] * (1.0 - t) + v1[0] * t),
                            r * (v0[1] * (1.0 - t) + v1[1] * t),
                            r * (v0[2] * (1.0 - t) + v1[2] * t),
                        ];
                        positions.push(pos);
                    }
                }
            }

            // Face atoms: interior of each triangular face
            if shell > 2 {
                for face in &faces {
                    let v0 = vertices[face[0]];
                    let v1 = vertices[face[1]];
                    let v2 = vertices[face[2]];
                    // Barycentric coordinates (i, j) with i >= 1, j >= 1, i+j <= shell-1
                    for i in 1..shell - 1 {
                        for j in 1..shell - i {
                            let fi = i as f64 / s;
                            let fj = j as f64 / s;
                            let fk = 1.0 - fi - fj;
                            let pos = [
                                r * (fk * v0[0] + fi * v1[0] + fj * v2[0]),
                                r * (fk * v0[1] + fi * v1[1] + fj * v2[1]),
                                r * (fk * v0[2] + fi * v1[2] + fj * v2[2]),
                            ];
                            positions.push(pos);
                        }
                    }
                }
            }
        }

        positions
    }

    /// Find the 30 edges of an icosahedron from its 12 vertices.
    /// Two vertices are connected if their dot product equals 1/√5
    /// (angular separation ≈ 63.43°, the edge angle).
    fn find_icosahedron_edges(vertices: &[[f64; 3]; 12]) -> Vec<[usize; 2]> {
        let mut edges = Vec::new();
        // For unit vectors on an icosahedron, connected vertices have
        // dot product = 1/√5 ≈ 0.4472
        let target_dot = 1.0 / 5.0_f64.sqrt();
        for i in 0..12 {
            for j in (i + 1)..12 {
                let dot = vertices[i][0] * vertices[j][0]
                    + vertices[i][1] * vertices[j][1]
                    + vertices[i][2] * vertices[j][2];
                if (dot - target_dot).abs() < 0.01 {
                    edges.push([i, j]);
                }
            }
        }
        edges
    }

    /// Find the 20 triangular faces from edges.
    /// A face is a triple (i, j, k) where all three pairs are edges.
    fn find_icosahedron_faces(
        vertices: &[[f64; 3]; 12],
        edges: &[[usize; 2]],
    ) -> Vec<[usize; 3]> {
        let _ = vertices; // used only for type context
        let mut adj = [[false; 12]; 12];
        for e in edges {
            adj[e[0]][e[1]] = true;
            adj[e[1]][e[0]] = true;
        }
        let mut faces = Vec::new();
        for i in 0..12 {
            for j in (i + 1)..12 {
                if !adj[i][j] {
                    continue;
                }
                for m in (j + 1)..12 {
                    if adj[i][m] && adj[j][m] {
                        faces.push([i, j, m]);
                    }
                }
            }
        }
        faces
    }

    /// Generate XYZ file string from atom positions.
    fn positions_to_xyz(element: &str, positions: &[[f64; 3]]) -> String {
        let mut s = format!("{}\nMackay icosahedron\n", positions.len());
        for p in positions {
            s.push_str(&format!("{} {:.6} {:.6} {:.6}\n", element, p[0], p[1], p[2]));
        }
        s
    }

    // ── Tests ───────────────────────────────────────────────────────────

    /// Verify Mackay atom counts: N(k) = (10k³ + 15k² + 11k + 3) / 3
    #[test]
    fn test_mackay_atom_counts() {
        // Known counts: k=1→13, k=2→55, k=3→147, k=4→309, k=5→561, k=6→923
        let expected = [(1, 13), (2, 55), (3, 147), (4, 309), (5, 561), (6, 923)];
        for (k, n) in expected {
            let positions = mackay_icosahedron(k);
            assert_eq!(
                positions.len(),
                n,
                "Mackay k={}: expected {} atoms, got {}",
                k,
                n,
                positions.len()
            );
        }
    }

    /// Parse a k=6 Au Mackay icosahedron (923 atoms) from generated XYZ.
    #[test]
    fn test_parse_au_icosahedron_k6() {
        let positions = mackay_icosahedron(6);
        assert_eq!(positions.len(), 923);

        let xyz_content = positions_to_xyz("Au", &positions);
        let points = parse_xyz(&xyz_content).unwrap();

        // Correct atom count
        assert_eq!(points.len(), 923);

        // All labels are Au
        assert!(points.iter().all(|p| p.label.as_deref() == Some("Au")));

        // First atom is the origin (shell 0)
        let origin = &points[0].position;
        assert!(
            origin[0].abs() < 1e-12 && origin[1].abs() < 1e-12 && origin[2].abs() < 1e-12,
            "First atom should be at origin, got {:?}",
            origin
        );

        // Positions are in nm (Å * 0.1)
        let nn_nm = 4.078 / std::f64::consts::SQRT_2 * 0.1; // ≈ 0.2884 nm
        // Shell-1 vertex atoms should be at distance nn_nm from origin
        let r1 = (points[1].position[0].powi(2)
            + points[1].position[1].powi(2)
            + points[1].position[2].powi(2))
        .sqrt();
        assert!(
            (r1 - nn_nm).abs() < 1e-6,
            "Shell-1 atom distance: expected {:.6} nm, got {:.6} nm",
            nn_nm,
            r1
        );

        // Outermost shell (k=6) atoms should be at distance 6 * nn_nm
        let r_max_expected = 6.0 * nn_nm;
        let r_max = points
            .iter()
            .map(|p| {
                (p.position[0].powi(2) + p.position[1].powi(2) + p.position[2].powi(2)).sqrt()
            })
            .fold(0.0_f64, f64::max);
        assert!(
            (r_max - r_max_expected).abs() < 1e-6,
            "Max radius: expected {:.6} nm, got {:.6} nm",
            r_max_expected,
            r_max
        );

        // Icosahedral symmetry check: centre of mass should be at origin
        let com: [f64; 3] = [
            points.iter().map(|p| p.position[0]).sum::<f64>() / 923.0,
            points.iter().map(|p| p.position[1]).sum::<f64>() / 923.0,
            points.iter().map(|p| p.position[2]).sum::<f64>() / 923.0,
        ];
        let com_r = (com[0].powi(2) + com[1].powi(2) + com[2].powi(2)).sqrt();
        assert!(
            com_r < 1e-10,
            "Centre of mass should be at origin, got {:?} (r={:.2e})",
            com,
            com_r
        );

        // No duplicate positions (all atoms should be distinct)
        let mut unique = points.iter().map(|p| p.position).collect::<Vec<_>>();
        unique.sort_by(|a, b| {
            a[0].partial_cmp(&b[0])
                .unwrap()
                .then(a[1].partial_cmp(&b[1]).unwrap())
                .then(a[2].partial_cmp(&b[2]).unwrap())
        });
        unique.dedup_by(|a, b| {
            (a[0] - b[0]).abs() < 1e-10
                && (a[1] - b[1]).abs() < 1e-10
                && (a[2] - b[2]).abs() < 1e-10
        });
        assert_eq!(
            unique.len(),
            923,
            "Found {} duplicate positions",
            923 - unique.len()
        );

        // Particle diameter check: ~3.46 nm for k=6 Au icosahedron
        let diameter = 2.0 * r_max;
        println!(
            "Au icosahedron k=6: {} atoms, diameter = {:.3} nm, "
            ,
            points.len(),
            diameter,
        );
        println!(
            "  nn distance = {:.4} nm, shell radii: {:.4}..{:.4} nm",
            nn_nm,
            nn_nm,
            r_max
        );
    }

    /// Test error handling: atom count mismatch.
    #[test]
    fn test_parse_xyz_count_mismatch() {
        let content = "5\nWrong count\nAu 0.0 0.0 0.0\nAu 1.0 0.0 0.0\n";
        let result = parse_xyz(content);
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(
            err.contains("5") && err.contains("2"),
            "Error should mention expected vs actual count: {}",
            err
        );
    }

    /// Test error handling: malformed coordinate line.
    #[test]
    fn test_parse_xyz_bad_coordinate() {
        let content = "1\nBad data\nAu 1.0 xyz 3.0\n";
        let result = parse_xyz(content);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("Invalid"));
    }

    /// Test empty lines are skipped correctly.
    #[test]
    fn test_parse_xyz_empty_lines() {
        let content = "2\nWith blanks\nAu 0.0 0.0 0.0\n\nAu 1.0 1.0 1.0\n\n";
        let points = parse_xyz(content).unwrap();
        assert_eq!(points.len(), 2);
    }
}
