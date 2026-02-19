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
}
