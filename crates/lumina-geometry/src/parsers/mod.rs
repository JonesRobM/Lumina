//! File format parsers for importing molecular and mesh geometries.
//!
//! Supported formats:
//! - [`.xyz`](xyz) — XYZ molecular coordinate files
//! - [`.cif`](cif) — Crystallographic Information Files
//! - [`.obj`](obj) — Wavefront OBJ mesh files

pub mod cif;
pub mod obj;
pub mod xyz;

use thiserror::Error;

/// Errors during geometry file parsing.
#[derive(Debug, Error)]
pub enum ParseError {
    #[error("Failed to read file: {0}")]
    IoError(#[from] std::io::Error),

    #[error("Parse error at line {line}: {message}")]
    FormatError { line: usize, message: String },

    #[error("Unsupported file format: {0}")]
    UnsupportedFormat(String),
}

/// A parsed atomic/vertex position with optional metadata.
#[derive(Debug, Clone)]
pub struct ParsedPoint {
    /// Position in 3D space (angstroms for molecular files, nm for meshes).
    pub position: [f64; 3],
    /// Atom type or vertex label, if available.
    pub label: Option<String>,
}
