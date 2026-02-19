//! Parser for Wavefront `.obj` mesh files.
//!
//! **Status: Stub for initial implementation.**
//!
//! OBJ files define 3D meshes as collections of vertices (`v`) and faces (`f`).
//! The mesh surface is used to define the boundary of a nanostructure, which
//! is then filled with dipoles by the discretisation module.

use super::{ParseError, ParsedPoint};

/// Parse an OBJ file, extracting vertex positions.
pub fn parse_obj(_content: &str) -> Result<Vec<ParsedPoint>, ParseError> {
    // TODO: Parse vertices and faces. Faces are needed for the
    // inside/outside test during discretisation.
    todo!("OBJ parser")
}
