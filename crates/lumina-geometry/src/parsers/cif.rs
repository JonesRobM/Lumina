//! Parser for `.cif` (Crystallographic Information File) format.
//!
//! **Status: Stub for initial implementation.**
//!
//! CIF files describe crystal structures with unit cell parameters and
//! fractional coordinates. Parsing requires handling:
//! - `_cell_length_a/b/c` and `_cell_angle_alpha/beta/gamma`
//! - `_atom_site_fract_x/y/z` or `_atom_site_Cartn_x/y/z`
//! - Symmetry operations from `_symmetry_equiv_pos_as_xyz`

use super::{ParseError, ParsedPoint};

/// Parse a CIF file from a string.
pub fn parse_cif(_content: &str) -> Result<Vec<ParsedPoint>, ParseError> {
    Err(ParseError::UnsupportedFormat(
        "CIF parser not yet implemented".to_string(),
    ))
}
