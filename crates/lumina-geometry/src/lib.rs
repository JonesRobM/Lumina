//! # LuminaCDA Geometry
//!
//! Geometry handling for the LuminaCDA framework. This crate provides:
//!
//! - **Parametric primitives** ([`primitives`]) — Spheres, cylinders, cuboids,
//!   helices, and ellipsoids defined by simple parameters.
//! - **Discretisation** ([`discretise`]) — Converts shapes into cubic dipole
//!   lattices at a given spacing.
//! - **File parsers** ([`parsers`]) — Import geometries from `.xyz`, `.cif`,
//!   and `.obj` files.
//! - **Transformations** ([`transform`]) — Scale, rotate, translate, and
//!   mirror operations.

pub mod discretise;
pub mod parsers;
pub mod primitives;
pub mod transform;
