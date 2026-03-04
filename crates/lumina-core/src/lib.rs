//! # Lumina Core
//!
//! The numerical backbone of the Lumina framework. This crate implements
//! optical solvers for computing the electromagnetic response of nanostructured
//! materials.
//!
//! ## Architecture
//!
//! All solvers implement the [`solver::OpticalSolver`] trait, which provides a
//! uniform interface for computing cross-sections, dipole moments, and
//! near-field maps. The primary implementation is the Coupled Dipole
//! Approximation ([`solver::cda::CdaSolver`]).
//!
//! ## Modules
//!
//! - [`types`] — Core data structures (dipoles, parameters, results).
//! - [`solver`] — Optical solver trait and CDA implementation.
//! - [`fields`] — Near-field and far-field computation.
//! - [`mie`] — Analytical Mie theory for validation.
//! - [`nonlinear`] — SHG/THG source terms (stub).
//! - [`solver::ewald`] — Ewald-accelerated periodic Green's function.
//! - [`time_domain`] — Time-domain field reconstruction (stub).

pub mod fields;
pub mod mie;
pub mod nonlinear;
pub mod periodic;
pub mod solver;
pub mod time_domain;
pub mod types;

// Re-export frequently used periodic-array types for convenience.
pub use solver::ewald::EwaldGreens;
pub use types::{BlochCrossSections, DispersionMap};
