//! # Lumina Materials
//!
//! Material property providers for the Lumina framework. All materials
//! implement the [`MaterialProvider`](provider::MaterialProvider) trait,
//! which provides frequency-dependent complex dielectric functions.
//!
//! ## Available data sources
//!
//! | Source | Module | Status |
//! |--------|--------|--------|
//! | Johnson & Christy (Au/Ag/Cu) | [`johnson_christy`] | Implemented |
//! | Palik handbook (TiO₂, SiO₂) | [`palik`] | Implemented |
//! | DFT (VASP/Gaussian) | [`dft`] | Stub |
//!
//! ## Interpolation
//!
//! Tabulated data is interpolated using natural cubic splines
//! ([`spline::CubicSpline`]) to provide smooth, continuous material
//! properties at arbitrary wavelengths within the data range.

pub mod dft;
pub mod johnson_christy;
pub mod palik;
pub mod provider;
pub mod spline;
