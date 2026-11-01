//! Shared helpers for the `lumina-core` integration tests.
//!
//! The dipole-build pipeline lives in [`lumina_core::build_dipoles`] (debt item D1).
//! These tests go through it rather than re-implementing it: when each test hand-rolled
//! `radiative_correction(clausius_mossotti(..), k) * f` for itself, the test baseline
//! silently drifted from production -- surface averaging landed and the tests kept
//! building every dipole with fill fraction 1, so the published accuracy figures
//! described a configuration nobody ran.

use num_complex::Complex64;

use lumina_core::build_dipoles;
use lumina_core::types::Dipole;
use lumina_geometry::discretise::LatticePoint;
use lumina_geometry::scene::SceneDipoles;

/// Build dipoles from a discretised lattice, giving every dipole the same eps, through
/// the production path.
///
/// `wavelength_nm` sets the wavenumber (k = 2*pi*n/lambda) exactly as production does; eps is
/// supplied directly instead of being resolved from a material table, which is why the
/// resolver here is infallible.
///
/// Uniform `spacing_nm` and no anisotropy hints match how these tests discretise a
/// single primitive. Per-dipole spacing or hints need a `SceneDipoles` built by hand.
pub fn dipoles_with_constant_eps(
    lattice: &[LatticePoint],
    spacing_nm: f64,
    wavelength_nm: f64,
    environment_n: f64,
    epsilon: Complex64,
) -> Vec<Dipole> {
    let geom = SceneDipoles {
        positions: lattice.iter().map(|p| p.position).collect(),
        materials: vec![String::new(); lattice.len()],
        spacings: vec![spacing_nm; lattice.len()],
        hints: vec![None; lattice.len()],
        fill_fractions: lattice.iter().map(|p| p.fill_fraction).collect(),
    };

    build_dipoles(&geom, wavelength_nm, environment_n, |_, _| {
        Ok::<_, std::convert::Infallible>(epsilon)
    })
    .expect("a constant-eps resolver cannot fail")
}
