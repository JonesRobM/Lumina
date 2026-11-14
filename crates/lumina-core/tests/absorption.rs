//! Absorption cross-section: radiative-reaction term and energy balance.
//!
//! With α in SI/ε₀ units (α = 4π α_Gauss), Draine's absorption formula
//! C_abs = k Σ [Im(p·(α⁻¹)*·p*) − k³/(6π)|p|²]. A lossless scatterer must
//! therefore have C_abs = 0 and C_sca = C_ext.

use lumina_core::mie::mie_cross_sections;
use lumina_core::solver::cda::CdaSolver;
use lumina_core::solver::OpticalSolver;
use lumina_core::types::{clausius_mossotti, radiative_correction, Dipole, SimulationParams};
use lumina_geometry::discretise::discretise_primitive;
use lumina_geometry::primitives::{Primitive, Sphere};
use num_complex::Complex64;

fn params(wavelength_nm: f64) -> SimulationParams {
    SimulationParams {
        wavelength_range: [wavelength_nm, wavelength_nm],
        num_wavelengths: 1,
        environment_n: 1.0,
        ..Default::default()
    }
}

#[test]
fn lossless_single_dipole_has_zero_absorption() {
    let wl = 500.0;
    let k = 2.0 * std::f64::consts::PI / wl;
    let alpha = radiative_correction(clausius_mossotti(1000.0, Complex64::new(4.0, 0.0), 1.0), k);
    let dipoles = vec![Dipole::isotropic([0.0, 0.0, 0.0], alpha)];
    let cs = CdaSolver::default()
        .compute_cross_sections(&dipoles, wl, &params(wl))
        .expect("single-dipole solve");

    let c_sca_expected = k.powi(4) * alpha.norm_sqr() / (6.0 * std::f64::consts::PI);
    assert!(
        cs.absorption.abs() < 1e-9 * cs.extinction,
        "lossless dipole absorbs: C_abs = {:e}, C_ext = {:e}",
        cs.absorption,
        cs.extinction
    );
    assert!(
        (cs.scattering - c_sca_expected).abs() < 1e-9 * c_sca_expected,
        "C_sca = {:e}, expected k⁴|α|²/(6π) = {:e}",
        cs.scattering,
        c_sca_expected
    );
}

#[test]
fn lossless_sphere_has_zero_absorption_and_mie_scattering() {
    // ε = 4 sphere, r = 10 nm, 2 nm spacing (~520 dipoles, direct solve).
    let (wl, radius, spacing) = (500.0, 10.0, 2.0);
    let eps = Complex64::new(4.0, 0.0);
    let k = 2.0 * std::f64::consts::PI / wl;
    let sphere = Primitive::Sphere(Sphere {
        centre: [0.0; 3],
        radius,
    });
    let dipoles: Vec<Dipole> = discretise_primitive(&sphere, spacing)
        .iter()
        .map(|p| {
            let alpha = radiative_correction(clausius_mossotti(spacing.powi(3), eps, 1.0), k);
            Dipole::isotropic(p.position, alpha)
        })
        .collect();
    let solver = CdaSolver::with_fcd(1000, true, spacing);
    let cs = solver
        .compute_cross_sections(&dipoles, wl, &params(wl))
        .expect("sphere solve");

    assert!(
        cs.absorption.abs() < 1e-6 * cs.extinction,
        "lossless sphere absorbs: C_abs = {:e}, C_ext = {:e}",
        cs.absorption,
        cs.extinction
    );
    let (_, mie_sca, _) = mie_cross_sections(radius, wl, eps, 1.0);
    let err = (cs.scattering - mie_sca).abs() / mie_sca;
    assert!(
        err < 0.15,
        "C_sca = {:e} vs Mie {:e} ({:.1}%)",
        cs.scattering,
        mie_sca,
        err * 100.0
    );
}
