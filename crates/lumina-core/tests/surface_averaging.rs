//! Validation: surface averaging (polarisability scaling) improves CDA accuracy.
//!
//! # Method
//!
//! Surface dipoles — those whose cubic cell straddles the particle boundary — have
//! their polarisability scaled by their fill fraction f ∈ (0, 1):
//!
//!   α_SA = f · radiative_correction(clausius_mossotti(V, ε, ε_m), k)
//!
//! Interior dipoles (f = 1) are unchanged.
//!
//! # Scope
//!
//! The test uses a TiO₂-like dielectric (ε = 4 + 0.5i) where the improvement is
//! monotonic and predictable. For Drude metals (Au at λ > 550 nm), polarisability
//! scaling reduces staircase artefacts asymptotically (as d → 0) but does not
//! guarantee strict improvement at a single coarse spacing due to non-monotonic
//! coupling effects.
//!
//! # Engineering note
//!
//! Polarisability scaling is not the Evans & Stephens χ-mixing approach described in
//! Yurkin & Hoekstra (2007). However, χ-mixing passes through the CM resonance
//! (ε_eff + 2ε_m = 0) for Drude metals at fill fractions f ≈ 0.2–0.3. No production
//! DDA code (ADDA, DDSCAT, OpenDDA) implements surface averaging; they use binary
//! inclusion/exclusion. Polarisability scaling is a deliberate engineering choice for
//! LuminaCDA that avoids pathological resonances for metallic nanostructures.

use lumina_core::mie::mie_cross_sections;
use lumina_core::solver::cda::CdaSolver;
use lumina_core::solver::OpticalSolver;
use lumina_core::types::{clausius_mossotti, radiative_correction, Dipole, SimulationParams};
use lumina_geometry::scene::{ObjectSpec, ObjectTransform, SceneSpec, ShapeSpec};
use num_complex::Complex64;

fn sphere_geom(radius_nm: f64, spacing_nm: f64) -> lumina_geometry::scene::SceneDipoles {
    SceneSpec {
        objects: vec![ObjectSpec {
            name: "sphere".into(),
            shape: ShapeSpec::Sphere { radius: radius_nm },
            material: "test".into(),
            dipole_spacing: spacing_nm,
            transform: ObjectTransform::default(),
            species_map: Default::default(),
        }],
    }
    .build_geometry()
    .expect("sphere geometry should build")
}

fn build_dipoles_scaled(
    geom: &lumina_geometry::scene::SceneDipoles,
    fill_fractions: &[f64],
    epsilon: Complex64,
    epsilon_m: f64,
    k: f64,
) -> Vec<Dipole> {
    geom.positions
        .iter()
        .zip(geom.spacings.iter())
        .zip(fill_fractions.iter())
        .map(|((pos, &spacing), &f)| {
            let volume = spacing.powi(3);
            let alpha = radiative_correction(clausius_mossotti(volume, epsilon, epsilon_m), k) * f;
            Dipole::isotropic(*pos, alpha)
        })
        .collect()
}

/// Verify the geometry pipeline populates fill fractions correctly.
#[test]
fn test_fill_fractions_are_populated() {
    let geom = sphere_geom(10.0, 3.0);
    assert_eq!(
        geom.fill_fractions.len(),
        geom.positions.len(),
        "fill_fractions must be parallel to positions"
    );
    for &f in &geom.fill_fractions {
        assert!(f > 0.0 && f <= 1.0, "fill_fraction out of range: {f}");
    }
    assert!(
        geom.fill_fractions.iter().any(|&f| f < 1.0),
        "Expected some surface dipoles with fill_fraction < 1"
    );
}

/// Verify that surface averaging reduces C_ext toward the Mie value for a
/// dielectric sphere where the no-SA case overestimates due to coarse discretisation.
///
/// TiO₂-like: ε = 4 + 0.5i, R = 10 nm, d = 3 nm, λ = 600 nm.
/// At this spacing the coarse cubic lattice over-counts volume (N·d³/V ≈ 1.10),
/// so no-SA overestimates C_ext. SA corrects surface dipoles downward.
#[test]
fn test_surface_averaging_improves_dielectric_sphere() {
    let radius_nm = 10.0_f64;
    let spacing_nm = 3.0_f64;
    let wavelength_nm = 600.0_f64;
    let epsilon_m = 1.0_f64;
    let epsilon = Complex64::new(4.0, 0.5); // TiO₂-like dielectric

    let k = 2.0 * std::f64::consts::PI * epsilon_m.sqrt() / wavelength_nm;

    let (c_mie, _, _) = mie_cross_sections(radius_nm, wavelength_nm, epsilon, epsilon_m.sqrt());

    let geom = sphere_geom(radius_nm, spacing_nm);

    let fracs_no_sa = vec![1.0_f64; geom.positions.len()];
    let dipoles_no_sa = build_dipoles_scaled(&geom, &fracs_no_sa, epsilon, epsilon_m, k);
    let dipoles_sa = build_dipoles_scaled(&geom, &geom.fill_fractions, epsilon, epsilon_m, k);

    let solver = CdaSolver::default();
    let params = SimulationParams {
        wavelength_range: [wavelength_nm, wavelength_nm],
        num_wavelengths: 1,
        environment_n: epsilon_m.sqrt(),
        ..Default::default()
    };

    let c_no_sa = solver
        .compute_cross_sections(&dipoles_no_sa, wavelength_nm, &params)
        .expect("no-SA solve")
        .extinction;
    let c_sa = solver
        .compute_cross_sections(&dipoles_sa, wavelength_nm, &params)
        .expect("SA solve")
        .extinction;

    let err_no_sa = (c_no_sa - c_mie).abs() / c_mie;
    let err_sa = (c_sa - c_mie).abs() / c_mie;

    println!(
        "Mie: {:.3e} nm²  |  no-SA: {:.3e} nm² ({:.1}%)  |  SA: {:.3e} nm² ({:.1}%)",
        c_mie, c_no_sa, err_no_sa * 100.0, c_sa, err_sa * 100.0,
    );

    // SA must be strictly closer to Mie
    assert!(
        err_sa < err_no_sa,
        "SA did not improve accuracy: {:.1}% → {:.1}%",
        err_no_sa * 100.0,
        err_sa * 100.0,
    );

    // Residual SA error < 30%
    assert!(err_sa < 0.30, "SA error {:.1}% exceeds 30%", err_sa * 100.0);

    // SA must not zero out polarisabilities
    assert!(c_sa > 0.1 * c_mie, "SA produced near-zero C_ext: {:.3e}", c_sa);
}
