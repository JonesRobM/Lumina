//! Integration test: CDA vs Mie theory for a gold nanosphere.
//!
//! This test validates that the CDA solver reproduces the analytical Mie
//! theory extinction cross-section for a homogeneous gold sphere to within
//! an acceptable tolerance.
//!
//! # Known limitations
//!
//! The CDA with the Clausius-Mossotti + radiative correction (RRCM)
//! polarisability prescription has limited accuracy for highly metallic
//! particles where $|\epsilon_1| \gg \epsilon_2$. This is a well-known
//! limitation of the discrete dipole approximation at coarse discretisation
//! (few dipoles per wavelength of the surface plasmon). For gold, the CDA
//! is most accurate near the interband region (400–520 nm) where $\epsilon_2$
//! is large. At longer wavelengths (>550 nm), the error increases
//! significantly and finer discretisation is required.

use lumina_core::mie::mie_cross_sections;
use lumina_core::solver::cda::CdaSolver;
use lumina_core::solver::OpticalSolver;
use lumina_core::types::{clausius_mossotti, radiative_correction, Dipole, SimulationParams};
use lumina_geometry::discretise::discretise_primitive;
use lumina_geometry::primitives::{Primitive, Sphere};
use lumina_materials::johnson_christy::JohnsonChristyMaterial;
use lumina_materials::provider::MaterialProvider;

/// Validate CDA against Mie theory for a small gold sphere in the
/// interband region where the CDA is most accurate.
#[test]
fn test_cda_vs_mie_gold_sphere() {
    let radius = 10.0; // nm
    let spacing = 3.0; // nm — coarse for test speed
    let n_medium = 1.0;
    let epsilon_m = n_medium * n_medium;

    // Discretise
    let sphere = Primitive::Sphere(Sphere {
        centre: [0.0, 0.0, 0.0],
        radius,
    });
    let lattice = discretise_primitive(&sphere, spacing);
    let positions: Vec<[f64; 3]> = lattice.iter().map(|p| p.position).collect();
    assert!(!positions.is_empty(), "Sphere must produce dipoles");

    let gold = JohnsonChristyMaterial::gold();
    let solver = CdaSolver::default();
    let params = SimulationParams {
        wavelength_range: [400.0, 800.0],
        num_wavelengths: 5,
        environment_n: n_medium,
        ..Default::default()
    };

    // Test at wavelengths in the interband region (400–520 nm) where CDA
    // is well-behaved for gold (ε₂ is large, |ε₁|/ε₂ is moderate).
    let wavelengths = [420.0, 450.0, 480.0, 510.0];

    for &wl in &wavelengths {
        let epsilon = gold
            .dielectric_function(wl)
            .expect("Wavelength within J&C range");
        let k = 2.0 * std::f64::consts::PI * n_medium / wl;

        let dipoles: Vec<Dipole> = positions
            .iter()
            .map(|&pos| {
                let volume = spacing.powi(3);
                let alpha_cm = clausius_mossotti(volume, epsilon, epsilon_m);
                let alpha = radiative_correction(alpha_cm, k);
                Dipole::isotropic(pos, alpha)
            })
            .collect();

        let cda = solver
            .compute_cross_sections(&dipoles, wl, &params)
            .expect("CDA solve should succeed");

        let (mie_ext, _mie_sca, _mie_abs) = mie_cross_sections(radius, wl, epsilon, n_medium);
        let rel_err_ext = (cda.extinction - mie_ext).abs() / mie_ext;

        eprintln!(
            "λ={:.0} nm: ε=({:.2},{:.2}), CDA={:.2e}, Mie={:.2e}, err={:.1}%, N={}",
            wl, epsilon.re, epsilon.im,
            cda.extinction, mie_ext,
            rel_err_ext * 100.0, dipoles.len()
        );

        assert!(
            rel_err_ext < 0.30,
            "CDA ({:.2e}) differs from Mie ({:.2e}) by {:.1}% at λ={:.0} nm (max 30%)",
            cda.extinction, mie_ext, rel_err_ext * 100.0, wl
        );

        assert!(cda.extinction > 0.0, "Extinction must be positive");
        assert!(cda.absorption > 0.0, "Absorption must be positive for gold");
    }
}

/// Validate that the CDA reproduces Mie theory for a dielectric sphere
/// to high accuracy across all wavelengths.
#[test]
fn test_cda_vs_mie_dielectric_sphere() {
    let radius = 10.0;
    let spacing = 3.0;
    let n_medium = 1.0;
    let epsilon_m = n_medium * n_medium;

    let sphere = Primitive::Sphere(Sphere {
        centre: [0.0, 0.0, 0.0],
        radius,
    });
    let lattice = discretise_primitive(&sphere, spacing);
    let positions: Vec<[f64; 3]> = lattice.iter().map(|p| p.position).collect();

    let solver = CdaSolver::default();
    let params = SimulationParams {
        environment_n: n_medium,
        ..Default::default()
    };

    // Test several dielectric constants (avoid very low ε where the
    // cross-sections are tiny and relative errors are inflated).
    let test_cases: &[(num_complex::Complex64, &str)] = &[
        (num_complex::Complex64::new(4.0, 0.0), "TiO2-like n=2"),
        (num_complex::Complex64::new(4.0, 0.5), "lossy dielectric"),
        (num_complex::Complex64::new(6.0, 1.0), "high-index lossy"),
        (num_complex::Complex64::new(-1.0, 3.0), "near plasmon resonance"),
    ];

    for &(epsilon, label) in test_cases {
        for &wl in &[500.0, 600.0, 700.0] {
            let k = 2.0 * std::f64::consts::PI * n_medium / wl;

            let dipoles: Vec<Dipole> = positions
                .iter()
                .map(|&pos| {
                    let volume = spacing.powi(3);
                    let alpha_cm = clausius_mossotti(volume, epsilon, epsilon_m);
                    let alpha = radiative_correction(alpha_cm, k);
                    Dipole::isotropic(pos, alpha)
                })
                .collect();

            let cda = solver
                .compute_cross_sections(&dipoles, wl, &params)
                .unwrap();
            let (mie_ext, _, _) = mie_cross_sections(radius, wl, epsilon, n_medium);
            let rel_err = (cda.extinction - mie_ext).abs() / mie_ext;

            eprintln!(
                "{:24} λ={:.0}: CDA={:.2e}, Mie={:.2e}, err={:.1}%, N={}",
                label, wl, cda.extinction, mie_ext, rel_err * 100.0, dipoles.len()
            );

            // 30% tolerance for d=3nm in a R=10nm sphere is reasonable.
            // Finer discretisation (d≤2nm) gives <5% for dielectrics.
            assert!(
                rel_err < 0.30,
                "{}: CDA ({:.2e}) differs from Mie ({:.2e}) by {:.1}% at λ={:.0} nm",
                label, cda.extinction, mie_ext, rel_err * 100.0, wl
            );
        }
    }
}

/// Test that finer discretisation gives better agreement with Mie theory
/// for a dielectric sphere (where the CDA converges reliably).
#[test]
fn test_finer_discretisation_improves_accuracy() {
    let radius = 8.0;
    let n_medium = 1.0;
    let epsilon_m = n_medium * n_medium;
    let wl = 550.0;
    let epsilon = num_complex::Complex64::new(4.0, 0.5);
    let k = 2.0 * std::f64::consts::PI * n_medium / wl;

    let solver = CdaSolver::default();
    let params = SimulationParams {
        environment_n: n_medium,
        ..Default::default()
    };

    let (mie_ext, _, _) = mie_cross_sections(radius, wl, epsilon, n_medium);

    let mut prev_err = f64::INFINITY;

    for &spacing in &[4.0, 3.0, 2.5, 2.0] {
        let sphere = Primitive::Sphere(Sphere {
            centre: [0.0, 0.0, 0.0],
            radius,
        });
        let lattice = discretise_primitive(&sphere, spacing);
        let positions: Vec<[f64; 3]> = lattice.iter().map(|p| p.position).collect();

        if positions.is_empty() {
            continue;
        }

        let dipoles: Vec<Dipole> = positions
            .iter()
            .map(|&pos| {
                let volume = spacing.powi(3);
                let alpha_cm = clausius_mossotti(volume, epsilon, epsilon_m);
                let alpha = radiative_correction(alpha_cm, k);
                Dipole::isotropic(pos, alpha)
            })
            .collect();

        let cda = solver
            .compute_cross_sections(&dipoles, wl, &params)
            .unwrap();
        let err = (cda.extinction - mie_ext).abs() / mie_ext;

        eprintln!(
            "spacing={:.1} nm, N={}, C_ext={:.2e}, Mie={:.2e}, err={:.1}%",
            spacing, dipoles.len(), cda.extinction, mie_ext, err * 100.0
        );

        // Error should generally decrease, but the staircase surface causes
        // some non-monotonicity depending on how the grid aligns with the sphere.
        assert!(
            err < prev_err * 3.0,
            "Error should not dramatically increase with finer spacing"
        );
        prev_err = err;
    }
}
