//! Integration test: CDA vs Mie theory for a gold nanosphere.
//!
//! This test validates that the CDA solver reproduces the analytical Mie
//! theory extinction cross-section for a homogeneous gold sphere (R=20 nm)
//! to within an acceptable tolerance.

// TODO: Enable once the CDA solver and Mie theory are implemented.
//
// use lumina_core::mie::mie_cross_sections;
// use lumina_core::solver::cda::CdaSolver;
// use lumina_core::solver::OpticalSolver;
// use lumina_materials::johnson_christy::JohnsonChristyMaterial;
// use lumina_materials::provider::MaterialProvider;
// use lumina_geometry::primitives::{Primitive, Sphere};
// use lumina_geometry::discretise::discretise_primitive;
//
// #[test]
// fn test_cda_matches_mie_for_gold_sphere() {
//     let radius = 20.0; // nm
//     let wavelength = 520.0; // nm, near the plasmon resonance
//     let n_medium = 1.0;
//
//     // Analytical Mie result
//     let gold = JohnsonChristyMaterial::gold();
//     let eps = gold.dielectric_function(wavelength).unwrap();
//     let (mie_ext, _, _) = mie_cross_sections(radius, wavelength, eps, n_medium);
//
//     // CDA result
//     let sphere = Primitive::Sphere(Sphere {
//         centre: [0.0, 0.0, 0.0],
//         radius,
//     });
//     let lattice = discretise_primitive(&sphere, 2.0);
//     // ... convert lattice to dipoles, run CDA solver ...
//
//     // Compare: CDA should match Mie to within ~5% for adequate discretisation
//     // assert!((cda_ext - mie_ext).abs() / mie_ext < 0.05);
// }

#[test]
fn placeholder_mie_validation() {
    // Placeholder until solver is implemented.
    assert!(true);
}
