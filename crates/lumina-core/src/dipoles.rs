//! Construction of per-wavelength dipole polarisabilities from scene geometry.
//!
//! Every solve -- linear, SHG, THG, CLI or GUI -- needs the same pipeline: resolve
//! eps for each dipole's material, form the Clausius-Mossotti polarisability from
//! the cell volume, apply the radiative correction, and scale by the surface
//! fill fraction. This module owns that pipeline so the formula lives in exactly
//! one place.

use num_complex::Complex64;

use lumina_geometry::scene::SceneDipoles;

use crate::types::{
    clausius_mossotti, ellipsoid_polarisability_tensor, radiative_correction, Dipole,
};

/// Build the dipole polarisabilities for one wavelength.
///
/// # Arguments
/// * `geom` -- Flat scene geometry from [`SceneSpec::build_geometry`].
/// * `material_wavelength_nm` -- Wavelength at which eps is looked up *and* at which
///   the wavenumber is evaluated. For the harmonic solves this is lambda/2 or lambda/3,
///   *not* the fundamental.
/// * `environment_n` -- Refractive index of the surrounding medium. The wavenumber
///   k = 2*pi*n/lambda and the environment permittivity eps_m = n^2 are both derived from it
///   here rather than taken as parameters: a caller passing a k belonging to a
///   different harmonic than the eps lookup would be quietly, slightly wrong, and no
///   test would catch it. The environment is non-dispersive, so there is no case
///   where lambda and k legitimately disagree.
/// * `resolve_eps` -- Maps (material id, wavelength) to eps. Its error type is the
///   caller's: the CLI propagates `anyhow::Error`, the GUI a `String`, and the
///   harmonic sweeps discard it with `.ok()?` to skip a wavelength where the
///   material data does not reach.
///
/// # Surface averaging
/// The fill fraction scales the **polarisability** (alpha_SA = f * alpha_RRCM), not eps.
/// eps-mixing causes a spurious Clausius-Mossotti resonance for Drude metals at
/// f ~= 0.2-0.3.
///
/// [`SceneSpec::build_geometry`]: lumina_geometry::scene::SceneSpec::build_geometry
pub fn build_dipoles<E>(
    geom: &SceneDipoles,
    material_wavelength_nm: f64,
    environment_n: f64,
    resolve_eps: impl Fn(&str, f64) -> Result<Complex64, E>,
) -> Result<Vec<Dipole>, E> {
    let k = 2.0 * std::f64::consts::PI * environment_n / material_wavelength_nm;
    let epsilon_m = environment_n * environment_n;

    let n = geom.positions.len();
    let mut dipoles = Vec::with_capacity(n);

    for i in 0..n {
        let epsilon = resolve_eps(&geom.materials[i], material_wavelength_nm)?;
        let cell_size = geom.spacings[i];
        let volume = cell_size.powi(3);
        let f = geom.fill_fractions[i];

        let dipole = if let Some(hint) = &geom.hints[i] {
            // Anisotropic path: ellipsoidal depolarisation factors.
            let tensor = ellipsoid_polarisability_tensor(
                volume,
                hint.depol_factors,
                epsilon,
                epsilon_m,
                k,
                hint.rotation,
            );
            Dipole {
                position: geom.positions[i],
                polarisability: tensor.map(|v| v * f),
                cell_size,
            }
        } else {
            // Isotropic path: standard spherical Clausius-Mossotti.
            let alpha_cm = clausius_mossotti(volume, epsilon, epsilon_m);
            let alpha = radiative_correction(alpha_cm, k) * f;
            Dipole::isotropic(geom.positions[i], alpha, cell_size)
        };

        dipoles.push(dipole);
    }

    Ok(dipoles)
}

#[cfg(test)]
mod tests {
    use super::*;
    use lumina_geometry::scene::{DipolePolarisabilityHint, SceneDipoles};

    fn scene(fill_fractions: Vec<f64>, hints: Vec<Option<DipolePolarisabilityHint>>) -> SceneDipoles {
        let n = fill_fractions.len();
        SceneDipoles {
            positions: (0..n).map(|i| [i as f64 * 10.0, 0.0, 0.0]).collect(),
            materials: vec!["Au_JC".to_string(); n],
            spacings: vec![2.0; n],
            hints,
            fill_fractions,
        }
    }

    fn eps_const(eps: Complex64) -> impl Fn(&str, f64) -> Result<Complex64, String> {
        move |_, _| Ok(eps)
    }

    /// Interior dipole (f = 1) must equal plain RRCM Clausius-Mossotti.
    #[test]
    fn isotropic_interior_dipole_matches_rrcm() {
        let eps = Complex64::new(-10.0, 1.2);
        let out = build_dipoles(&scene(vec![1.0], vec![None]), 600.0, 1.0, eps_const(eps)).unwrap();

        let k = 2.0 * std::f64::consts::PI / 600.0;
        let expected = radiative_correction(clausius_mossotti(8.0, eps, 1.0), k);
        assert!((out[0].polarisability[0] - expected).norm() < 1e-12);
        assert_eq!(out[0].position, [0.0, 0.0, 0.0]);
    }

    /// Surface averaging is polarisability scaling, NOT epsilon mixing:
    /// halving the fill fraction must halve alpha exactly.
    #[test]
    fn fill_fraction_scales_polarisability_linearly() {
        let eps = Complex64::new(-10.0, 1.2);
        let full = build_dipoles(&scene(vec![1.0], vec![None]), 600.0, 1.0, eps_const(eps)).unwrap();
        let half = build_dipoles(&scene(vec![0.5], vec![None]), 600.0, 1.0, eps_const(eps)).unwrap();

        let ratio = half[0].polarisability[0] / full[0].polarisability[0];
        assert!(
            (ratio - Complex64::new(0.5, 0.0)).norm() < 1e-12,
            "expected exactly 0.5 (polarisability scaling), got {ratio} -- \
             an epsilon-mixing implementation would not give exactly 0.5"
        );
    }

    /// The anisotropic branch must also be fill-fraction scaled, every element.
    #[test]
    fn fill_fraction_scales_anisotropic_tensor() {
        let eps = Complex64::new(4.0, 0.5);
        let hint = DipolePolarisabilityHint {
            depol_factors: [0.2, 0.3, 0.5],
            rotation: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        };
        let full = build_dipoles(&scene(vec![1.0], vec![Some(hint.clone())]), 600.0, 1.0, eps_const(eps)).unwrap();
        let half = build_dipoles(&scene(vec![0.5], vec![Some(hint)]), 600.0, 1.0, eps_const(eps)).unwrap();

        for c in 0..9 {
            let f = full[0].polarisability[c];
            if f.norm() < 1e-30 { continue; }
            let ratio = half[0].polarisability[c] / f;
            assert!(
                (ratio - Complex64::new(0.5, 0.0)).norm() < 1e-12,
                "tensor element {c} not scaled by f: ratio {ratio}"
            );
        }
    }

    /// Each dipole must be resolved at its own material id, in order.
    #[test]
    fn resolves_each_dipole_material_in_order() {
        let mut geom = scene(vec![1.0, 1.0], vec![None, None]);
        geom.materials = vec!["Au_JC".into(), "Ag_JC".into()];
        let seen = std::cell::RefCell::new(Vec::new());
        let out = build_dipoles(&geom, 600.0, 1.0, |mat: &str, _wl| {
            seen.borrow_mut().push(mat.to_string());
            Ok::<_, String>(Complex64::new(-10.0, 1.2))
        })
        .unwrap();
        assert_eq!(out.len(), 2);
        assert_eq!(*seen.borrow(), vec!["Au_JC".to_string(), "Ag_JC".to_string()]);
    }

    /// The resolver's error is the caller's, and must propagate unchanged --
    /// this is what lets the harmonic sweeps skip an out-of-range wavelength.
    #[test]
    fn resolver_error_propagates_to_caller() {
        let err = build_dipoles(&scene(vec![1.0], vec![None]), 600.0, 1.0, |_: &str, _| {
            Err::<Complex64, _>("out of range".to_string())
        })
        .unwrap_err();
        assert_eq!(err, "out of range");
    }

    /// The material wavelength is independent of k -- the harmonic solves
    /// evaluate epsilon at lambda/2 while correcting at k(lambda/2).
    #[test]
    fn material_wavelength_is_passed_through_to_resolver() {
        let seen_wl = std::cell::Cell::new(0.0);
        build_dipoles(&scene(vec![1.0], vec![None]), 300.0, 1.0, |_: &str, wl| {
            seen_wl.set(wl);
            Ok::<_, String>(Complex64::new(4.0, 0.1))
        })
        .unwrap();
        assert_eq!(seen_wl.get(), 300.0);
    }

    /// k and eps_m are derived internally from (lambda, n_env), so they cannot be
    /// passed inconsistently. This pins the derivation: k = 2*pi*n/lambda.
    #[test]
    fn wavenumber_is_derived_from_wavelength_and_environment_index() {
        let eps = Complex64::new(-10.0, 1.2);
        let n_env = 1.5_f64;
        let lambda = 600.0_f64;
        let out = build_dipoles(&scene(vec![1.0], vec![None]), lambda, n_env, eps_const(eps)).unwrap();

        let k = 2.0 * std::f64::consts::PI * n_env / lambda;
        let expected = radiative_correction(clausius_mossotti(8.0, eps, n_env * n_env), k);
        assert!(
            (out[0].polarisability[0] - expected).norm() < 1e-12,
            "expected alpha from k = 2*pi*n/lambda with eps_m = n^2, got {:?}",
            out[0].polarisability[0]
        );
    }


    /// Every other test here uses a uniform scene, so an implementation that read
    /// `positions[0]` / `spacings[0]` / `fill_fractions[0]` / `hints[0]` for *every*
    /// dipole would pass them all -- the fixture would be degenerate in exactly the
    /// dimension under test, which is the flaw that allowed the original
    /// surface-averaging near-miss. This scene varies every per-dipole field at once
    /// and checks each dipole on its own.
    #[test]
    fn each_dipole_uses_its_own_position_spacing_fill_fraction_and_hint() {
        let eps = Complex64::new(4.0, 0.5);
        let n_env = 1.0_f64;
        let lambda = 600.0_f64;
        let k = 2.0 * std::f64::consts::PI * n_env / lambda;

        let hint = DipolePolarisabilityHint {
            depol_factors: [0.2, 0.3, 0.5],
            rotation: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        };

        let geom = SceneDipoles {
            positions: vec![[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]],
            materials: vec!["Au_JC".into(), "Ag_JC".into(), "Cu_JC".into()],
            spacings: vec![2.0, 3.0, 4.0],
            hints: vec![None, Some(hint.clone()), None],
            fill_fractions: vec![1.0, 0.5, 0.25],
        };

        let out = build_dipoles(&geom, lambda, n_env, eps_const(eps)).unwrap();
        assert_eq!(out.len(), 3);

        // Positions must be carried through per dipole, not broadcast from [0].
        assert_eq!(out[0].position, [1.0, 0.0, 0.0]);
        assert_eq!(out[1].position, [0.0, 2.0, 0.0]);
        assert_eq!(out[2].position, [0.0, 0.0, 3.0]);

        // Dipole 0: isotropic, spacing 2 (V = 8), f = 1.
        let expect0 = radiative_correction(clausius_mossotti(8.0, eps, 1.0), k);
        assert!(
            (out[0].polarisability[0] - expect0).norm() < 1e-12,
            "dipole 0 wrong: got {:?}, want {expect0:?}",
            out[0].polarisability[0]
        );

        // Dipole 2: isotropic, spacing 4 (V = 64), f = 0.25 -- a different volume
        // *and* a different fill fraction from dipole 0.
        let expect2 = radiative_correction(clausius_mossotti(64.0, eps, 1.0), k) * 0.25;
        assert!(
            (out[2].polarisability[0] - expect2).norm() < 1e-12,
            "dipole 2 wrong (spacing or fill fraction taken from another dipole): got {:?}, want {expect2:?}",
            out[2].polarisability[0]
        );

        // Dipole 1 is the only one with a hint, so it must be the only anisotropic
        // one: an ellipsoid tensor has unequal diagonal entries here, whereas the
        // isotropic path gives equal ones and exactly zero off-diagonals.
        assert!(
            (out[1].polarisability[0] - out[1].polarisability[4]).norm() > 1e-12,
            "dipole 1 should be anisotropic (its hint was ignored)"
        );
        assert!(
            (out[0].polarisability[0] - out[0].polarisability[4]).norm() < 1e-30,
            "dipole 0 should be isotropic (a hint leaked from another dipole)"
        );
        assert!(
            (out[2].polarisability[0] - out[2].polarisability[4]).norm() < 1e-30,
            "dipole 2 should be isotropic (a hint leaked from another dipole)"
        );

        // And dipole 1's anisotropic tensor must carry its own f = 0.5.
        let full = ellipsoid_polarisability_tensor(27.0, hint.depol_factors, eps, 1.0, k, hint.rotation);
        assert!(
            (out[1].polarisability[0] - full[0] * 0.5).norm() < 1e-12,
            "dipole 1 anisotropic path used the wrong spacing or fill fraction"
        );
    }


    /// FCD averages over the *source* cell, so each dipole must carry its own cell
    /// size -- a scene mixing spacings cannot be described by one scalar.
    #[test]
    fn each_dipole_carries_its_own_cell_size() {
        let geom = SceneDipoles {
            positions: vec![[0.0; 3], [10.0, 0.0, 0.0]],
            materials: vec!["Au_JC".to_string(); 2],
            spacings: vec![2.0, 3.5],
            hints: vec![None, None],
            fill_fractions: vec![1.0, 1.0],
        };
        let out = build_dipoles(&geom, 600.0, 1.0, eps_const(Complex64::new(4.0, 0.5))).unwrap();
        assert_eq!(out[0].cell_size, 2.0);
        assert_eq!(out[1].cell_size, 3.5);
    }

}
