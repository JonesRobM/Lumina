//! Substrate image-dipole correction for the CDA.
//!
//! Models the reflection of dipole fields from a planar, non-magnetic substrate
//! at z = z_interface using the quasi-static image-dipole approximation.
//!
//! # Physical model
//!
//! For a nanostructure above a substrate with dielectric function ε_sub(ω),
//! each real dipole **p**_j at position **r**_j induces an image dipole at the
//! mirror position **r**_j′ = (x_j, y_j, 2z_interface − z_j):
//!
//! **p**_j′ = Δε · T · **p**_j,  T = diag(−1, −1, +1)
//!
//! where Δε = (ε_sub − ε_env)/(ε_sub + ε_env) is the quasi-static Fresnel
//! factor. The image contribution to the interaction matrix is:
//!
//! A_ij += −Δε · T · **G**(r_i, r_j′)
//!
//! (including i = j, the self-image term).
//!
//! # Limitations
//!
//! This is the quasi-static approximation. Full retarded half-space Green's
//! functions (Sommerfeld integrals) are deferred to v0.4.

use num_complex::Complex64;
use serde::{Deserialize, Serialize};

/// Specifies a planar substrate below the nanostructure.
///
/// All dipoles must satisfy `z > z_interface` for the image-dipole formula to
/// be physically valid.
///
/// # TOML example
/// ```toml
/// [simulation.substrate]
/// z_interface = 0.0
/// material    = "SiO2_Palik"
/// ```
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SubstrateSpec {
    /// z-coordinate of the substrate–medium interface in nm.
    /// Defaults to 0.0 (interface at the origin).
    #[serde(default)]
    pub z_interface: f64,
    /// Material identifier for the substrate. Resolved via the same material
    /// providers used for particle materials (e.g. `"SiO2_Palik"`).
    pub material: String,
}

/// Compute the quasi-static Fresnel reflection factor.
///
/// Δε = (ε_sub − ε_env) / (ε_sub + ε_env)
///
/// This factor scales the image-dipole Green's function contribution. It is
/// real and positive for a dense substrate (ε_sub > ε_env), indicating an
/// in-phase reflection that enhances the local field near the surface.
///
/// # Arguments
/// * `eps_sub`  — Substrate dielectric function at the current wavelength.
/// * `eps_env`  — Environment dielectric constant ε_env = n_env².
pub fn fresnel_delta_eps(eps_sub: Complex64, eps_env: f64) -> Complex64 {
    let env = Complex64::from(eps_env);
    (eps_sub - env) / (eps_sub + env)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fresnel_delta_eps_glass_in_air() {
        // SiO2 in air: ε_sub ≈ 2.25, ε_env = 1.0 → Δε = 1.25/3.25 ≈ 0.385
        let eps_sub = Complex64::new(2.25, 0.0);
        let delta = fresnel_delta_eps(eps_sub, 1.0);
        let expected = 1.25 / 3.25;
        assert!(
            (delta.re - expected).abs() < 1e-10,
            "Δε.re = {:.4}, expected {:.4}",
            delta.re, expected
        );
        assert!(delta.im.abs() < 1e-12, "Δε for lossless dielectric should be real");
    }

    #[test]
    fn test_fresnel_delta_eps_same_medium() {
        // Substrate and environment identical → Δε = 0
        let eps = Complex64::new(2.0, 0.0);
        let delta = fresnel_delta_eps(eps, 2.0);
        assert!(delta.norm() < 1e-12, "Δε should be 0 when substrate = environment");
    }

    #[test]
    fn test_fresnel_delta_eps_metallic_substrate() {
        // Metal substrate: large |ε_sub| → Δε → 1
        let eps_sub = Complex64::new(-100.0, 5.0); // metal-like
        let delta = fresnel_delta_eps(eps_sub, 1.0);
        // As ε_sub → ∞, Δε → 1
        assert!(delta.norm() > 0.95, "metallic substrate should give |Δε| close to 1");
    }
}
