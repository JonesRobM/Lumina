//! Lorentz–Drude (LD) analytic dielectric functions for metals.
//!
//! Parameters from A. D. Rakić, A. B. Djurišić, J. M. Elazar and
//! M. L. Majewski, "Optical properties of metallic films for vertical-cavity
//! optoelectronic devices", *Appl. Opt.* **37**, 5271 (1998),
//! doi:10.1364/AO.37.005271, Table 2 (LD model).
//!
//! The values below were taken from the transcription used by
//! refractiveindex.info (M. Polyanskiy, `refractiveindex.info-scripts`,
//! "Rakic 1998 - <metal> (LD model).py"), and agree with the widely used
//! `LD.m` transcription by B. Ung. The model is
//!
//! $$\varepsilon(\omega) = 1 - \frac{f_0\,\omega_p^2}{\omega(\omega + i\Gamma_0)}
//!   + \sum_{j\ge1} \frac{f_j\,\omega_p^2}{\omega_j^2 - \omega^2 - i\omega\Gamma_j}$$
//!
//! with all energies in eV. The sign convention gives Im ε > 0 for absorbing
//! media, matching the rest of Lumina.
//!
//! # Validity range
//!
//! The paper states that the fits cover 0.1–6 eV overall. Each metal's range
//! here is that interval intersected with the per-metal range used by
//! refractiveindex.info, so wavelengths outside the data the fit was made
//! against are refused with [`MaterialError::OutOfRange`].

use num_complex::Complex64;

use crate::provider::{MaterialError, MaterialProvider};

/// Photon energy (eV) × wavelength (nm): $hc$ in eV·nm.
const HC_EV_NM: f64 = 1_239.841_984;

/// One Lorentz oscillator: strength $f_j$, damping $\Gamma_j$ (eV) and
/// resonance energy $\omega_j$ (eV).
#[derive(Debug, Clone, Copy)]
struct Oscillator {
    f: f64,
    gamma_ev: f64,
    omega_ev: f64,
}

const fn osc(f: f64, gamma_ev: f64, omega_ev: f64) -> Oscillator {
    Oscillator {
        f,
        gamma_ev,
        omega_ev,
    }
}

/// A metal described by the Rakić 1998 Lorentz–Drude model.
#[derive(Debug, Clone)]
pub struct LorentzDrudeMaterial {
    name: &'static str,
    plasma_ev: f64,
    f0: f64,
    gamma0_ev: f64,
    oscillators: Vec<Oscillator>,
    /// Validity range in photon energy (eV).
    energy_range_ev: (f64, f64),
}

impl LorentzDrudeMaterial {
    /// Silver. Valid 0.1–5 eV (248–12 398 nm).
    pub fn silver() -> Self {
        Self {
            name: "Ag (Rakić 1998, Lorentz–Drude)",
            plasma_ev: 9.01,
            f0: 0.845,
            gamma0_ev: 0.048,
            oscillators: vec![
                osc(0.065, 3.886, 0.816),
                osc(0.124, 0.452, 4.481),
                osc(0.011, 0.065, 8.185),
                osc(0.840, 0.916, 9.083),
                osc(5.646, 2.419, 20.29),
            ],
            energy_range_ev: (0.1, 5.0),
        }
    }

    /// Aluminium. Valid 0.1–6 eV (207–12 398 nm).
    pub fn aluminium() -> Self {
        Self {
            name: "Al (Rakić 1998, Lorentz–Drude)",
            plasma_ev: 14.98,
            f0: 0.523,
            gamma0_ev: 0.047,
            oscillators: vec![
                osc(0.227, 0.333, 0.162),
                osc(0.050, 0.312, 1.544),
                osc(0.166, 1.351, 1.808),
                osc(0.030, 3.382, 3.473),
            ],
            energy_range_ev: (0.1, 6.0),
        }
    }

    /// Gold. Valid 0.2–5 eV (248–6 199 nm).
    pub fn gold() -> Self {
        Self {
            name: "Au (Rakić 1998, Lorentz–Drude)",
            plasma_ev: 9.03,
            f0: 0.760,
            gamma0_ev: 0.053,
            oscillators: vec![
                osc(0.024, 0.241, 0.415),
                osc(0.010, 0.345, 0.830),
                osc(0.071, 0.870, 2.969),
                osc(0.601, 2.494, 4.304),
                osc(4.384, 2.214, 13.32),
            ],
            energy_range_ev: (0.2, 5.0),
        }
    }

    /// Beryllium. Valid 0.1–5 eV (248–12 398 nm).
    pub fn beryllium() -> Self {
        Self {
            name: "Be (Rakić 1998, Lorentz–Drude)",
            plasma_ev: 18.51,
            f0: 0.084,
            gamma0_ev: 0.035,
            oscillators: vec![
                osc(0.031, 1.664, 0.100),
                osc(0.140, 3.395, 1.032),
                osc(0.530, 4.454, 3.183),
                osc(0.130, 1.802, 4.604),
            ],
            energy_range_ev: (0.1, 5.0),
        }
    }

    /// Chromium. Valid 0.1–5 eV (248–12 398 nm).
    pub fn chromium() -> Self {
        Self {
            name: "Cr (Rakić 1998, Lorentz–Drude)",
            plasma_ev: 10.75,
            f0: 0.168,
            gamma0_ev: 0.047,
            oscillators: vec![
                osc(0.151, 3.175, 0.121),
                osc(0.150, 1.305, 0.543),
                osc(1.149, 2.676, 1.970),
                osc(0.825, 1.335, 8.775),
            ],
            energy_range_ev: (0.1, 5.0),
        }
    }

    /// Copper. Valid 0.1–6 eV (207–12 398 nm).
    pub fn copper() -> Self {
        Self {
            name: "Cu (Rakić 1998, Lorentz–Drude)",
            plasma_ev: 10.83,
            f0: 0.575,
            gamma0_ev: 0.030,
            oscillators: vec![
                osc(0.061, 0.378, 0.291),
                osc(0.104, 1.056, 2.957),
                osc(0.723, 3.213, 5.300),
                osc(0.638, 4.305, 11.18),
            ],
            energy_range_ev: (0.1, 6.0),
        }
    }

    /// Nickel. Valid 0.2–5 eV (248–6 199 nm).
    pub fn nickel() -> Self {
        Self {
            name: "Ni (Rakić 1998, Lorentz–Drude)",
            plasma_ev: 15.92,
            f0: 0.096,
            gamma0_ev: 0.048,
            oscillators: vec![
                osc(0.100, 4.511, 0.174),
                osc(0.135, 1.334, 0.582),
                osc(0.106, 2.178, 1.597),
                osc(0.729, 6.292, 6.089),
            ],
            energy_range_ev: (0.2, 5.0),
        }
    }

    /// Palladium. Valid 0.1–5 eV (248–12 398 nm).
    pub fn palladium() -> Self {
        Self {
            name: "Pd (Rakić 1998, Lorentz–Drude)",
            plasma_ev: 9.72,
            f0: 0.330,
            gamma0_ev: 0.008,
            oscillators: vec![
                osc(0.649, 2.950, 0.336),
                osc(0.121, 0.555, 0.501),
                osc(0.638, 4.621, 1.659),
                osc(0.453, 3.236, 5.715),
            ],
            energy_range_ev: (0.1, 5.0),
        }
    }

    /// Platinum. Valid 0.1–5 eV (248–12 398 nm).
    // ω₃ = 3.141 eV is a fitted resonance energy from the paper, not π.
    #[allow(clippy::approx_constant)]
    pub fn platinum() -> Self {
        Self {
            name: "Pt (Rakić 1998, Lorentz–Drude)",
            plasma_ev: 9.59,
            f0: 0.333,
            gamma0_ev: 0.080,
            oscillators: vec![
                osc(0.191, 0.517, 0.780),
                osc(0.659, 1.838, 1.314),
                osc(0.547, 3.668, 3.141),
                osc(3.576, 8.517, 9.249),
            ],
            energy_range_ev: (0.1, 5.0),
        }
    }

    /// Titanium. Valid 0.1–5 eV (248–12 398 nm).
    pub fn titanium() -> Self {
        Self {
            name: "Ti (Rakić 1998, Lorentz–Drude)",
            plasma_ev: 7.29,
            f0: 0.148,
            gamma0_ev: 0.082,
            oscillators: vec![
                osc(0.899, 2.276, 0.777),
                osc(0.393, 2.518, 1.545),
                osc(0.187, 1.663, 2.509),
                osc(0.001, 1.762, 19.43),
            ],
            energy_range_ev: (0.1, 5.0),
        }
    }

    /// Tungsten. Valid 0.1–5 eV (248–12 398 nm).
    pub fn tungsten() -> Self {
        Self {
            name: "W (Rakić 1998, Lorentz–Drude)",
            plasma_ev: 13.22,
            f0: 0.206,
            gamma0_ev: 0.064,
            oscillators: vec![
                osc(0.054, 0.530, 1.004),
                osc(0.166, 1.281, 1.917),
                osc(0.706, 3.332, 3.580),
                osc(2.590, 5.836, 7.498),
            ],
            energy_range_ev: (0.1, 5.0),
        }
    }

    /// Model dielectric function at photon energy `omega_ev`, without a range check.
    fn epsilon_at_energy(&self, omega_ev: f64) -> Complex64 {
        let w = omega_ev;
        let wp2 = self.plasma_ev * self.plasma_ev;
        // Intraband (Drude) term with plasma frequency Ω_p = √f₀ ω_p.
        let mut eps =
            Complex64::new(1.0, 0.0) - self.f0 * wp2 / (w * Complex64::new(w, self.gamma0_ev));
        // Interband (Lorentz) terms.
        for o in &self.oscillators {
            eps += o.f * wp2 / Complex64::new(o.omega_ev * o.omega_ev - w * w, -w * o.gamma_ev);
        }
        eps
    }
}

impl MaterialProvider for LorentzDrudeMaterial {
    fn name(&self) -> &str {
        self.name
    }

    fn wavelength_range(&self) -> (f64, f64) {
        let (e_min, e_max) = self.energy_range_ev;
        (HC_EV_NM / e_max, HC_EV_NM / e_min)
    }

    fn dielectric_function(&self, wavelength_nm: f64) -> Result<Complex64, MaterialError> {
        let (min, max) = self.wavelength_range();
        // Compare in energy so the documented eV limits are the ones applied.
        let omega_ev = HC_EV_NM / wavelength_nm;
        // The small tolerance lets the range end points (converted to nm and
        // back) pass despite rounding; NaN fails both comparisons.
        let (e_min, e_max) = self.energy_range_ev;
        if !(omega_ev >= e_min * (1.0 - 1e-12) && omega_ev <= e_max * (1.0 + 1e-12)) {
            return Err(MaterialError::OutOfRange {
                wavelength_nm,
                min,
                max,
            });
        }
        Ok(self.epsilon_at_energy(omega_ev))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Check ε against $(n + ik)^2$ from refractiveindex.info's Rakić-LD
    /// tables (n, k given to 5 significant figures).
    fn assert_matches_nk(m: &LorentzDrudeMaterial, wavelength_nm: f64, n: f64, k: f64) {
        let expected = Complex64::new(n, k).powi(2);
        let eps = m.dielectric_function(wavelength_nm).expect("in range");
        let rel = (eps - expected).norm() / expected.norm();
        assert!(
            rel < 1e-3,
            "{} at {wavelength_nm} nm: ε = {eps}, expected {expected} (rel. err {rel:.2e})",
            m.name()
        );
    }

    #[test]
    fn matches_published_values() {
        // (material, [(λ/nm, n, k); 2]) taken from refractiveindex.info
        // database/data/main/<metal>/nk/Rakic-LD.yml.
        let cases = [
            (
                LorentzDrudeMaterial::silver(),
                [(397.46, 0.16007, 1.799), (1001.3, 0.21994, 6.4341)],
            ),
            (
                LorentzDrudeMaterial::aluminium(),
                [(401.43, 0.45683, 4.7236), (1000.5, 1.4765, 9.2887)],
            ),
            (
                LorentzDrudeMaterial::gold(),
                [(402.85, 1.4112, 1.749), (996.62, 0.25469, 5.9622)],
            ),
            (
                LorentzDrudeMaterial::beryllium(),
                [(401.08, 2.8749, 3.1161), (998.35, 3.2487, 3.8405)],
            ),
            (
                LorentzDrudeMaterial::chromium(),
                [(401.08, 1.7714, 3.5652), (998.35, 4.4537, 4.2425)],
            ),
            (
                LorentzDrudeMaterial::copper(),
                [(399.16, 1.2496, 2.2459), (1007.5, 0.25557, 6.5372)],
            ),
            (
                LorentzDrudeMaterial::nickel(),
                [(399.48, 1.6285, 2.3684), (1000.7, 2.8011, 5.0401)],
            ),
            (
                LorentzDrudeMaterial::palladium(),
                [(399.83, 1.2932, 2.9193), (999.62, 2.4421, 5.8583)],
            ),
            (
                LorentzDrudeMaterial::platinum(),
                [(399.83, 1.7237, 2.8515), (999.62, 3.4497, 5.7307)],
            ),
            (
                LorentzDrudeMaterial::titanium(),
                [(400.13, 1.0647, 2.5401), (1002.3, 2.9718, 3.8098)],
            ),
            (
                LorentzDrudeMaterial::tungsten(),
                [(399.83, 3.3153, 2.52), (999.62, 3.0308, 3.4667)],
            ),
        ];
        for (m, points) in &cases {
            for &(wl, n, k) in points {
                assert_matches_nk(m, wl, n, k);
            }
        }
    }

    #[test]
    fn approaches_drude_limit_at_low_energy() {
        // Well below every interband resonance the free-electron term dominates:
        // ε ≈ 1 − Ω_p²/(ω(ω + iΓ₀)) + Σ f_j ω_p²/ω_j² (static Lorentz background).
        // Use Ag at its long-wavelength limit (0.1 eV), where ω ≪ ω_j for j ≥ 2.
        let m = LorentzDrudeMaterial::silver();
        let w = 0.1;
        let wp2 = m.plasma_ev * m.plasma_ev;
        let drude = Complex64::new(1.0, 0.0) - m.f0 * wp2 / (w * Complex64::new(w, m.gamma0_ev));
        let eps = m.dielectric_function(HC_EV_NM / w).expect("in range");
        assert!(
            eps.re < -5000.0,
            "Ag should be strongly metallic at 0.1 eV: {eps}"
        );
        assert!(eps.im > 0.0, "Absorbing medium needs Im ε > 0: {eps}");
        let rel = (eps - drude).norm() / drude.norm();
        assert!(
            rel < 0.02,
            "Ag at 0.1 eV should be Drude-dominated: ε = {eps}, Drude = {drude}"
        );
        // Free-carrier absorption grows towards low energy.
        let eps_1ev = m.dielectric_function(HC_EV_NM).expect("in range");
        assert!(eps.im > eps_1ev.im && eps.re < eps_1ev.re);
    }

    #[test]
    fn imaginary_part_positive_across_range() {
        for m in [
            LorentzDrudeMaterial::silver(),
            LorentzDrudeMaterial::aluminium(),
            LorentzDrudeMaterial::gold(),
            LorentzDrudeMaterial::beryllium(),
            LorentzDrudeMaterial::chromium(),
            LorentzDrudeMaterial::copper(),
            LorentzDrudeMaterial::nickel(),
            LorentzDrudeMaterial::palladium(),
            LorentzDrudeMaterial::platinum(),
            LorentzDrudeMaterial::titanium(),
            LorentzDrudeMaterial::tungsten(),
        ] {
            let (lo, hi) = m.wavelength_range();
            for i in 0..=200 {
                let wl = lo * (hi / lo).powf(i as f64 / 200.0);
                let eps = m.dielectric_function(wl).expect("in range");
                assert!(eps.im > 0.0, "{} at {wl:.1} nm: ε = {eps}", m.name());
            }
        }
    }

    #[test]
    fn refuses_wavelengths_outside_fitted_range() {
        let pt = LorentzDrudeMaterial::platinum();
        let (lo, hi) = pt.wavelength_range();
        assert!(
            (lo - 247.968).abs() < 1e-3 && (hi - 12_398.42).abs() < 1e-2,
            "{lo}, {hi}"
        );
        assert!(pt.dielectric_function(lo).is_ok());
        assert!(pt.dielectric_function(hi).is_ok());
        assert!(matches!(
            pt.dielectric_function(200.0),
            Err(MaterialError::OutOfRange { .. })
        ));
        assert!(matches!(
            pt.dielectric_function(20_000.0),
            Err(MaterialError::OutOfRange { .. })
        ));
        assert!(pt.dielectric_function(f64::NAN).is_err());
    }
}
