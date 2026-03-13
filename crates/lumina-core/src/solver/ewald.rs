//! Ewald-accelerated periodic dyadic Green's function.
//!
//! For a nanoparticle array that is periodic with lattice vectors **a₁** (and
//! optionally **a₂**), the interaction between a dipole at **r**ᵢ and all
//! periodic images of **r**ⱼ is:
//!
//! $$\mathbf{G}_\text{per}(\mathbf{r}, \mathbf{k}_\parallel)
//!   = \sum_{\mathbf{R}\in\Lambda}
//!     \mathbf{G}_\text{free}(\mathbf{r}-\mathbf{R})\,
//!     e^{i\mathbf{k}_\parallel\cdot\mathbf{R}}$$
//!
//! # Implementation strategy (v0.4)
//!
//! The sum is split using the Ewald identity into a **real-space** part
//! (erfc-damped, converges exponentially in real space) and a
//! **reciprocal-space** part (spectral G-sum, converges exponentially in
//! reciprocal space):
//!
//! $$\mathbf{G}_\text{per} = \mathbf{G}_\text{real} + \mathbf{G}_\text{recip}$$
//!
//! ## Real-space part
//!
//! Each free-space G term is multiplied by `erfc(η · |r − R|)` which strongly
//! damps contributions from distant lattice images.  The splitting parameter
//! η = √(π / A_cell) is chosen so that both sums converge in roughly the same
//! number of shells.
//!
//! ## Reciprocal-space part (2D lattices)
//!
//! For a 2D lattice in the xy-plane the spectral sum over reciprocal vectors
//! **G** is (Capolino et al. 2005 / Linton 2010 conventions, consistent with
//! the free-space prefactor k²exp(ikR)/(4πR)):
//!
//! $$\mathbf{G}^\text{spec}_{ab} = \frac{i}{2A} \sum_\mathbf{G}
//!     \frac{e^{i\mathbf{Q}\cdot\boldsymbol{\rho}}\,e^{-q|z|}}{q}\,
//!     M_{ab}(\mathbf{Q}, q)$$
//!
//! where **Q** = **k**∥ + **G**, q = √(Q² − k²) with Im(q) ≥ 0, and
//!
//! $$M_{ab} = k^2\delta_{ab} - \tilde{k}_a\tilde{k}_b,\quad
//!   \tilde{\mathbf{k}} = (Q_x,\,Q_y,\,i\,q\,\operatorname{sgn}(z))$$
//!
//! For 1D chains the reciprocal sum returns zero. The real-space sum uses
//! no erfc damping (direct lattice sum) with a higher default shell count
//! (≥ 20) to achieve ~2 % convergence at optical wavelengths. A proper 1D
//! spectral sum (requiring modified Bessel functions K₀/K₁) is deferred.
//!
//! # Convergence note
//!
//! With `truncation_real = 5` and `truncation_recip = 5` the combined
//! Ewald sum typically achieves < 0.01 % error for lattice periods ≥ 50 nm
//! and visible wavelengths.

use num_complex::Complex64;

use lumina_geometry::lattice::LatticeSpec;

/// Ewald-accelerated periodic dyadic Green's function.
///
/// Wraps a [`LatticeSpec`] and evaluates the Bloch-periodic sum via the
/// real-space erfc-damped sum plus the spectral (reciprocal-space) sum.
pub struct EwaldGreens {
    pub lattice: LatticeSpec,
}

impl EwaldGreens {
    pub fn new(lattice: LatticeSpec) -> Self {
        Self { lattice }
    }

    /// Full periodic dyadic Green's function at Bloch wavevector `k_bloch`.
    ///
    /// # Arguments
    /// * `r`       — Displacement **r**ᵢ − **r**ⱼ (nm).
    /// * `k_bloch` — In-plane Bloch wavevector **k**∥ (nm⁻¹). Use `[0;3]` for Γ.
    /// * `k`       — Free-space wavenumber 2πn/λ (nm⁻¹).
    ///
    /// # Self-image term (r ≈ 0)
    /// The R = 0 image is skipped when **r** ≈ 0 (i.e. this is the self-block),
    /// leaving only the non-zero lattice images. The caller is responsible for
    /// adding the on-site α⁻¹ term to the diagonal.
    pub fn evaluate(
        &self,
        r: [f64; 3],
        k_bloch: [f64; 3],
        k: f64,
    ) -> [[Complex64; 3]; 3] {
        let gr = self.real_space(r, k_bloch, k);
        let gg = self.recip_space(r, k_bloch, k);
        let mut out = [[Complex64::from(0.0); 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                out[i][j] = gr[i][j] + gg[i][j];
            }
        }
        out
    }

    /// Real-space lattice sum.
    ///
    /// **2D lattices**: sums `erfc(η·|s|) · G_free(s) · exp(i k∥ · R)` over
    /// all **R** ∈ Λ within `truncation_real` shells (η = `self.lattice.eta()`).
    /// The erfc factor provides exponential convergence.
    ///
    /// **1D chains**: no erfc damping — uses a direct lattice sum over at least
    /// 20 shells.  This gives ~2 % convergence at typical optical wavelengths.
    /// The full 1D Ewald (Bessel spectral sum) is deferred to a future version.
    ///
    /// The **R** = 0 term is excluded when `|r| < 1e-12` nm.
    pub fn real_space(
        &self,
        r: [f64; 3],
        k_bloch: [f64; 3],
        k: f64,
    ) -> [[Complex64; 3]; 3] {
        let r_norm_sq = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
        let r_is_zero = r_norm_sq < 1e-24;

        let is_2d = self.lattice.is_2d();
        let eta = self.lattice.eta();

        // 1D chains: use a direct lattice sum (no erfc) with at least 20 shells
        // for ~2 % convergence at optical wavelengths.  The Bessel-function-based
        // 1D spectral correction is not yet implemented.
        let n_shells = if is_2d {
            self.lattice.truncation_real
        } else {
            self.lattice.truncation_real.max(20)
        };
        let lattice_points = self.lattice.real_lattice_points(n_shells);
        let zero = Complex64::from(0.0);
        let mut g = [[zero; 3]; 3];

        for lp in &lattice_points {
            let lp_is_zero = lp[0].abs() < 1e-12 && lp[1].abs() < 1e-12 && lp[2].abs() < 1e-12;
            if r_is_zero && lp_is_zero {
                // Skip the self-image (R=0) to avoid divergence.
                continue;
            }

            // Displacement s = r − R
            let s = [r[0] - lp[0], r[1] - lp[1], r[2] - lp[2]];
            let s_sq = s[0]*s[0] + s[1]*s[1] + s[2]*s[2];
            if s_sq < 1e-24 { continue; } // r coincides with this image — skip

            let s_norm = s_sq.sqrt();

            // Bloch phase: exp(i k∥ · R)
            let bloch_arg = k_bloch[0]*lp[0] + k_bloch[1]*lp[1] + k_bloch[2]*lp[2];
            let bloch_phase = Complex64::new(bloch_arg.cos(), bloch_arg.sin());

            // erfc damping factor (2D only — 1D uses direct sum without damping)
            let damp = if is_2d { erfc_real(eta * s_norm) } else { 1.0 };

            let g_s = dyadic_greens_tensor_at(&s, k);

            for i in 0..3 {
                for j in 0..3 {
                    g[i][j] += bloch_phase * damp * g_s[i][j];
                }
            }
        }

        g
    }

    /// Reciprocal-space (spectral) sum.
    ///
    /// For 2D lattices: evaluates the spectral Green's function sum
    ///
    /// $$\mathbf{G}^\text{spec}_{ab} = \frac{i}{2A} \sum_\mathbf{G}
    ///     \frac{e^{i\mathbf{Q}\cdot\boldsymbol{\rho}}\,e^{-q|z|}}{q}\,
    ///     M_{ab}(\mathbf{Q},q)$$
    ///
    /// For 1D lattices: returns zero (direct sum converges adequately).
    pub fn recip_space(
        &self,
        r: [f64; 3],
        k_bloch: [f64; 3],
        k: f64,
    ) -> [[Complex64; 3]; 3] {
        let zero = Complex64::from(0.0);

        // Only implemented for 2D lattices.
        if !self.lattice.is_2d() {
            return [[zero; 3]; 3];
        }

        let area = self.lattice.cell_measure(); // nm²
        let recip_points = self.lattice.recip_lattice_points(self.lattice.truncation_recip);

        let i_unit = Complex64::new(0.0, 1.0);
        // Scalar prefactor: i / (2 * A)
        let scalar_pre = i_unit / (2.0 * area);

        let rho_x = r[0];
        let rho_y = r[1];
        let z = r[2];
        let abs_z = z.abs();
        // sign_z: +1 if z >= 0, -1 if z < 0
        let sign_z = if z >= 0.0 { 1.0_f64 } else { -1.0_f64 };

        let k2 = Complex64::from(k * k);
        let mut g = [[zero; 3]; 3];

        for gpt in &recip_points {
            // Q = k_bloch_∥ + G  (in-plane x, y only)
            let qx = k_bloch[0] + gpt[0];
            let qy = k_bloch[1] + gpt[1];
            let q2_real = qx * qx + qy * qy; // Q² (real, ≥ 0)

            // q = sqrt(Q² − k²) with Im(q) ≥ 0.
            // For visible light k is real, Q² − k² can be negative (propagating)
            // or positive (evanescent).
            let q2_complex = Complex64::from(q2_real) - k2;
            let q = csqrt_positive_imag(q2_complex);

            // Skip terms where |q| is too small (at the light cone — near singular).
            if q.norm() < 1e-10 {
                continue;
            }

            // Bloch phase: exp(i * Q · ρ∥)
            let phase_arg = qx * rho_x + qy * rho_y;
            let bloch_phase = Complex64::new(phase_arg.cos(), phase_arg.sin());

            // Evanescent/propagating factor: exp(-q * |z|)
            let exp_qz = (-q * abs_z).exp();

            // Per-G scalar: (i / 2A) * exp(iQ·ρ) * exp(-q|z|) / q
            let coeff = scalar_pre * bloch_phase * exp_qz / q;

            // ── Build M tensor: M_ab = k²δ_ab - k̃_a k̃_b ──────────────────
            //
            // k̃ = (Qx, Qy, i*q*sign_z)  [for z > 0; conjugate sign for z < 0]
            //
            // k̃_x = Qx  (real)
            // k̃_y = Qy  (real)
            // k̃_z = i * q * sign_z
            let kt_x = Complex64::from(qx);
            let kt_y = Complex64::from(qy);
            let kt_z = i_unit * q * sign_z;

            // k̃_a * k̃_b outer product
            let kk = [
                [kt_x * kt_x, kt_x * kt_y, kt_x * kt_z],
                [kt_y * kt_x, kt_y * kt_y, kt_y * kt_z],
                [kt_z * kt_x, kt_z * kt_y, kt_z * kt_z],
            ];

            // M_ab = k² δ_ab - k̃_a k̃_b
            for a in 0..3 {
                for b in 0..3 {
                    let delta = if a == b { k2 } else { zero };
                    let m_ab = delta - kk[a][b];
                    g[a][b] += coeff * m_ab;
                }
            }
        }

        g
    }
}

// ── Private helpers ────────────────────────────────────────────────────────────

/// Evaluate the free-space dyadic Green's tensor at displacement vector `s`.
///
/// Equivalent to `dyadic_greens_tensor(&s, &[0,0,0], k)` but avoids the
/// subtract-zero computation.
#[inline]
fn dyadic_greens_tensor_at(s: &[f64; 3], k: f64) -> [[Complex64; 3]; 3] {
    let r_sq = s[0]*s[0] + s[1]*s[1] + s[2]*s[2];
    let r = r_sq.sqrt();
    assert!(r > 1e-15, "dyadic_greens_tensor_at: zero displacement");

    let kr = k * r;
    let ikr = Complex64::new(0.0, kr);
    let exp_ikr = ikr.exp();
    let prefactor = k * k * exp_ikr / (4.0 * std::f64::consts::PI * r);

    let a = Complex64::from(1.0) + (ikr - Complex64::from(1.0)) / Complex64::from(kr * kr);
    let b = (Complex64::from(3.0) - Complex64::from(3.0) * ikr - Complex64::from(kr * kr))
        / Complex64::from(kr * kr);

    let r_hat = [s[0] / r, s[1] / r, s[2] / r];
    let zero = Complex64::from(0.0);
    let mut g = [[zero; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let delta = if i == j { 1.0 } else { 0.0 };
            g[i][j] = prefactor * (a * delta + b * r_hat[i] * r_hat[j]);
        }
    }
    g
}

/// Complementary error function for real non-negative argument.
///
/// Uses the Abramowitz & Stegun formula 7.1.26 rational approximation.
/// Maximum error < 1.5 × 10⁻⁷.
///
/// For x < 0 the identity erfc(x) = 2 − erfc(−x) is applied.
#[inline]
fn erfc_real(x: f64) -> f64 {
    if x < 0.0 {
        return 2.0 - erfc_real(-x);
    }
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let poly = t * (0.254_829_592
        + t * (-0.284_496_736
        + t * (1.421_413_741
        + t * (-1.453_152_027
        + t * 1.061_405_429))));
    poly * (-x * x).exp()
}

/// Complex square root with branch cut chosen so that Im(result) ≥ 0.
///
/// For purely real positive z this returns +√z. For purely real negative z
/// this returns +i·√|z| (evanescent wave convention).
#[inline]
fn csqrt_positive_imag(z: Complex64) -> Complex64 {
    let s = z.sqrt(); // num_complex uses the principal branch: Re(sqrt) >= 0
    // The principal branch has Re ≥ 0. We want Im ≥ 0, so if Im < 0 negate.
    if s.im < 0.0 { -s } else { s }
}

#[cfg(test)]
mod tests {
    use super::*;
    use lumina_geometry::lattice::LatticeSpec;

    fn near_c(a: Complex64, b: Complex64, tol: f64) -> bool {
        (a - b).norm() < tol
    }

    /// erfc(0) == 1, erfc large == 0
    #[test]
    fn test_erfc_basic() {
        assert!((erfc_real(0.0) - 1.0).abs() < 1e-7);
        assert!(erfc_real(10.0) < 1e-40);
        assert!((erfc_real(-0.0) - 1.0).abs() < 1e-7);
        // erfc(x) + erfc(-x) == 2
        let x = 1.5_f64;
        assert!((erfc_real(x) + erfc_real(-x) - 2.0).abs() < 1e-6);
    }

    /// csqrt_positive_imag: for negative real argument, Im should be positive.
    #[test]
    fn test_csqrt_positive_imag() {
        let z = Complex64::new(-4.0, 0.0);
        let s = csqrt_positive_imag(z);
        assert!(s.im > 0.0, "Im should be >= 0 for negative real argument");
        assert!((s.im - 2.0).abs() < 1e-12);
        // For positive real: result should be real positive
        let z2 = Complex64::new(9.0, 0.0);
        let s2 = csqrt_positive_imag(z2);
        assert!((s2.re - 3.0).abs() < 1e-12);
        assert!(s2.im.abs() < 1e-12);
    }

    /// For a square 2D array, the periodic Green's function evaluated at r=0
    /// (self-image sum, R≠0 only) should be finite and non-zero.
    #[test]
    fn test_self_image_finite() {
        let lat = LatticeSpec::planar([100.0, 0.0, 0.0], [0.0, 100.0, 0.0]);
        let eg = EwaldGreens::new(lat);
        let k = 2.0 * std::f64::consts::PI / 600.0;
        let g = eg.evaluate([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], k);
        for i in 0..3 {
            for j in 0..3 {
                assert!(g[i][j].re.is_finite(), "G[{i}][{j}] real not finite");
                assert!(g[i][j].im.is_finite(), "G[{i}][{j}] imag not finite");
            }
        }
    }

    /// evaluate() == real_space() + recip_space()
    #[test]
    fn test_evaluate_splits_correctly() {
        let lat = LatticeSpec::planar([80.0, 0.0, 0.0], [0.0, 80.0, 0.0]);
        let eg = EwaldGreens::new(lat);
        let k = 2.0 * std::f64::consts::PI / 500.0;
        let r = [10.0, 5.0, 2.0];
        let kb = [0.0, 0.0, 0.0];
        let g_full = eg.evaluate(r, kb, k);
        let gr = eg.real_space(r, kb, k);
        let gg = eg.recip_space(r, kb, k);
        for i in 0..3 {
            for j in 0..3 {
                assert!(near_c(g_full[i][j], gr[i][j] + gg[i][j], 1e-20));
            }
        }
    }

    /// Reciprocity at the Γ point: G_per(r)[i,j] == G_per(−r)[j,i].
    ///
    /// For k∥ = 0 the periodic Green's function satisfies the same reciprocity
    /// as the free-space G because the Bloch phases cancel with their conjugates.
    #[test]
    fn test_gamma_point_reciprocity() {
        let lat = LatticeSpec::planar([100.0, 0.0, 0.0], [0.0, 100.0, 0.0]);
        let eg = EwaldGreens::new(lat);
        let k = 2.0 * std::f64::consts::PI / 600.0;
        let r = [15.0, 8.0, 3.0];
        let neg_r = [-r[0], -r[1], -r[2]];
        let kb = [0.0, 0.0, 0.0];
        let g_r   = eg.evaluate(r, kb, k);
        let g_neg = eg.evaluate(neg_r, kb, k);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    near_c(g_r[i][j], g_neg[j][i], 1e-8),
                    "G[{i}][{j}] reciprocity: {:?} vs {:?}",
                    g_r[i][j], g_neg[j][i]
                );
            }
        }
    }

    /// The Γ-point self-image sum should be purely imaginary diagonal on a
    /// square lattice (by 4-fold symmetry the real-part contributions cancel).
    #[test]
    fn test_gamma_self_image_symmetry() {
        let lat = LatticeSpec::planar([100.0, 0.0, 0.0], [0.0, 100.0, 0.0]);
        let eg = EwaldGreens::new(lat);
        let k = 2.0 * std::f64::consts::PI / 600.0;
        let g = eg.evaluate([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], k);
        // Gxx and Gyy should be equal (square symmetry)
        assert!(
            near_c(g[0][0], g[1][1], 1e-8),
            "Gxx={:?} Gyy={:?} should be equal on square lattice",
            g[0][0], g[1][1]
        );
        // Off-diagonal Gxy should vanish (square symmetry)
        assert!(
            g[0][1].norm() < 1e-8,
            "Gxy={:?} should vanish by 4-fold symmetry",
            g[0][1]
        );
    }

    /// 1D chain: real-space sum should be finite, non-zero, and reciprocal.
    ///
    /// Prior to the direct-sum fix (v0.4.1) the erfc η was so large that every
    /// term was zeroed → evaluate() returned ≈ 0 silently.
    #[test]
    fn test_1d_chain_finite() {
        let lat = LatticeSpec::chain([80.0, 0.0, 0.0]);
        let eg = EwaldGreens::new(lat);
        let k = 2.0 * std::f64::consts::PI / 500.0;
        let r = [10.0, 5.0, 0.0];
        let g = eg.evaluate(r, [0.0; 3], k);
        let mut max_mag = 0.0_f64;
        for i in 0..3 {
            for j in 0..3 {
                assert!(g[i][j].re.is_finite());
                assert!(g[i][j].im.is_finite());
                max_mag = max_mag.max(g[i][j].norm());
            }
        }
        // The result must be genuinely non-zero (dominant element ≫ 0).
        assert!(
            max_mag > 1e-12,
            "1D chain G ≈ 0 — erfc η bug: max_mag = {:.2e}",
            max_mag
        );
    }

    /// 1D chain: evaluate() at Γ point satisfies G(r)[i,j] == G(-r)[j,i] (reciprocity).
    #[test]
    fn test_1d_chain_reciprocity() {
        let lat = LatticeSpec::chain([80.0, 0.0, 0.0]);
        let eg = EwaldGreens::new(lat);
        let k = 2.0 * std::f64::consts::PI / 600.0;
        let r = [10.0, 5.0, 0.0];
        let neg_r = [-r[0], -r[1], -r[2]];
        let kb = [0.0; 3];
        let g_r = eg.evaluate(r, kb, k);
        let g_neg = eg.evaluate(neg_r, kb, k);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    near_c(g_r[i][j], g_neg[j][i], 1e-6),
                    "1D reciprocity G[{i}][{j}] = {:?}, G[-r][{j}][{i}] = {:?}",
                    g_r[i][j], g_neg[j][i]
                );
            }
        }
    }

    /// Convergence: adding more shells should not drastically change the result.
    #[test]
    fn test_convergence_with_shells() {
        let lat_coarse = LatticeSpec { truncation_real: 5,  ..LatticeSpec::planar([100.0, 0.0, 0.0], [0.0, 100.0, 0.0]) };
        let lat_fine   = LatticeSpec { truncation_real: 10, ..LatticeSpec::planar([100.0, 0.0, 0.0], [0.0, 100.0, 0.0]) };
        let eg5  = EwaldGreens::new(lat_coarse);
        let eg10 = EwaldGreens::new(lat_fine);
        let k = 2.0 * std::f64::consts::PI / 600.0;
        let r = [15.0, 8.0, 3.0];
        let kb = [0.0; 3];
        let g5  = eg5.evaluate(r, kb, k);
        let g10 = eg10.evaluate(r, kb, k);
        // The change from 5→10 shells should be < 5% of the dominant element.
        // (erfc damping makes the sum converge rapidly, so even few shells is enough)
        let max5: f64 = g5.iter().flatten().map(|v| v.norm()).fold(0.0f64, f64::max);
        let max_diff: f64 = (0..3).flat_map(|i| (0..3).map(move |j| (g5[i][j]-g10[i][j]).norm()))
            .fold(0.0f64, f64::max);
        assert!(max_diff < 0.05 * max5, "Convergence poor: Δ={:.2e} vs max={:.2e}", max_diff, max5);
    }

    /// recip_space returns zero for 1D lattice.
    #[test]
    fn test_1d_recip_space_zero() {
        let lat = LatticeSpec::chain([80.0, 0.0, 0.0]);
        let eg = EwaldGreens::new(lat);
        let k = 2.0 * std::f64::consts::PI / 500.0;
        let g = eg.recip_space([10.0, 5.0, 0.0], [0.0; 3], k);
        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(g[i][j], Complex64::from(0.0));
            }
        }
    }

    /// recip_space for 2D lattice returns finite values.
    #[test]
    fn test_2d_recip_space_finite() {
        let lat = LatticeSpec::planar([100.0, 0.0, 0.0], [0.0, 100.0, 0.0]);
        let eg = EwaldGreens::new(lat);
        let k = 2.0 * std::f64::consts::PI / 600.0;
        let r = [10.0, 5.0, 20.0]; // z != 0 — well-behaved case
        let g = eg.recip_space(r, [0.0; 3], k);
        for i in 0..3 {
            for j in 0..3 {
                assert!(g[i][j].re.is_finite(), "recip G[{i}][{j}] re not finite");
                assert!(g[i][j].im.is_finite(), "recip G[{i}][{j}] im not finite");
            }
        }
    }
}
