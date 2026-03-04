//! Lattice specification for periodic nanostructure simulations.
//!
//! [`LatticeSpec`] describes a 1D chain or 2D periodic array by its
//! real-space lattice vectors. It also provides helpers for computing
//! reciprocal lattice vectors and the optimal Ewald splitting parameter
//! used by [`lumina_core::solver::ewald::EwaldGreens`].

use serde::{Deserialize, Serialize};

fn default_truncation() -> usize { 5 }

/// Lattice geometry for a 1D chain or 2D periodic array.
///
/// # TOML Examples
///
/// **2D square array, 100 nm pitch:**
/// ```toml
/// [simulation.lattice]
/// a1 = [100.0, 0.0, 0.0]
/// a2 = [0.0, 100.0, 0.0]
/// ```
///
/// **1D chain, 80 nm pitch along x:**
/// ```toml
/// [simulation.lattice]
/// a1 = [80.0, 0.0, 0.0]
/// ```
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LatticeSpec {
    /// First lattice vector in nm (required for 1D and 2D).
    pub a1: [f64; 3],
    /// Second lattice vector in nm. `None` for 1D chains.
    pub a2: Option<[f64; 3]>,
    /// Ewald splitting parameter η (nm⁻¹). Auto-computed from cell geometry
    /// if not supplied: η = √(π / A_cell) for 2D, √(π / |a₁|) for 1D.
    #[serde(default)]
    pub ewald_eta: Option<f64>,
    /// Number of real-space shells to sum over (default 5).
    #[serde(default = "default_truncation")]
    pub truncation_real: usize,
    /// Number of reciprocal-space shells to sum over (default 5).
    #[serde(default = "default_truncation")]
    pub truncation_recip: usize,
}

impl LatticeSpec {
    /// Construct a 1D chain with the given lattice vector.
    pub fn chain(a1: [f64; 3]) -> Self {
        Self {
            a1,
            a2: None,
            ewald_eta: None,
            truncation_real: 5,
            truncation_recip: 5,
        }
    }

    /// Construct a 2D array with the given lattice vectors.
    pub fn planar(a1: [f64; 3], a2: [f64; 3]) -> Self {
        Self {
            a1,
            a2: Some(a2),
            ewald_eta: None,
            truncation_real: 5,
            truncation_recip: 5,
        }
    }

    /// Returns true if this is a 2D (planar) lattice.
    pub fn is_2d(&self) -> bool {
        self.a2.is_some()
    }

    /// Unit cell area (nm²) for 2D, or |a₁| (nm) for 1D.
    pub fn cell_measure(&self) -> f64 {
        match self.a2 {
            Some(a2) => {
                // |a1 × a2|  (z-component for xy-plane lattice)
                let cx = self.a1[1] * a2[2] - self.a1[2] * a2[1];
                let cy = self.a1[2] * a2[0] - self.a1[0] * a2[2];
                let cz = self.a1[0] * a2[1] - self.a1[1] * a2[0];
                (cx * cx + cy * cy + cz * cz).sqrt()
            }
            None => {
                let [x, y, z] = self.a1;
                (x * x + y * y + z * z).sqrt()
            }
        }
    }

    /// Optimal Ewald splitting parameter η (nm⁻¹).
    ///
    /// Uses the supplied value if set, otherwise auto-computes:
    /// - 2D: η = √(π / A_cell)
    /// - 1D: η = √(π / |a₁|)
    pub fn eta(&self) -> f64 {
        if let Some(eta) = self.ewald_eta {
            return eta;
        }
        let m = self.cell_measure();
        (std::f64::consts::PI / m).sqrt()
    }

    /// Reciprocal lattice vectors for a 2D planar array.
    ///
    /// For lattice vectors **a₁**, **a₂** in the xy plane, the reciprocal
    /// vectors satisfy **b**ᵢ · **a**ⱼ = 2π δᵢⱼ:
    ///
    /// **b₁** = 2π (**a₂ × ẑ**) / (**a₁ · (a₂ × ẑ)**)
    /// **b₂** = 2π (**ẑ × a₁**) / (**a₂ · (ẑ × a₁)**)
    ///
    /// Returns `None` for 1D lattices.
    pub fn reciprocal_vectors(&self) -> Option<([f64; 3], [f64; 3])> {
        let a2 = self.a2?;
        let a1 = self.a1;
        let two_pi = 2.0 * std::f64::consts::PI;

        // Cross product a1 × a2 to get the area normal
        let cx = a1[1] * a2[2] - a1[2] * a2[1];
        let cy = a1[2] * a2[0] - a1[0] * a2[2];
        let cz = a1[0] * a2[1] - a1[1] * a2[0];
        let area = (cx * cx + cy * cy + cz * cz).sqrt();

        // b1 = 2π (a2 × ẑ) / (a1 · (a2 × ẑ))
        // For a lattice in xy: ẑ = [0,0,1]
        // a2 × ẑ = [a2y, -a2x, 0]
        let a2_x_z = [a2[1], -a2[0], 0.0];
        let denom1 = a1[0] * a2_x_z[0] + a1[1] * a2_x_z[1];
        let b1 = [
            two_pi * a2_x_z[0] / denom1,
            two_pi * a2_x_z[1] / denom1,
            0.0,
        ];

        // b2 = 2π (ẑ × a1) / (a2 · (ẑ × a1))
        // ẑ × a1 = [-a1y, a1x, 0]
        let z_x_a1 = [-a1[1], a1[0], 0.0];
        let denom2 = a2[0] * z_x_a1[0] + a2[1] * z_x_a1[1];
        let b2 = [
            two_pi * z_x_a1[0] / denom2,
            two_pi * z_x_a1[1] / denom2,
            0.0,
        ];

        let _ = area; // used implicitly via cross product
        Some((b1, b2))
    }

    /// Reciprocal lattice vector for a 1D chain.
    ///
    /// Returns **b₁** = 2π **a₁** / |**a₁**|²
    /// Returns `None` for 2D lattices.
    pub fn reciprocal_vector_1d(&self) -> Option<[f64; 3]> {
        if self.a2.is_some() {
            return None;
        }
        let [x, y, z] = self.a1;
        let len_sq = x * x + y * y + z * z;
        let scale = 2.0 * std::f64::consts::PI / len_sq;
        Some([x * scale, y * scale, z * scale])
    }

    /// Generate the real-space lattice points within `n_shells` shells.
    ///
    /// For 2D: returns all **R** = m₁**a₁** + m₂**a₂** with |m₁|,|m₂| ≤ n_shells.
    /// For 1D: returns all **R** = m₁**a₁** with |m₁| ≤ n_shells.
    /// The origin (m₁ = m₂ = 0) **is** included (caller filters it if needed).
    pub fn real_lattice_points(&self, n_shells: usize) -> Vec<[f64; 3]> {
        let ns = n_shells as i64;
        match self.a2 {
            Some(a2) => {
                let mut pts = Vec::with_capacity((2 * ns as usize + 1).pow(2));
                for m1 in -ns..=ns {
                    for m2 in -ns..=ns {
                        pts.push([
                            m1 as f64 * self.a1[0] + m2 as f64 * a2[0],
                            m1 as f64 * self.a1[1] + m2 as f64 * a2[1],
                            m1 as f64 * self.a1[2] + m2 as f64 * a2[2],
                        ]);
                    }
                }
                pts
            }
            None => {
                let mut pts = Vec::with_capacity(2 * ns as usize + 1);
                for m1 in -ns..=ns {
                    pts.push([
                        m1 as f64 * self.a1[0],
                        m1 as f64 * self.a1[1],
                        m1 as f64 * self.a1[2],
                    ]);
                }
                pts
            }
        }
    }

    /// Generate the reciprocal-space lattice points within `n_shells` shells.
    ///
    /// For 2D: returns all **G** = n₁**b₁** + n₂**b₂** with |n₁|,|n₂| ≤ n_shells.
    /// For 1D: returns all **G** = n₁**b₁** with |n₁| ≤ n_shells.
    /// The origin **is** included.
    pub fn recip_lattice_points(&self, n_shells: usize) -> Vec<[f64; 3]> {
        let ns = n_shells as i64;
        match self.reciprocal_vectors() {
            Some((b1, b2)) => {
                let mut pts = Vec::with_capacity((2 * ns as usize + 1).pow(2));
                for n1 in -ns..=ns {
                    for n2 in -ns..=ns {
                        pts.push([
                            n1 as f64 * b1[0] + n2 as f64 * b2[0],
                            n1 as f64 * b1[1] + n2 as f64 * b2[1],
                            0.0,
                        ]);
                    }
                }
                pts
            }
            None => {
                // 1D
                let b1 = self.reciprocal_vector_1d().unwrap_or([0.0; 3]);
                let mut pts = Vec::with_capacity(2 * ns as usize + 1);
                for n1 in -ns..=ns {
                    pts.push([
                        n1 as f64 * b1[0],
                        n1 as f64 * b1[1],
                        n1 as f64 * b1[2],
                    ]);
                }
                pts
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn near(a: f64, b: f64) -> bool { (a - b).abs() < 1e-10 }

    #[test]
    fn test_square_lattice_reciprocal() {
        let lat = LatticeSpec::planar([100.0, 0.0, 0.0], [0.0, 100.0, 0.0]);
        let (b1, b2) = lat.reciprocal_vectors().unwrap();
        let two_pi = 2.0 * std::f64::consts::PI;
        assert!(near(b1[0], two_pi / 100.0));
        assert!(near(b1[1], 0.0));
        assert!(near(b2[0], 0.0));
        assert!(near(b2[1], two_pi / 100.0));
    }

    #[test]
    fn test_hexagonal_lattice_reciprocal() {
        let a = 100.0_f64;
        let a2 = [a / 2.0, a * (3.0_f64).sqrt() / 2.0, 0.0];
        let lat = LatticeSpec::planar([a, 0.0, 0.0], a2);
        let (b1, b2) = lat.reciprocal_vectors().unwrap();
        let two_pi = 2.0 * std::f64::consts::PI;
        assert!(near(b1[0] * a, two_pi));
        assert!(near(b1[0] * a2[0] + b1[1] * a2[1], 0.0));
        assert!(near(b2[0] * a, 0.0));
        assert!(near(b2[0] * a2[0] + b2[1] * a2[1], two_pi));
    }

    #[test]
    fn test_cell_area_square() {
        let lat = LatticeSpec::planar([100.0, 0.0, 0.0], [0.0, 100.0, 0.0]);
        assert!(near(lat.cell_measure(), 10000.0));
    }

    #[test]
    fn test_eta_auto() {
        let lat = LatticeSpec::planar([100.0, 0.0, 0.0], [0.0, 100.0, 0.0]);
        let expected = (std::f64::consts::PI / 10000.0_f64).sqrt();
        assert!(near(lat.eta(), expected));
    }

    #[test]
    fn test_1d_chain() {
        let lat = LatticeSpec::chain([80.0, 0.0, 0.0]);
        assert!(!lat.is_2d());
        assert!(lat.reciprocal_vectors().is_none());
        let b1 = lat.reciprocal_vector_1d().unwrap();
        let two_pi = 2.0 * std::f64::consts::PI;
        assert!(near(b1[0], two_pi / 80.0));
    }

    #[test]
    fn test_real_lattice_points_count_2d() {
        let lat = LatticeSpec::planar([100.0, 0.0, 0.0], [0.0, 100.0, 0.0]);
        let pts = lat.real_lattice_points(2);
        assert_eq!(pts.len(), 25); // (2*2+1)^2
    }

    #[test]
    fn test_real_lattice_points_count_1d() {
        let lat = LatticeSpec::chain([80.0, 0.0, 0.0]);
        let pts = lat.real_lattice_points(3);
        assert_eq!(pts.len(), 7); // 2*3+1
    }
}
