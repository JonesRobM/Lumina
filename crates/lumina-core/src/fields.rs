//! Near-field and far-field computation from solved dipole moments.
//!
//! Once the self-consistent dipole moments are known, the electric field
//! at any observation point $\mathbf{r}_{\text{obs}}$ is:
//!
//! $$
//! \mathbf{E}(\mathbf{r}_{\text{obs}}) = \mathbf{E}_{\text{inc}}(\mathbf{r}_{\text{obs}})
//! + \sum_{i=1}^{N} \mathbf{G}(\mathbf{r}_{\text{obs}}, \mathbf{r}_i) \cdot \mathbf{p}_i
//! $$

use num_complex::Complex64;

/// Compute the electric field at a single observation point due to all dipoles.
///
/// # Arguments
/// * `obs_point` - Observation position (nm).
/// * `dipole_positions` - Positions of all dipoles (nm).
/// * `dipole_moments` - Complex dipole moments, shape (N, 3).
/// * `k` - Wavenumber in the medium (nm^{-1}).
///
/// # Returns
/// The complex electric field vector [Ex, Ey, Ez] at the observation point.
pub fn field_at_point(
    _obs_point: &[f64; 3],
    _dipole_positions: &[[f64; 3]],
    _dipole_moments: &[[Complex64; 3]],
    _k: f64,
) -> [Complex64; 3] {
    // TODO: Sum Green's function contributions from each dipole.
    todo!("Field computation at observation point")
}
