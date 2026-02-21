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

use crate::solver::cda::greens::dyadic_greens_tensor;
use crate::solver::{NearFieldPlane, SolverError};
use crate::types::{Dipole, DipoleResponse, FarFieldMap, IncidentField, NearFieldMap};

/// Compute the electric field at a single observation point due to all dipoles.
pub fn field_at_point(
    obs_point: &[f64; 3],
    dipoles: &[Dipole],
    response: &DipoleResponse,
    k: f64,
    incident: &IncidentField,
) -> [Complex64; 3] {
    // Start with the incident field
    let mut e = incident.at_position(obs_point, k);

    // Add contribution from each dipole via the Green's function
    for (i, dipole) in dipoles.iter().enumerate() {
        // Skip if observation point is too close to the dipole
        let dx = obs_point[0] - dipole.position[0];
        let dy = obs_point[1] - dipole.position[1];
        let dz = obs_point[2] - dipole.position[2];
        let r_sq = dx * dx + dy * dy + dz * dz;
        if r_sq < 1e-20 {
            continue;
        }

        let g = dyadic_greens_tensor(obs_point, &dipole.position, k);
        for a in 0..3 {
            for b in 0..3 {
                e[a] += g[[a, b]] * response.moments[[i, b]];
            }
        }
    }

    e
}

/// Compute a 2D near-field intensity map on an observation plane.
pub fn compute_near_field_map(
    dipoles: &[Dipole],
    response: &DipoleResponse,
    plane: &NearFieldPlane,
    k: f64,
    incident: &IncidentField,
) -> Result<NearFieldMap, SolverError> {
    let nx = plane.nx;
    let ny = plane.ny;

    // Build local coordinate system on the plane
    let (u_hat, v_hat) = plane_basis(&plane.normal);

    let mut positions = Vec::with_capacity(nx * ny);
    let mut field_intensity = Vec::with_capacity(nx * ny);

    let x_min = -plane.half_width;
    let x_max = plane.half_width;
    let y_min = -plane.half_height;
    let y_max = plane.half_height;

    let dx = if nx > 1 { (x_max - x_min) / (nx - 1) as f64 } else { 0.0 };
    let dy = if ny > 1 { (y_max - y_min) / (ny - 1) as f64 } else { 0.0 };

    for iy in 0..ny {
        for ix in 0..nx {
            let u = x_min + ix as f64 * dx;
            let v = y_min + iy as f64 * dy;

            let obs = [
                plane.centre[0] + u * u_hat[0] + v * v_hat[0],
                plane.centre[1] + u * u_hat[1] + v * v_hat[1],
                plane.centre[2] + u * u_hat[2] + v * v_hat[2],
            ];

            let e = field_at_point(&obs, dipoles, response, k, incident);
            let intensity = e[0].norm_sqr() + e[1].norm_sqr() + e[2].norm_sqr();

            positions.push(obs);
            field_intensity.push(intensity);
        }
    }

    Ok(NearFieldMap {
        positions,
        field_intensity,
        nx,
        ny,
        extent: [
            plane.centre[0] + x_min * u_hat[0] + y_min * v_hat[0],
            plane.centre[0] + x_max * u_hat[0] + y_max * v_hat[0],
            plane.centre[1] + x_min * u_hat[1] + y_min * v_hat[1],
            plane.centre[1] + x_max * u_hat[1] + y_max * v_hat[1],
        ],
    })
}

/// Construct an orthonormal basis (u, v) for a plane with the given normal.
fn plane_basis(normal: &[f64; 3]) -> ([f64; 3], [f64; 3]) {
    let n = normalise(normal);

    // Choose a vector not parallel to n
    let seed = if n[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };

    // u = normalise(seed - (seed . n) * n)
    let dot = seed[0] * n[0] + seed[1] * n[1] + seed[2] * n[2];
    let u_raw = [
        seed[0] - dot * n[0],
        seed[1] - dot * n[1],
        seed[2] - dot * n[2],
    ];
    let u = normalise(&u_raw);

    // v = n x u
    let v = [
        n[1] * u[2] - n[2] * u[1],
        n[2] * u[0] - n[0] * u[2],
        n[0] * u[1] - n[1] * u[0],
    ];

    (u, v)
}

fn normalise(v: &[f64; 3]) -> [f64; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    [v[0] / len, v[1] / len, v[2] / len]
}

/// Compute the far-field radiation pattern on a spherical grid.
///
/// Samples the differential scattering intensity
/// $\frac{dC_\text{sca}}{d\Omega} \propto |(\mathbf{I} - \hat{\mathbf{r}}\hat{\mathbf{r}}) \cdot \mathbf{F}(\hat{\mathbf{r}})|^2$
/// at `n_theta × n_phi` directions, where
/// $\mathbf{F}(\hat{\mathbf{r}}) = \sum_i \mathbf{p}_i e^{-ik\hat{\mathbf{r}}\cdot\mathbf{r}_i}$.
///
/// # Arguments
/// * `dipoles` — Dipole positions.
/// * `response` — Solved dipole moments at a single wavelength.
/// * `k` — Wavenumber in the medium (nm⁻¹).
/// * `n_theta` — Number of polar angle samples in [0, π].
/// * `n_phi` — Number of azimuthal angle samples in [0, 2π).
pub fn compute_far_field(
    dipoles: &[Dipole],
    response: &DipoleResponse,
    k: f64,
    n_theta: usize,
    n_phi: usize,
) -> FarFieldMap {
    use std::f64::consts::PI;

    let mut theta_vals = Vec::with_capacity(n_theta * n_phi);
    let mut phi_vals = Vec::with_capacity(n_theta * n_phi);
    let mut intensity = Vec::with_capacity(n_theta * n_phi);

    for it in 0..n_theta {
        let theta = PI * it as f64 / (n_theta - 1).max(1) as f64;
        let sin_t = theta.sin();
        let cos_t = theta.cos();

        for ip in 0..n_phi {
            let phi = 2.0 * PI * ip as f64 / n_phi as f64;
            let sin_p = phi.sin();
            let cos_p = phi.cos();

            // Unit direction vector r̂
            let rhat = [sin_t * cos_p, sin_t * sin_p, cos_t];

            // Structure factor F(r̂) = Σᵢ pᵢ exp(-ik r̂·rᵢ)
            let mut f = [Complex64::new(0.0, 0.0); 3];
            for (i, dipole) in dipoles.iter().enumerate() {
                let kdotr = k * (
                    rhat[0] * dipole.position[0]
                    + rhat[1] * dipole.position[1]
                    + rhat[2] * dipole.position[2]
                );
                let phase = Complex64::new(0.0, -kdotr).exp();
                for a in 0..3 {
                    f[a] += response.moments[[i, a]] * phase;
                }
            }

            // Transverse projection: (I - r̂r̂) · F
            // (f_transverse)_a = f_a - (r̂ · f) r̂_a
            let rdotf = rhat[0] * f[0] + rhat[1] * f[1] + rhat[2] * f[2];
            let ft = [
                f[0] - rdotf * rhat[0],
                f[1] - rdotf * rhat[1],
                f[2] - rdotf * rhat[2],
            ];

            let intens = ft[0].norm_sqr() + ft[1].norm_sqr() + ft[2].norm_sqr();

            theta_vals.push(theta);
            phi_vals.push(phi);
            intensity.push(intens);
        }
    }

    FarFieldMap {
        wavelength_nm: response.wavelength_nm,
        theta: theta_vals,
        phi: phi_vals,
        intensity,
        n_theta,
        n_phi,
    }
}

/// Compute the circular dichroism ΔC_ext = C_ext(LCP) − C_ext(RCP).
///
/// Uses the two-response formula (Draine & Goodman approach extended for CD):
/// $\Delta C_\text{ext} = \frac{k}{|E_0|^2} \sum_i
///   \bigl[\operatorname{Im}(p_{x,i}^{(y)}) - \operatorname{Im}(p_{y,i}^{(x)})\bigr]$
///
/// where $p_{x,i}^{(y)}$ is the y-component of dipole $i$'s moment when the incident
/// field is x-polarised, and vice versa.
///
/// # Arguments
/// * `response_x` — Solved dipole moments for x-polarised incidence.
/// * `response_y` — Solved dipole moments for y-polarised incidence.
/// * `k` — Wavenumber in the medium (nm⁻¹).
pub fn compute_circular_dichroism(
    response_x: &DipoleResponse,
    response_y: &DipoleResponse,
    k: f64,
) -> f64 {
    let n = response_x.moments.nrows();
    let mut delta: f64 = 0.0;

    for i in 0..n {
        // y-component of dipole i under x-pol excitation
        let py_under_x = response_x.moments[[i, 1]].im;
        // x-component of dipole i under y-pol excitation
        let px_under_y = response_y.moments[[i, 0]].im;
        delta += py_under_x - px_under_y;
    }

    k * delta
}
