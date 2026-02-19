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
use crate::types::{Dipole, DipoleResponse, IncidentField, NearFieldMap};

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
