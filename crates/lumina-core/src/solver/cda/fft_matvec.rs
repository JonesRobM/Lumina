//! FFT-accelerated block-Toeplitz matvec for the CDA.
//!
//! When dipoles lie on a regular cubic grid and no periodic BC, substrate,
//! or FCD corrections are active, the interaction matrix is 3-level block-
//! Toeplitz. Embedding in a (2Nx × 2Ny × 2Nz) circulant reduces the matvec
//! from O(N²) to O(N log N) via 3D FFT convolution.

use ndarray::{Array1, ArrayView1};
use num_complex::Complex64;
use rustfft::FftPlanner;

use crate::types::Dipole;
use super::assembly::invert_3x3_pub;
use super::greens::dyadic_greens_tensor;

/// Check whether `dipoles` lie on a regular cubic grid with spacing `cell_size`.
///
/// Returns `Some((Nx, Ny, Nz))` if all dipoles are grid-aligned, `None` otherwise.
pub fn detect_regular_grid(dipoles: &[Dipole], cell_size: f64) -> Option<(usize, usize, usize)> {
    if dipoles.is_empty() {
        return None;
    }

    // Handle single dipole edge case.
    if dipoles.len() == 1 {
        return Some((1, 1, 1));
    }

    // Find bounding box.
    let mut x_min = dipoles[0].position[0];
    let mut y_min = dipoles[0].position[1];
    let mut z_min = dipoles[0].position[2];
    let mut x_max = x_min;
    let mut y_max = y_min;
    let mut z_max = z_min;

    for d in dipoles.iter() {
        x_min = x_min.min(d.position[0]);
        y_min = y_min.min(d.position[1]);
        z_min = z_min.min(d.position[2]);
        x_max = x_max.max(d.position[0]);
        y_max = y_max.max(d.position[1]);
        z_max = z_max.max(d.position[2]);
    }

    // Compute grid dimensions.
    let nx = ((x_max - x_min) / cell_size).round() as usize + 1;
    let ny = ((y_max - y_min) / cell_size).round() as usize + 1;
    let nz = ((z_max - z_min) / cell_size).round() as usize + 1;

    // Verify all dipoles lie on grid points.
    for d in dipoles.iter() {
        let fx = (d.position[0] - x_min) / cell_size;
        let fy = (d.position[1] - y_min) / cell_size;
        let fz = (d.position[2] - z_min) / cell_size;

        if (fx - fx.round()).abs() >= 1e-6
            || (fy - fy.round()).abs() >= 1e-6
            || (fz - fz.round()).abs() >= 1e-6
        {
            return None;
        }
    }

    Some((nx, ny, nz))
}

/// Pre-computed plan for an FFT-accelerated block-Toeplitz matvec.
pub struct FftMatvecPlan {
    nx: usize,
    ny: usize,
    nz: usize,
    /// Flat index into the (2Nx × 2Ny × 2Nz) padded grid for each dipole.
    grid_index: Vec<usize>,
    /// Precomputed FFT of 9 kernel components G_αβ (α,β ∈ {0,1,2}).
    /// Shape: 9 × fft_len where fft_len = 2Nx · 2Ny · 2Nz.
    kernel_fft: Vec<Vec<rustfft::num_complex::Complex<f64>>>,
    fft_len: usize,
}

impl FftMatvecPlan {
    /// Attempt to build an FFT matvec plan for the given dipole set.
    ///
    /// Returns `None` if dipoles are not on a regular cubic grid.
    pub fn try_build(dipoles: &[Dipole], cell_size: f64, k: f64) -> Option<Self> {
        let (nx, ny, nz) = detect_regular_grid(dipoles, cell_size)?;

        let pad_nx = 2 * nx;
        let pad_ny = 2 * ny;
        let pad_nz = 2 * nz;
        let fft_len = pad_nx * pad_ny * pad_nz;

        // Find bounding box origin to map dipole positions to grid indices.
        let x_min = dipoles.iter().map(|d| d.position[0]).fold(f64::INFINITY, f64::min);
        let y_min = dipoles.iter().map(|d| d.position[1]).fold(f64::INFINITY, f64::min);
        let z_min = dipoles.iter().map(|d| d.position[2]).fold(f64::INFINITY, f64::min);

        // Map each dipole to its flat index in the padded grid.
        let mut grid_index = Vec::with_capacity(dipoles.len());
        for d in dipoles.iter() {
            let li = ((d.position[0] - x_min) / cell_size).round() as usize;
            let mi = ((d.position[1] - y_min) / cell_size).round() as usize;
            let ni = ((d.position[2] - z_min) / cell_size).round() as usize;
            grid_index.push(li * (pad_ny * pad_nz) + mi * pad_nz + ni);
        }

        // Build and FFT the 9 kernel buffers.
        let mut kernel_fft: Vec<Vec<rustfft::num_complex::Complex<f64>>> =
            vec![vec![rustfft::num_complex::Complex::new(0.0f64, 0.0f64); fft_len]; 9];

        for component_idx in 0..9usize {
            let alpha = component_idx / 3; // output component
            let beta = component_idx % 3;  // input component

            let mut buf: Vec<rustfft::num_complex::Complex<f64>> =
                vec![rustfft::num_complex::Complex::new(0.0f64, 0.0f64); fft_len];

            // Fill the kernel buffer using circulant embedding.
            let nx_i = nx as isize;
            let ny_i = ny as isize;
            let nz_i = nz as isize;

            for dl in (-(nx_i - 1))..=(nx_i - 1) {
                for dm in (-(ny_i - 1))..=(ny_i - 1) {
                    for dn in (-(nz_i - 1))..=(nz_i - 1) {
                        // Skip the self-interaction (origin).
                        if dl == 0 && dm == 0 && dn == 0 {
                            continue;
                        }

                        let r1 = [dl as f64 * cell_size, dm as f64 * cell_size, dn as f64 * cell_size];
                        let r2 = [0.0f64; 3];
                        let g = dyadic_greens_tensor(&r1, &r2, k);

                        // Circulant embedding: use rem_euclid for correct negative-displacement mapping.
                        let cl = dl.rem_euclid(pad_nx as isize) as usize;
                        let cm = dm.rem_euclid(pad_ny as isize) as usize;
                        let cn = dn.rem_euclid(pad_nz as isize) as usize;
                        let flat = cl * (pad_ny * pad_nz) + cm * pad_nz + cn;

                        let g_val = g[alpha][beta];
                        buf[flat] = rustfft::num_complex::Complex::new(g_val.re, g_val.im);
                    }
                }
            }

            // 3D FFT of kernel buffer.
            fft3d_forward(&mut buf, pad_nx, pad_ny, pad_nz);

            kernel_fft[component_idx] = buf;
        }

        Some(FftMatvecPlan {
            nx,
            ny,
            nz,
            grid_index,
            kernel_fft,
            fft_len,
        })
    }

    /// Compute the matrix-vector product A·x using FFT convolution.
    ///
    /// A = α⁻¹ (diagonal) − G (off-diagonal), so:
    /// - Diagonal contribution: α⁻¹·x per dipole.
    /// - Off-diagonal contribution: subtract FFT-convolution result (G·x).
    pub fn matvec(&self, dipoles: &[Dipole], x: ArrayView1<Complex64>) -> Array1<Complex64> {
        let n = dipoles.len();
        let mut y = Array1::<Complex64>::zeros(3 * n);

        // ── Diagonal: α⁻¹ · x ──────────────────────────────────────────────
        for i in 0..n {
            let inv_alpha = invert_3x3_pub(&dipoles[i].polarisability);
            for a in 0..3 {
                for b in 0..3 {
                    y[3 * i + a] += inv_alpha[3 * a + b] * x[3 * i + b];
                }
            }
        }

        // ── Off-diagonal: −G·x via FFT convolution ─────────────────────────
        // For each source component β, scatter x into buffer, FFT, multiply
        // with each kernel row α, accumulate, IFFT, gather.
        let pad_ny = 2 * self.ny;
        let pad_nz = 2 * self.nz;

        // Accumulators for each output component α (in frequency domain).
        // We accumulate pointwise multiplications before doing inverse FFT.
        let mut acc: Vec<Vec<rustfft::num_complex::Complex<f64>>> = vec![
            vec![rustfft::num_complex::Complex::new(0.0f64, 0.0f64); self.fft_len]; 3
        ];

        for beta in 0..3usize {
            // Scatter: place x[3j+β] at grid position of dipole j.
            let mut buf: Vec<rustfft::num_complex::Complex<f64>> =
                vec![rustfft::num_complex::Complex::new(0.0f64, 0.0f64); self.fft_len];

            for j in 0..n {
                let xj = x[3 * j + beta];
                buf[self.grid_index[j]] = rustfft::num_complex::Complex::new(xj.re, xj.im);
            }

            // Forward FFT of the scattered input.
            fft3d_forward(&mut buf, 2 * self.nx, pad_ny, pad_nz);

            // Accumulate into each output component α.
            for alpha in 0..3usize {
                let k_fft = &self.kernel_fft[alpha * 3 + beta];
                for idx in 0..self.fft_len {
                    acc[alpha][idx] += k_fft[idx] * buf[idx];
                }
            }
        }

        // Inverse FFT for each output component and gather results.
        let norm = 1.0 / self.fft_len as f64;

        for alpha in 0..3usize {
            fft3d_inverse(&mut acc[alpha], 2 * self.nx, pad_ny, pad_nz);

            // Gather: subtract (A = α⁻¹ − G, so off-diagonal part is −G).
            for i in 0..n {
                let gathered = acc[alpha][self.grid_index[i]];
                let val = Complex64::new(gathered.re * norm, gathered.im * norm);
                y[3 * i + alpha] -= val;
            }
        }

        y
    }
}

/// Perform a 3D forward FFT on a flat row-major buffer of shape (nx, ny, nz).
///
/// Uses three passes of 1D FFTs:
/// - x-pass: ny*nz transforms of length nx
/// - y-pass: nx*nz transforms of length ny
/// - z-pass: nx*ny transforms of length nz
fn fft3d_forward(
    buf: &mut Vec<rustfft::num_complex::Complex<f64>>,
    nx: usize,
    ny: usize,
    nz: usize,
) {
    let mut planner = FftPlanner::<f64>::new();

    // z-pass: nx*ny contiguous transforms of length nz.
    {
        let fft = planner.plan_fft_forward(nz);
        for i in 0..nx {
            for j in 0..ny {
                let offset = i * ny * nz + j * nz;
                fft.process(&mut buf[offset..offset + nz]);
            }
        }
    }

    // y-pass: nx*nz transforms of length ny with stride nz.
    {
        let fft = planner.plan_fft_forward(ny);
        let mut tmp = vec![rustfft::num_complex::Complex::new(0.0f64, 0.0f64); ny];
        for i in 0..nx {
            for k in 0..nz {
                // Gather stride-nz elements.
                for j in 0..ny {
                    tmp[j] = buf[i * ny * nz + j * nz + k];
                }
                fft.process(&mut tmp);
                // Scatter back.
                for j in 0..ny {
                    buf[i * ny * nz + j * nz + k] = tmp[j];
                }
            }
        }
    }

    // x-pass: ny*nz transforms of length nx with stride ny*nz.
    {
        let fft = planner.plan_fft_forward(nx);
        let mut tmp = vec![rustfft::num_complex::Complex::new(0.0f64, 0.0f64); nx];
        for j in 0..ny {
            for k in 0..nz {
                // Gather stride-ny*nz elements.
                for i in 0..nx {
                    tmp[i] = buf[i * ny * nz + j * nz + k];
                }
                fft.process(&mut tmp);
                // Scatter back.
                for i in 0..nx {
                    buf[i * ny * nz + j * nz + k] = tmp[i];
                }
            }
        }
    }
}

/// Perform a 3D inverse FFT on a flat row-major buffer of shape (nx, ny, nz).
///
/// Uses three passes of 1D inverse FFTs. Note: rustfft does NOT normalise the
/// inverse FFT; the caller must divide by `nx * ny * nz` after the call.
fn fft3d_inverse(
    buf: &mut Vec<rustfft::num_complex::Complex<f64>>,
    nx: usize,
    ny: usize,
    nz: usize,
) {
    let mut planner = FftPlanner::<f64>::new();

    // z-pass.
    {
        let fft = planner.plan_fft_inverse(nz);
        for i in 0..nx {
            for j in 0..ny {
                let offset = i * ny * nz + j * nz;
                fft.process(&mut buf[offset..offset + nz]);
            }
        }
    }

    // y-pass.
    {
        let fft = planner.plan_fft_inverse(ny);
        let mut tmp = vec![rustfft::num_complex::Complex::new(0.0f64, 0.0f64); ny];
        for i in 0..nx {
            for k in 0..nz {
                for j in 0..ny {
                    tmp[j] = buf[i * ny * nz + j * nz + k];
                }
                fft.process(&mut tmp);
                for j in 0..ny {
                    buf[i * ny * nz + j * nz + k] = tmp[j];
                }
            }
        }
    }

    // x-pass.
    {
        let fft = planner.plan_fft_inverse(nx);
        let mut tmp = vec![rustfft::num_complex::Complex::new(0.0f64, 0.0f64); nx];
        for j in 0..ny {
            for k in 0..nz {
                for i in 0..nx {
                    tmp[i] = buf[i * ny * nz + j * nz + k];
                }
                fft.process(&mut tmp);
                for i in 0..nx {
                    buf[i * ny * nz + j * nz + k] = tmp[i];
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Dipole;
    use crate::solver::cda::assembly;

    fn make_dipole(pos: [f64; 3], alpha: Complex64) -> Dipole {
        let zero = Complex64::from(0.0);
        Dipole {
            position: pos,
            polarisability: [alpha, zero, zero, zero, alpha, zero, zero, zero, alpha],
        }
    }

    /// Build an N×N×N regular cubic grid of dipoles with the given spacing.
    fn make_grid(n: usize, spacing: f64, alpha: Complex64) -> Vec<Dipole> {
        let mut dipoles = Vec::with_capacity(n * n * n);
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    let pos = [
                        i as f64 * spacing,
                        j as f64 * spacing,
                        k as f64 * spacing,
                    ];
                    dipoles.push(make_dipole(pos, alpha));
                }
            }
        }
        dipoles
    }

    #[test]
    fn test_detect_regular_grid_known() {
        let alpha = Complex64::new(1000.0, 200.0);
        let dipoles = make_grid(4, 1.0, alpha);
        let result = detect_regular_grid(&dipoles, 1.0);
        assert_eq!(result, Some((4, 4, 4)));
    }

    #[test]
    fn test_detect_regular_grid_irregular() {
        let alpha = Complex64::new(1000.0, 200.0);
        let mut dipoles = make_grid(4, 1.0, alpha);
        // Displace one dipole by 0.1 nm.
        dipoles[7].position[0] += 0.1;
        let result = detect_regular_grid(&dipoles, 1.0);
        assert_eq!(result, None);
    }

    #[test]
    fn test_fft_matvec_agrees_with_on_the_fly() {
        let wavelength = 600.0_f64;
        let k = 2.0 * std::f64::consts::PI / wavelength;
        let spacing = 3.0_f64;
        let alpha = Complex64::new(1000.0, 200.0);

        let dipoles = make_grid(5, spacing, alpha);
        let n = dipoles.len(); // 125

        let plan = FftMatvecPlan::try_build(&dipoles, spacing, k)
            .expect("5×5×5 regular grid should build a plan");

        // Build a deterministic test vector.
        let x: Array1<Complex64> = Array1::from_iter((0..3 * n).map(|i| {
            Complex64::new(
                ((i * 7 + 1) % 17) as f64 / 17.0,
                ((i * 11 + 3) % 13) as f64 / 13.0,
            )
        }));

        let y_fft = plan.matvec(&dipoles, x.view());

        let y_otf = assembly::matvec_on_the_fly(
            &dipoles,
            &x,
            k,
            false,
            spacing,
            None,
            [0.0; 3],
            None,
        );

        // Compute relative error ||y_fft - y_otf|| / ||y_otf||.
        let diff_norm: f64 = y_fft
            .iter()
            .zip(y_otf.iter())
            .map(|(a, b)| (a - b).norm_sqr())
            .sum::<f64>()
            .sqrt();
        let ref_norm: f64 = y_otf.iter().map(|v| v.norm_sqr()).sum::<f64>().sqrt();

        let rel_err = diff_norm / ref_norm;
        assert!(
            rel_err < 1e-10,
            "FFT matvec relative error {:.2e} exceeds 1e-10 for 5×5×5 grid",
            rel_err
        );
    }

    #[test]
    fn test_fft_matvec_10x10x10() {
        let wavelength = 600.0_f64;
        let k = 2.0 * std::f64::consts::PI / wavelength;
        let spacing = 3.0_f64;
        let alpha = Complex64::new(1000.0, 200.0);

        let dipoles = make_grid(10, spacing, alpha);
        let n = dipoles.len(); // 1000

        let plan = FftMatvecPlan::try_build(&dipoles, spacing, k)
            .expect("10×10×10 regular grid should build a plan");

        let x: Array1<Complex64> = Array1::from_iter((0..3 * n).map(|i| {
            Complex64::new(
                ((i * 7 + 1) % 17) as f64 / 17.0,
                ((i * 11 + 3) % 13) as f64 / 13.0,
            )
        }));

        let y_fft = plan.matvec(&dipoles, x.view());

        let y_otf = assembly::matvec_on_the_fly(
            &dipoles,
            &x,
            k,
            false,
            spacing,
            None,
            [0.0; 3],
            None,
        );

        let diff_norm: f64 = y_fft
            .iter()
            .zip(y_otf.iter())
            .map(|(a, b)| (a - b).norm_sqr())
            .sum::<f64>()
            .sqrt();
        let ref_norm: f64 = y_otf.iter().map(|v| v.norm_sqr()).sum::<f64>().sqrt();

        let rel_err = diff_norm / ref_norm;
        assert!(
            rel_err < 1e-10,
            "FFT matvec relative error {:.2e} exceeds 1e-10 for 10×10×10 grid",
            rel_err
        );
    }
}
