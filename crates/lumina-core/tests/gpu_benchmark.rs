//! GPU vs CPU performance benchmark.
//!
//! Run with:
//!   cargo test -p lumina-core --features gpu --release -- gpu_benchmark --nocapture
//!
//! This is NOT a unit test — it prints timing comparisons. All assertions are
//! for correctness (results must agree within f32 tolerance), not performance.

use lumina_compute::{ComputeBackend, CpuBackend};
use ndarray::Array2;
use num_complex::Complex64;
use std::sync::Arc;
use std::time::Instant;

/// Build a random complex matrix of size `dim x dim` seeded deterministically.
fn random_matrix(dim: usize) -> Array2<Complex64> {
    let mut m = Array2::zeros((dim, dim));
    for i in 0..dim {
        for j in 0..dim {
            let re = ((i * 137 + j * 251 + 31) % 997) as f64 / 997.0 - 0.5;
            let im = ((i * 173 + j * 311 + 59) % 991) as f64 / 991.0 - 0.5;
            m[[i, j]] = Complex64::new(re, im);
        }
    }
    // Make diagonally dominant so GMRES converges
    for i in 0..dim {
        m[[i, i]] += Complex64::new(dim as f64, 0.0);
    }
    m
}

/// Build a random complex vector of length `dim`.
fn random_vector(dim: usize) -> ndarray::Array1<Complex64> {
    let mut v = ndarray::Array1::zeros(dim);
    for i in 0..dim {
        let re = ((i * 193 + 41) % 983) as f64 / 983.0 - 0.5;
        let im = ((i * 229 + 67) % 977) as f64 / 977.0 - 0.5;
        v[i] = Complex64::new(re, im);
    }
    v
}

/// Measure how long a closure takes and return (result, elapsed_ms).
fn timed<T>(f: impl FnOnce() -> T) -> (T, f64) {
    let start = Instant::now();
    let result = f();
    let elapsed = start.elapsed().as_secs_f64() * 1000.0;
    (result, elapsed)
}

#[test]
fn gpu_benchmark() {
    use lumina_core::solver::cda::iterative::solve_gmres;

    // Initialise backends
    let cpu: Arc<dyn ComputeBackend> = Arc::new(CpuBackend::new());

    #[cfg(feature = "gpu")]
    let gpu: Option<Arc<dyn ComputeBackend>> = {
        match lumina_compute::GpuBackend::new_blocking() {
            Ok(g) => {
                let info = g.device_info();
                println!("\n=== GPU Device: {} ===", info.name);
                Some(Arc::new(g))
            }
            Err(e) => {
                println!("\n=== GPU not available: {} — CPU-only results ===", e);
                None
            }
        }
    };
    #[cfg(not(feature = "gpu"))]
    let gpu: Option<Arc<dyn ComputeBackend>> = {
        println!("\n=== GPU feature not enabled — CPU-only results ===");
        None
    };

    // ── Matvec benchmark ──
    println!("\n--- Matvec (A*x) ---");
    println!(
        "{:<8} {:>12} {:>12} {:>10}",
        "Dim", "CPU (ms)", "GPU (ms)", "Speedup"
    );
    println!("{}", "-".repeat(46));

    for &dim in &[300, 1000, 3000] {
        let matrix = random_matrix(dim);
        let vector = random_vector(dim);

        // Warm up
        let _ = cpu.matvec(&matrix, &vector);
        if let Some(ref g) = gpu {
            let _ = g.matvec(&matrix, &vector);
        }

        // CPU timing
        let cpu_runs = if dim <= 1000 { 10 } else { 3 };
        let mut cpu_total = 0.0;
        let mut cpu_result = None;
        for _ in 0..cpu_runs {
            let (res, t) = timed(|| cpu.matvec(&matrix, &vector).unwrap());
            cpu_total += t;
            cpu_result = Some(res);
        }
        let cpu_ms = cpu_total / cpu_runs as f64;

        // GPU timing
        if let Some(ref g) = gpu {
            let gpu_runs = if dim <= 1000 { 10 } else { 3 };
            let mut gpu_total = 0.0;
            let mut gpu_result = None;
            for _ in 0..gpu_runs {
                let (res, t) = timed(|| g.matvec(&matrix, &vector).unwrap());
                gpu_total += t;
                gpu_result = Some(res);
            }
            let gpu_ms = gpu_total / gpu_runs as f64;
            let speedup = cpu_ms / gpu_ms;
            println!(
                "{:<8} {:>12.3} {:>12.3} {:>9.2}x",
                dim, cpu_ms, gpu_ms, speedup
            );

            // Verify correctness
            let cpu_v = cpu_result.unwrap();
            let gpu_v = gpu_result.unwrap();
            let max_diff: f64 = cpu_v
                .iter()
                .zip(gpu_v.iter())
                .map(|(a, b)| (a - b).norm())
                .fold(0.0_f64, f64::max);
            let norm: f64 = cpu_v.iter().map(|c| c.norm()).fold(0.0_f64, f64::max);
            let rel_err = max_diff / norm;
            println!(
                "         relative error: {:.2e} (f32 limit: ~1e-6)",
                rel_err
            );
            assert!(rel_err < 1e-4, "GPU/CPU matvec diverge at dim={dim}");
        } else {
            println!("{:<8} {:>12.3} {:>12} {:>10}", dim, cpu_ms, "N/A", "N/A");
        }
    }

    // ── GMRES solve benchmark ──
    println!("\n--- GMRES Solve ---");
    println!(
        "{:<8} {:>12} {:>12} {:>10}",
        "Dim", "CPU (ms)", "GPU (ms)", "Speedup"
    );
    println!("{}", "-".repeat(46));

    for &dim in &[300, 1000, 3000] {
        let matrix = random_matrix(dim);
        let rhs = random_vector(dim);

        // CPU GMRES (f64 precision → tight tolerance)
        let cpu_backend = Arc::clone(&cpu);
        let (cpu_result, cpu_ms) = timed(|| {
            let matvec_fn = |x: &ndarray::Array1<Complex64>| {
                cpu_backend
                    .matvec(&matrix, x)
                    .map_err(|e| lumina_core::solver::SolverError::LinAlgError(e.to_string()))
            };
            solve_gmres(&matvec_fn, &rhs, 1e-10, 1000).unwrap()
        });

        // GPU GMRES (f32 matvec → relaxed tolerance, ~1e-6 is the floor)
        if let Some(ref g) = gpu {
            let gpu_backend = Arc::clone(g);
            let (gpu_result, gpu_ms) = timed(|| {
                let matvec_fn = |x: &ndarray::Array1<Complex64>| {
                    gpu_backend
                        .matvec(&matrix, x)
                        .map_err(|e| {
                            lumina_core::solver::SolverError::LinAlgError(e.to_string())
                        })
                };
                solve_gmres(&matvec_fn, &rhs, 1e-6, 1000).unwrap()
            });
            let speedup = cpu_ms / gpu_ms;
            println!(
                "{:<8} {:>12.1} {:>12.1} {:>9.2}x",
                dim, cpu_ms, gpu_ms, speedup
            );

            // Verify solutions agree within f32 tolerance
            let max_diff: f64 = cpu_result
                .iter()
                .zip(gpu_result.iter())
                .map(|(a, b)| (a - b).norm())
                .fold(0.0_f64, f64::max);
            let norm: f64 = cpu_result
                .iter()
                .map(|c| c.norm())
                .fold(0.0_f64, f64::max);
            let rel_err = max_diff / norm;
            println!("         relative error: {:.2e}", rel_err);
            assert!(rel_err < 1e-3, "GPU/CPU GMRES diverge at dim={dim}");
        } else {
            println!("{:<8} {:>12.1} {:>12} {:>10}", dim, cpu_ms, "N/A", "N/A");
        }
    }

    println!();
}
