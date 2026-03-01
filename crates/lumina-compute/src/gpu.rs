//! GPU compute backend via wgpu with persistent session buffers.
//!
//! Implements [`ComputeBackend`] using wgpu compute shaders for the
//! performance-critical matrix-vector product used in GMRES iterations.
//!
//! # Persistent sessions
//!
//! For iterative (GMRES) solves, obtain a [`GpuSession`] via
//! [`ComputeBackend::create_session`]. The session pre-allocates all five GPU
//! buffers (matrix, input, output, staging, params) and builds the bind group
//! once. The interaction matrix is uploaded once via
//! [`MatvecSession::upload_matrix`]; each subsequent
//! [`MatvecSession::matvec`] call reuses all buffers without further
//! device allocations, eliminating the 7 per-iteration allocations of the
//! legacy path.
//!
//! # Precision
//!
//! WGSL has no native f64 support, so all GPU computation uses f32 (Complex32).
//! The f64 ↔ f32 conversion happens at the Rust boundary. This limits the
//! GMRES residual floor to ~1e-7, which is adequate for iterative solvers
//! where Arnoldi/Givens orthogonalisation remains in f64 on the CPU.
//!
//! # Thread safety
//!
//! `GpuBackend` and `GpuSession` are both `Send + Sync`. Multiple Rayon
//! threads (parallel wavelength sweeps) may each hold their own `GpuSession`
//! with independent buffers. Access to the shared wgpu device and command
//! queue is serialised through `Arc<Mutex<GpuState>>`, so GPU commands are
//! submitted one at a time — correct and race-free.

use std::borrow::Cow;
use std::sync::{Arc, Mutex};

use bytemuck::{Pod, Zeroable};
use ndarray::{Array1, Array2, ArrayView1};
use num_complex::Complex64;
use wgpu;

use crate::backend::{BackendType, ComputeBackend, ComputeError, DeviceInfo, MatvecSession};

// ─── Shader parameters ─────────────────────────────────────────────────────

/// Uniform buffer parameters passed to the WGSL shader.
#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
struct ShaderParams {
    dim: u32,
    _pad0: u32,
    _pad1: u32,
    _pad2: u32,
}

// ─── Internal state ─────────────────────────────────────────────────────────

/// Internal wgpu state shared between [`GpuBackend`] and all [`GpuSession`]s.
///
/// Protected by `Arc<Mutex<_>>` to serialise device/queue access across Rayon
/// threads running parallel wavelength sweeps.
struct GpuState {
    device: wgpu::Device,
    queue: wgpu::Queue,
    pipeline: wgpu::ComputePipeline,
    bind_group_layout: wgpu::BindGroupLayout,
}

// ─── Public types ───────────────────────────────────────────────────────────

/// GPU compute backend using wgpu.
///
/// Manages the wgpu device, queue, and compiled compute pipeline. For
/// iterative solves, use [`ComputeBackend::create_session`] to obtain a
/// [`GpuSession`] with pre-allocated persistent buffers.
pub struct GpuBackend {
    state: Arc<Mutex<GpuState>>,
    device_name: String,
}

/// A persistent GPU solve session with pre-allocated device buffers.
///
/// Created via [`ComputeBackend::create_session`]. Holds five GPU buffers
/// (matrix, input, output, staging, params) and a pre-built bind group.
/// Upload the matrix once with [`MatvecSession::upload_matrix`], then call
/// [`MatvecSession::matvec`] for each GMRES iteration without further
/// allocations.
///
/// Dropping the session frees all GPU buffers.
pub struct GpuSession {
    matrix_buf: wgpu::Buffer,
    input_buf: wgpu::Buffer,
    output_buf: wgpu::Buffer,
    staging_buf: wgpu::Buffer,
    // params_buf is owned here to keep the GPU buffer alive for the bind group.
    #[allow(dead_code)]
    params_buf: wgpu::Buffer,
    bind_group: wgpu::BindGroup,
    dim: usize,
    state: Arc<Mutex<GpuState>>,
}

// ─── GpuBackend implementation ──────────────────────────────────────────────

impl GpuBackend {
    /// Create a new GPU backend asynchronously.
    ///
    /// Requests a high-performance adapter, compiles the WGSL shader, and
    /// builds the compute pipeline. Falls back to any available adapter if
    /// no discrete GPU is found.
    pub async fn new_async() -> Result<Self, ComputeError> {
        let instance = wgpu::Instance::new(&wgpu::InstanceDescriptor {
            backends: wgpu::Backends::all(),
            ..Default::default()
        });

        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                force_fallback_adapter: false,
                compatible_surface: None,
            })
            .await
            .ok_or_else(|| ComputeError::Unavailable("No GPU adapter found".into()))?;

        let device_name = adapter.get_info().name.clone();

        let (device, queue) = adapter
            .request_device(
                &wgpu::DeviceDescriptor {
                    label: Some("lumina-compute"),
                    required_features: wgpu::Features::empty(),
                    required_limits: wgpu::Limits::default(),
                    memory_hints: wgpu::MemoryHints::Performance,
                },
                None,
            )
            .await
            .map_err(|e| ComputeError::DeviceError(format!("Failed to create device: {}", e)))?;

        let shader_source = include_str!("shaders/matvec.wgsl");
        let shader_module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("matvec_shader"),
            source: wgpu::ShaderSource::Wgsl(Cow::Borrowed(shader_source)),
        });

        let bind_group_layout =
            device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("matvec_bind_group_layout"),
                entries: &[
                    // binding 0: matrix (storage, read-only)
                    wgpu::BindGroupLayoutEntry {
                        binding: 0,
                        visibility: wgpu::ShaderStages::COMPUTE,
                        ty: wgpu::BindingType::Buffer {
                            ty: wgpu::BufferBindingType::Storage { read_only: true },
                            has_dynamic_offset: false,
                            min_binding_size: None,
                        },
                        count: None,
                    },
                    // binding 1: input vector (storage, read-only)
                    wgpu::BindGroupLayoutEntry {
                        binding: 1,
                        visibility: wgpu::ShaderStages::COMPUTE,
                        ty: wgpu::BindingType::Buffer {
                            ty: wgpu::BufferBindingType::Storage { read_only: true },
                            has_dynamic_offset: false,
                            min_binding_size: None,
                        },
                        count: None,
                    },
                    // binding 2: output vector (storage, read-write)
                    wgpu::BindGroupLayoutEntry {
                        binding: 2,
                        visibility: wgpu::ShaderStages::COMPUTE,
                        ty: wgpu::BindingType::Buffer {
                            ty: wgpu::BufferBindingType::Storage { read_only: false },
                            has_dynamic_offset: false,
                            min_binding_size: None,
                        },
                        count: None,
                    },
                    // binding 3: params uniform
                    wgpu::BindGroupLayoutEntry {
                        binding: 3,
                        visibility: wgpu::ShaderStages::COMPUTE,
                        ty: wgpu::BindingType::Buffer {
                            ty: wgpu::BufferBindingType::Uniform,
                            has_dynamic_offset: false,
                            min_binding_size: None,
                        },
                        count: None,
                    },
                ],
            });

        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("matvec_pipeline_layout"),
            bind_group_layouts: &[&bind_group_layout],
            push_constant_ranges: &[],
        });

        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("matvec_pipeline"),
            layout: Some(&pipeline_layout),
            module: &shader_module,
            entry_point: Some("main"),
            compilation_options: wgpu::PipelineCompilationOptions::default(),
            cache: None,
        });

        Ok(Self {
            state: Arc::new(Mutex::new(GpuState {
                device,
                queue,
                pipeline,
                bind_group_layout,
            })),
            device_name,
        })
    }

    /// Create a new GPU backend, blocking the current thread.
    ///
    /// Convenience wrapper around [`new_async`](Self::new_async) using
    /// `pollster` to block on the async initialisation.
    pub fn new_blocking() -> Result<Self, ComputeError> {
        pollster::block_on(Self::new_async())
    }

    /// Allocate a persistent solve session for `n` unknowns (typically `3 * N_dipoles`).
    ///
    /// Acquires the state lock briefly to allocate five GPU buffers and build
    /// the bind group, then releases the lock. All subsequent GPU operations
    /// on the session take the lock only for the duration of each matvec.
    pub fn new_session(&self, n: usize) -> Result<GpuSession, ComputeError> {
        let state = self.state.lock().unwrap();

        let vec_size = (n * std::mem::size_of::<[f32; 2]>()) as u64;
        let mat_size = (n * n * std::mem::size_of::<[f32; 2]>()) as u64;

        let matrix_buf = state.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("session_matrix"),
            size: mat_size,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let input_buf = state.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("session_input"),
            size: vec_size,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let output_buf = state.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("session_output"),
            size: vec_size,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        let staging_buf = state.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("session_staging"),
            size: vec_size,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        // params_buf holds the fixed dimension; written once at creation.
        let params = ShaderParams { dim: n as u32, _pad0: 0, _pad1: 0, _pad2: 0 };
        let params_buf = state.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("session_params"),
            size: std::mem::size_of::<ShaderParams>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        state.queue.write_buffer(&params_buf, 0, bytemuck::bytes_of(&params));

        // Pre-build the bind group — reused for every matvec in this session.
        let bind_group = state.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("session_bind_group"),
            layout: &state.bind_group_layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: matrix_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: input_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: output_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: params_buf.as_entire_binding(),
                },
            ],
        });

        drop(state);

        Ok(GpuSession {
            matrix_buf,
            input_buf,
            output_buf,
            staging_buf,
            params_buf,
            bind_group,
            dim: n,
            state: Arc::clone(&self.state),
        })
    }
}

// ─── MatvecSession for GpuSession ───────────────────────────────────────────

impl MatvecSession for GpuSession {
    /// Upload the interaction matrix to the GPU matrix buffer.
    ///
    /// Converts the matrix from `Complex64` to `f32` pairs and writes it to
    /// the pre-allocated matrix buffer via `queue.write_buffer`. Called once
    /// before the GMRES iteration loop.
    fn upload_matrix(&self, matrix: &Array2<Complex64>) -> Result<(), ComputeError> {
        if matrix.nrows() != self.dim || matrix.ncols() != self.dim {
            return Err(ComputeError::DeviceError(format!(
                "Session expects {}×{} matrix, got {}×{}",
                self.dim, self.dim, matrix.nrows(), matrix.ncols()
            )));
        }
        let data = matrix_to_f32(matrix);
        let state = self.state.lock().unwrap();
        state.queue.write_buffer(&self.matrix_buf, 0, bytemuck::cast_slice(&data));
        Ok(())
    }

    /// Perform $\mathbf{y} = \mathbf{A}\mathbf{x}$ using the pre-uploaded matrix.
    ///
    /// Per-call operations: write `x` to `input_buf`, dispatch the compute
    /// shader, copy `output_buf` to `staging_buf`, map and read the result.
    /// All buffer allocations and the bind group are reused from the session.
    fn matvec(&self, x: ArrayView1<Complex64>) -> Result<Array1<Complex64>, ComputeError> {
        if x.len() != self.dim {
            return Err(ComputeError::DeviceError(format!(
                "Session dim={}, got vector len={}",
                self.dim,
                x.len()
            )));
        }

        let input_data = vector_to_f32(x);
        let vec_size = (self.dim * std::mem::size_of::<[f32; 2]>()) as u64;

        // Lock for the entire GPU operation — serialises concurrent wavelength threads.
        let state = self.state.lock().unwrap();

        // Write input vector (reuses pre-allocated buffer).
        state.queue.write_buffer(&self.input_buf, 0, bytemuck::cast_slice(&input_data));

        // Encode compute pass using the pre-built bind group.
        let mut encoder = state
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("session_matvec_encoder"),
            });
        {
            let mut cpass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("session_matvec_pass"),
                timestamp_writes: None,
            });
            cpass.set_pipeline(&state.pipeline);
            cpass.set_bind_group(0, &self.bind_group, &[]);
            let num_workgroups = (self.dim as u32).div_ceil(256);
            cpass.dispatch_workgroups(num_workgroups, 1, 1);
        }
        encoder.copy_buffer_to_buffer(&self.output_buf, 0, &self.staging_buf, 0, vec_size);
        state.queue.submit(std::iter::once(encoder.finish()));

        // Map the staging buffer and wait for completion.
        let staging_slice = self.staging_buf.slice(..);
        let (sender, receiver) = std::sync::mpsc::channel();
        staging_slice.map_async(wgpu::MapMode::Read, move |result| {
            sender.send(result).unwrap();
        });
        state.device.poll(wgpu::Maintain::Wait);

        // Release the lock before blocking on channel recv (poll is complete).
        drop(state);

        receiver
            .recv()
            .map_err(|e| ComputeError::DeviceError(format!("Map recv failed: {}", e)))?
            .map_err(|e| ComputeError::DeviceError(format!("Buffer map failed: {}", e)))?;

        let data = staging_slice.get_mapped_range();
        let result_f32: &[[f32; 2]] = bytemuck::cast_slice(&data);
        let result = f32_to_vector(result_f32, self.dim);
        drop(data);
        self.staging_buf.unmap();

        Ok(result)
    }
}

// ─── ComputeBackend for GpuBackend ──────────────────────────────────────────

impl ComputeBackend for GpuBackend {
    fn device_info(&self) -> DeviceInfo {
        DeviceInfo {
            name: format!("GPU ({})", self.device_name),
            backend_type: BackendType::Gpu,
            memory_bytes: None,
            compute_units: None,
        }
    }

    /// Create a persistent GPU session for `n` unknowns.
    fn create_session(&self, n: usize) -> Option<Box<dyn MatvecSession>> {
        self.new_session(n).ok().map(|s| Box::new(s) as Box<dyn MatvecSession>)
    }

    fn parallel_matrix_fill(
        &self,
        rows: usize,
        cols: usize,
        fill_fn: &(dyn Fn(usize, usize) -> Complex64 + Send + Sync),
    ) -> Result<Array2<Complex64>, ComputeError> {
        // Matrix assembly stays on CPU (Rayon). The GPU backend only accelerates
        // matvec — assembly is branching-heavy and GPU-unfriendly.
        use rayon::prelude::*;

        let data: Vec<Complex64> = (0..rows * cols)
            .into_par_iter()
            .map(|idx| {
                let i = idx / cols;
                let j = idx % cols;
                fill_fn(i, j)
            })
            .collect();

        Array2::from_shape_vec((rows, cols), data)
            .map_err(|e| ComputeError::DeviceError(e.to_string()))
    }

    /// Legacy one-off matvec. For GMRES, prefer [`create_session`](ComputeBackend::create_session).
    ///
    /// Creates a temporary session, uploads the matrix, runs the matvec, and
    /// drops the session. Used by tests and the direct-solve compatibility path.
    fn matvec(
        &self,
        matrix: &Array2<Complex64>,
        vector: &Array1<Complex64>,
    ) -> Result<Array1<Complex64>, ComputeError> {
        let dim = matrix.nrows();
        if dim != matrix.ncols() || dim != vector.len() {
            return Err(ComputeError::DeviceError(format!(
                "Dimension mismatch: matrix {}×{}, vector {}",
                matrix.nrows(),
                matrix.ncols(),
                vector.len()
            )));
        }
        let session = self.new_session(dim)?;
        session.upload_matrix(matrix)?;
        session.matvec(vector.view())
    }
}

// ─── Conversion helpers ─────────────────────────────────────────────────────

/// Convert a Complex64 matrix to flat `[f32; 2]` pairs in row-major order.
fn matrix_to_f32(matrix: &Array2<Complex64>) -> Vec<[f32; 2]> {
    let (rows, cols) = matrix.dim();
    let mut out = Vec::with_capacity(rows * cols);
    for i in 0..rows {
        for j in 0..cols {
            let c = matrix[[i, j]];
            out.push([c.re as f32, c.im as f32]);
        }
    }
    out
}

/// Convert a Complex64 vector view to `[f32; 2]` pairs.
fn vector_to_f32(v: ArrayView1<Complex64>) -> Vec<[f32; 2]> {
    v.iter().map(|c| [c.re as f32, c.im as f32]).collect()
}

/// Convert `[f32; 2]` pairs back to a Complex64 ndarray vector.
fn f32_to_vector(data: &[[f32; 2]], dim: usize) -> Array1<Complex64> {
    Array1::from_vec(
        data[..dim]
            .iter()
            .map(|pair| Complex64::new(pair[0] as f64, pair[1] as f64))
            .collect(),
    )
}

// ─── Tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rayon::prelude::*;

    fn try_gpu() -> Option<GpuBackend> {
        GpuBackend::new_blocking().ok()
    }

    #[test]
    fn test_gpu_backend_creation() {
        match GpuBackend::new_blocking() {
            Ok(gpu) => {
                let info = gpu.device_info();
                assert_eq!(info.backend_type, BackendType::Gpu);
                println!("GPU backend created: {}", info.name);
            }
            Err(e) => {
                println!("GPU not available (expected in CI): {}", e);
            }
        }
    }

    #[test]
    fn test_gpu_matvec_identity() {
        let gpu = match try_gpu() {
            Some(g) => g,
            None => { println!("Skipping: no GPU"); return; }
        };

        let dim = 6;
        let mut matrix = Array2::<Complex64>::zeros((dim, dim));
        for i in 0..dim { matrix[[i, i]] = Complex64::new(1.0, 0.0); }
        let vector = Array1::from_vec(
            (0..dim).map(|i| Complex64::new(i as f64 + 1.0, 0.5 * i as f64)).collect(),
        );

        let result = gpu.matvec(&matrix, &vector).unwrap();
        for i in 0..dim {
            let err = (result[i] - vector[i]).norm();
            assert!(err < 1e-4, "Identity matvec err[{}] = {:.2e}", i, err);
        }
    }

    #[test]
    fn test_gpu_matvec_complex_matrix() {
        let gpu = match try_gpu() {
            Some(g) => g,
            None => { println!("Skipping: no GPU"); return; }
        };

        let dim = 3;
        let matrix = Array2::from_shape_vec(
            (dim, dim),
            vec![
                Complex64::new(2.0, 1.0),  Complex64::new(-1.0, 0.5), Complex64::new(0.0, 0.0),
                Complex64::new(0.3, -0.2), Complex64::new(3.0, 0.0),  Complex64::new(-0.5, 0.1),
                Complex64::new(0.0, 0.0),  Complex64::new(0.7, -0.3), Complex64::new(4.0, -1.0),
            ],
        ).unwrap();
        let vector = Array1::from_vec(vec![
            Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 1.0),
            Complex64::new(-1.0, 0.5),
        ]);
        let expected = matrix.dot(&vector);
        let result = gpu.matvec(&matrix, &vector).unwrap();
        for i in 0..dim {
            let err = (result[i] - expected[i]).norm();
            assert!(err < 1e-4, "Matvec mismatch [{}] err={:.2e}", i, err);
        }
    }

    #[test]
    fn test_gpu_matvec_larger_system() {
        let gpu = match try_gpu() {
            Some(g) => g,
            None => { println!("Skipping: no GPU"); return; }
        };

        let dim = 300;
        let mut matrix = Array2::<Complex64>::zeros((dim, dim));
        for i in 0..dim {
            matrix[[i, i]] = Complex64::new(4.0, 0.5);
            if i + 1 < dim { matrix[[i, i + 1]] = Complex64::new(-1.0, 0.1); }
            if i > 0 { matrix[[i, i - 1]] = Complex64::new(-1.0, -0.1); }
        }
        let vector = Array1::from_vec(
            (0..dim).map(|i| Complex64::new((i as f64).sin(), (i as f64).cos())).collect(),
        );
        let expected = matrix.dot(&vector);
        let result = gpu.matvec(&matrix, &vector).unwrap();
        let expected_norm: f64 = expected.iter().map(|c| c.norm_sqr()).sum::<f64>().sqrt();
        let max_err = (0..dim).map(|i| (result[i] - expected[i]).norm()).fold(0.0_f64, f64::max);
        let rel_err = max_err / expected_norm;
        assert!(rel_err < 1e-4, "300×300 matvec rel_err={:.2e}", rel_err);
    }

    // ─── New session tests ───────────────────────────────────────────────

    /// Session matvec must match CPU reference to within f32 precision floor.
    #[test]
    fn test_persistent_session_matches_direct() {
        let gpu = match try_gpu() {
            Some(g) => g,
            None => { println!("Skipping: no GPU"); return; }
        };

        let dim = 6;
        let matrix = Array2::from_shape_fn((dim, dim), |(i, j)| {
            Complex64::new((i as f64 + 1.0) * 0.5 + j as f64 * 0.1,
                           (i as f64 - j as f64) * 0.2)
        });
        let vector = Array1::from_vec(
            (0..dim).map(|i| Complex64::new(i as f64 * 0.3 + 1.0, i as f64 * 0.1)).collect(),
        );
        let expected = matrix.dot(&vector);

        let session = gpu.new_session(dim).expect("session creation failed");
        session.upload_matrix(&matrix).expect("upload failed");
        let result = session.matvec(vector.view()).expect("matvec failed");

        for i in 0..dim {
            let err = (result[i] - expected[i]).norm();
            assert!(
                err < 1e-5,
                "Session matvec mismatch [{}]: GPU={:?}, CPU={:?}, err={:.2e}",
                i, result[i], expected[i], err
            );
        }
    }

    /// Multiple matvecs on the same session must all agree with the CPU reference.
    #[test]
    fn test_session_matvec_repeated() {
        let gpu = match try_gpu() {
            Some(g) => g,
            None => { println!("Skipping: no GPU"); return; }
        };

        let dim = 9;
        let matrix = Array2::from_shape_fn((dim, dim), |(i, j)| {
            if i == j { Complex64::new(3.0, 0.5) }
            else { Complex64::new(0.1 * (i as f64 - j as f64).signum(), 0.05) }
        });
        let session = gpu.new_session(dim).expect("session creation failed");
        session.upload_matrix(&matrix).expect("upload failed");

        for trial in 0..5 {
            let vector = Array1::from_vec(
                (0..dim).map(|i| Complex64::new((i + trial) as f64 * 0.4, (trial as f64 - i as f64) * 0.1)).collect(),
            );
            let expected = matrix.dot(&vector);
            let result = session.matvec(vector.view()).expect("matvec failed");
            for i in 0..dim {
                let err = (result[i] - expected[i]).norm();
                assert!(err < 1e-4, "Trial {trial} mismatch [{}] err={:.2e}", i, err);
            }
        }
    }

    /// Parallel wavelength threads each get their own session — no panics, results match CPU.
    #[test]
    fn test_parallel_wavelengths_gpu() {
        use std::sync::Arc;

        let gpu = match try_gpu() {
            Some(g) => g,
            None => { println!("Skipping: no GPU"); return; }
        };
        let gpu = Arc::new(gpu);

        let dim = 12;
        let matrix = Array2::from_shape_fn((dim, dim), |(i, j)| {
            if i == j { Complex64::new(5.0, 1.0) }
            else { Complex64::new(-0.1, 0.05 * (i as f64 - j as f64)) }
        });
        let expected_base = {
            let v = Array1::from_vec((0..dim).map(|i| Complex64::new(i as f64, 0.0)).collect());
            matrix.dot(&v)
        };

        let results: Vec<_> = (0_usize..4).into_par_iter().map(|thread_id| {
            let gpu = Arc::clone(&gpu);
            let session = gpu.new_session(dim).expect("session creation");
            session.upload_matrix(&matrix).expect("upload");
            let v = Array1::from_vec(
                (0..dim).map(|i| Complex64::new(i as f64, thread_id as f64 * 0.1)).collect(),
            );
            let res = session.matvec(v.view()).expect("matvec");
            let expected = matrix.dot(&v);
            (res, expected)
        }).collect();

        for (i, (result, expected)) in results.iter().enumerate() {
            let max_err = (0..dim).map(|j| (result[j] - expected[j]).norm()).fold(0.0_f64, f64::max);
            assert!(max_err < 1e-4, "Thread {i} max_err={:.2e}", max_err);
        }
        let _ = expected_base;
    }
}
