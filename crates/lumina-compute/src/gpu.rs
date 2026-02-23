//! GPU compute backend via wgpu.
//!
//! Implements [`ComputeBackend`] using wgpu compute shaders for the
//! performance-critical matrix-vector product used in GMRES iterations.
//!
//! # Precision
//!
//! WGSL has no native f64 support, so all GPU computation uses f32
//! (Complex32). The f64 ↔ f32 conversion happens at the Rust boundary.
//! This is adequate for iterative solvers where the Arnoldi/Givens
//! orthogonalisation remains in f64 on the CPU.
//!
//! # Matrix caching
//!
//! The interaction matrix is uploaded to the GPU once and reused for all
//! GMRES iterations (~300 matvecs per solve). The matrix buffer is cached
//! and only re-uploaded when the dimension changes.

use std::borrow::Cow;
use std::sync::Mutex;

use bytemuck::{Pod, Zeroable};
use ndarray::{Array1, Array2};
use num_complex::Complex64;
use wgpu;

use crate::backend::{BackendType, ComputeBackend, ComputeError, DeviceInfo};

/// Uniform buffer parameters passed to the WGSL shader.
#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
struct ShaderParams {
    dim: u32,
    _pad0: u32,
    _pad1: u32,
    _pad2: u32,
}

/// Cached GPU matrix state for reuse across GMRES iterations.
struct CachedMatrix {
    /// The GPU buffer holding the matrix data.
    buffer: wgpu::Buffer,
    /// Dimension of the matrix (number of rows = number of columns).
    dim: usize,
}

/// GPU compute backend using wgpu.
///
/// The backend holds the wgpu device, queue, and compiled compute pipeline.
/// It caches the interaction matrix on the GPU for efficient reuse during
/// iterative solves.
pub struct GpuBackend {
    device: wgpu::Device,
    queue: wgpu::Queue,
    pipeline: wgpu::ComputePipeline,
    bind_group_layout: wgpu::BindGroupLayout,
    device_name: String,
    /// Cached matrix buffer from the last `matvec` call.
    cached_matrix: Mutex<Option<CachedMatrix>>,
}

impl GpuBackend {
    /// Create a new GPU backend asynchronously.
    ///
    /// Requests a discrete GPU adapter if available, otherwise falls back
    /// to any available adapter.
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

        // Compile the WGSL compute shader.
        let shader_source = include_str!("shaders/matvec.wgsl");
        let shader_module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("matvec_shader"),
            source: wgpu::ShaderSource::Wgsl(Cow::Borrowed(shader_source)),
        });

        // Create explicit bind group layout matching the shader bindings.
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
            device,
            queue,
            pipeline,
            bind_group_layout,
            device_name,
            cached_matrix: Mutex::new(None),
        })
    }

    /// Create a new GPU backend, blocking the current thread.
    ///
    /// This is a convenience wrapper around [`new_async`](Self::new_async)
    /// using `pollster` to block on the async initialisation.
    pub fn new_blocking() -> Result<Self, ComputeError> {
        pollster::block_on(Self::new_async())
    }

    /// Upload a matrix to the GPU, reusing the cached buffer if the
    /// dimension matches.
    fn ensure_matrix_uploaded(&self, matrix: &Array2<Complex64>) {
        let dim = matrix.nrows();
        let mut cache = self.cached_matrix.lock().unwrap();

        // Reuse existing buffer if dimension matches.
        if let Some(ref cached) = *cache {
            if cached.dim == dim {
                // Same dimension — overwrite the buffer contents.
                let data = matrix_to_f32(matrix);
                self.queue
                    .write_buffer(&cached.buffer, 0, bytemuck::cast_slice(&data));
                return;
            }
        }

        // Dimension changed or first upload — create a new buffer.
        let data = matrix_to_f32(matrix);
        let buffer = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("matrix_buffer"),
            size: (data.len() * std::mem::size_of::<[f32; 2]>()) as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        self.queue
            .write_buffer(&buffer, 0, bytemuck::cast_slice(&data));

        *cache = Some(CachedMatrix { buffer, dim });
    }

    /// Perform the GPU matvec and read back the result.
    fn gpu_matvec_inner(
        &self,
        vector: &Array1<Complex64>,
        dim: usize,
    ) -> Result<Array1<Complex64>, ComputeError> {
        let cache = self.cached_matrix.lock().unwrap();
        let cached = cache.as_ref().ok_or_else(|| {
            ComputeError::DeviceError("Matrix not uploaded to GPU".into())
        })?;

        // Convert input vector to f32.
        let input_data = vector_to_f32(vector);

        let vec_byte_size = (dim * std::mem::size_of::<[f32; 2]>()) as u64;

        // Create input vector buffer.
        let input_buffer = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("input_vec"),
            size: vec_byte_size,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        self.queue
            .write_buffer(&input_buffer, 0, bytemuck::cast_slice(&input_data));

        // Create output vector buffer.
        let output_buffer = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("output_vec"),
            size: vec_byte_size,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });

        // Create staging buffer for readback.
        let staging_buffer = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("staging"),
            size: vec_byte_size,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        // Create params uniform buffer.
        let params = ShaderParams {
            dim: dim as u32,
            _pad0: 0,
            _pad1: 0,
            _pad2: 0,
        };
        let params_buffer = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("params"),
            size: std::mem::size_of::<ShaderParams>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        self.queue
            .write_buffer(&params_buffer, 0, bytemuck::bytes_of(&params));

        // Create bind group.
        let bind_group = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("matvec_bind_group"),
            layout: &self.bind_group_layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: cached.buffer.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: input_buffer.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: output_buffer.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: params_buffer.as_entire_binding(),
                },
            ],
        });

        // Encode compute pass.
        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("matvec_encoder"),
            });

        {
            let mut cpass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("matvec_pass"),
                timestamp_writes: None,
            });
            cpass.set_pipeline(&self.pipeline);
            cpass.set_bind_group(0, &bind_group, &[]);

            // Dispatch: one thread per row, workgroup size = 256.
            let num_workgroups = (dim as u32).div_ceil(256);
            cpass.dispatch_workgroups(num_workgroups, 1, 1);
        }

        // Copy output to staging buffer for readback.
        encoder.copy_buffer_to_buffer(&output_buffer, 0, &staging_buffer, 0, vec_byte_size);

        // Submit and wait.
        self.queue.submit(std::iter::once(encoder.finish()));

        // Drop the lock before blocking on device poll.
        drop(cache);

        // Map the staging buffer for reading.
        let staging_slice = staging_buffer.slice(..);

        let (sender, receiver) = std::sync::mpsc::channel();
        staging_slice.map_async(wgpu::MapMode::Read, move |result| {
            sender.send(result).unwrap();
        });

        self.device.poll(wgpu::Maintain::Wait);

        receiver
            .recv()
            .map_err(|e| ComputeError::DeviceError(format!("Map recv failed: {}", e)))?
            .map_err(|e| ComputeError::DeviceError(format!("Buffer map failed: {}", e)))?;

        // Read the mapped data and convert back to f64.
        let data = staging_slice.get_mapped_range();
        let result_f32: &[[f32; 2]] = bytemuck::cast_slice(&data);
        let result = f32_to_vector(result_f32, dim);

        drop(data);
        staging_buffer.unmap();

        Ok(result)
    }
}

impl ComputeBackend for GpuBackend {
    fn device_info(&self) -> DeviceInfo {
        DeviceInfo {
            name: format!("GPU ({})", self.device_name),
            backend_type: BackendType::Gpu,
            memory_bytes: None,
            compute_units: None,
        }
    }

    fn parallel_matrix_fill(
        &self,
        rows: usize,
        cols: usize,
        fill_fn: &(dyn Fn(usize, usize) -> Complex64 + Send + Sync),
    ) -> Result<Array2<Complex64>, ComputeError> {
        // Matrix assembly stays on CPU (Rayon). The GPU backend only
        // accelerates matvec — assembly is I/O bound, not compute bound,
        // and the element function uses branching that GPUs handle poorly.
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

    fn matvec(
        &self,
        matrix: &Array2<Complex64>,
        vector: &Array1<Complex64>,
    ) -> Result<Array1<Complex64>, ComputeError> {
        let dim = matrix.nrows();

        if dim != matrix.ncols() || dim != vector.len() {
            return Err(ComputeError::DeviceError(format!(
                "Dimension mismatch: matrix {}x{}, vector {}",
                matrix.nrows(),
                matrix.ncols(),
                vector.len()
            )));
        }

        // Upload matrix (reuses cache if dimension matches).
        self.ensure_matrix_uploaded(matrix);

        // Run the GPU matvec.
        self.gpu_matvec_inner(vector, dim)
    }
}

// ─── Conversion helpers ────────────────────────────────────────────────

/// Convert an ndarray Complex64 matrix to a flat Vec of [f32; 2] pairs
/// (interleaved real/imag) in row-major order.
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

/// Convert an ndarray Complex64 vector to a Vec of [f32; 2] pairs.
fn vector_to_f32(vector: &Array1<Complex64>) -> Vec<[f32; 2]> {
    vector.iter().map(|c| [c.re as f32, c.im as f32]).collect()
}

/// Convert a slice of [f32; 2] pairs back to an ndarray Complex64 vector.
fn f32_to_vector(data: &[[f32; 2]], dim: usize) -> Array1<Complex64> {
    Array1::from_vec(
        data[..dim]
            .iter()
            .map(|pair| Complex64::new(pair[0] as f64, pair[1] as f64))
            .collect(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Try to create a GpuBackend. Skip the test if no GPU is available.
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
            None => {
                println!("Skipping: no GPU available");
                return;
            }
        };

        // Identity matrix: A*x = x
        let dim = 6;
        let mut matrix = Array2::<Complex64>::zeros((dim, dim));
        for i in 0..dim {
            matrix[[i, i]] = Complex64::new(1.0, 0.0);
        }
        let vector = Array1::from_vec(
            (0..dim)
                .map(|i| Complex64::new(i as f64 + 1.0, 0.5 * i as f64))
                .collect(),
        );

        let result = gpu.matvec(&matrix, &vector).unwrap();

        for i in 0..dim {
            let err = (result[i] - vector[i]).norm();
            assert!(
                err < 1e-4,
                "Identity matvec: result[{}] = {:?}, expected {:?}, err = {:.2e}",
                i, result[i], vector[i], err
            );
        }
    }

    #[test]
    fn test_gpu_matvec_complex_matrix() {
        let gpu = match try_gpu() {
            Some(g) => g,
            None => {
                println!("Skipping: no GPU available");
                return;
            }
        };

        // Small dense complex matrix
        let dim = 3;
        let matrix = Array2::from_shape_vec(
            (dim, dim),
            vec![
                Complex64::new(2.0, 1.0),
                Complex64::new(-1.0, 0.5),
                Complex64::new(0.0, 0.0),
                Complex64::new(0.3, -0.2),
                Complex64::new(3.0, 0.0),
                Complex64::new(-0.5, 0.1),
                Complex64::new(0.0, 0.0),
                Complex64::new(0.7, -0.3),
                Complex64::new(4.0, -1.0),
            ],
        )
        .unwrap();

        let vector = Array1::from_vec(vec![
            Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 1.0),
            Complex64::new(-1.0, 0.5),
        ]);

        // CPU reference
        let expected = matrix.dot(&vector);

        let result = gpu.matvec(&matrix, &vector).unwrap();

        for i in 0..dim {
            let err = (result[i] - expected[i]).norm();
            assert!(
                err < 1e-4,
                "Matvec mismatch at [{}]: GPU={:?}, CPU={:?}, err={:.2e}",
                i, result[i], expected[i], err
            );
        }
    }

    #[test]
    fn test_gpu_matvec_larger_system() {
        let gpu = match try_gpu() {
            Some(g) => g,
            None => {
                println!("Skipping: no GPU available");
                return;
            }
        };

        // 300x300 system (100 dipoles × 3 components)
        let dim = 300;
        let mut matrix = Array2::<Complex64>::zeros((dim, dim));
        for i in 0..dim {
            matrix[[i, i]] = Complex64::new(4.0, 0.5);
            if i + 1 < dim {
                matrix[[i, i + 1]] = Complex64::new(-1.0, 0.1);
            }
            if i > 0 {
                matrix[[i, i - 1]] = Complex64::new(-1.0, -0.1);
            }
        }

        let vector = Array1::from_vec(
            (0..dim)
                .map(|i| Complex64::new((i as f64).sin(), (i as f64).cos()))
                .collect(),
        );

        let expected = matrix.dot(&vector);
        let result = gpu.matvec(&matrix, &vector).unwrap();

        let mut max_err: f64 = 0.0;
        for i in 0..dim {
            let err = (result[i] - expected[i]).norm();
            max_err = max_err.max(err);
        }

        // f32 precision: relative error should be < 1e-3 for well-conditioned systems
        let expected_norm: f64 = expected.iter().map(|c| c.norm_sqr()).sum::<f64>().sqrt();
        let rel_err = max_err / expected_norm;
        assert!(
            rel_err < 1e-4,
            "300x300 matvec relative error {:.2e} too large (max abs err = {:.2e})",
            rel_err, max_err
        );
    }
}
