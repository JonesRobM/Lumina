//! GPU compute backend via wgpu.
//!
//! **Status: Stub for future implementation.**
//!
//! This backend will use `wgpu` compute shaders for embarrassingly parallel
//! operations (Green's function matrix assembly) and delegate dense linear
//! algebra to platform-specific libraries (cuSOLVER on CUDA, Metal Performance
//! Shaders on macOS) as they become available.
//!
//! # Design
//!
//! - Matrix assembly: WGSL compute shader that fills the 3x3 blocks of the
//!   interaction matrix in parallel.
//! - Matrix-vector product: Shader-based reduction for GMRES iterations.
//! - Dense solve: Delegated to cuSOLVER/cuBLAS via feature-gated FFI when
//!   available; falls back to CPU otherwise.

// This module is intentionally left as a stub.
// Gated behind the `gpu` feature flag.
