//! # Lumina Compute
//!
//! Compute backend abstraction for the Lumina framework. This crate
//! provides a [`ComputeBackend`](backend::ComputeBackend) trait that isolates
//! the physics code from device-specific execution details.
//!
//! ## Available backends
//!
//! | Backend | Feature flag | Status |
//! |---------|-------------|--------|
//! | CPU (Rayon) | `cpu` (default) | Implemented |
//! | GPU (wgpu) | `gpu` | Stub |
//! | Distributed (MPI) | `distributed` | Stub |

pub mod backend;

#[cfg(feature = "cpu")]
pub mod cpu;

#[cfg(feature = "gpu")]
pub mod gpu;

#[cfg(feature = "distributed")]
pub mod distributed;

pub use backend::{BackendType, ComputeBackend, ComputeError, DeviceInfo};

#[cfg(feature = "cpu")]
pub use cpu::CpuBackend;
