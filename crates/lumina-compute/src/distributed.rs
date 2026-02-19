//! Distributed (MPI) compute backend for HPC clusters.
//!
//! **Status: Stub for future implementation.**
//!
//! This backend will distribute the interaction matrix assembly and solve
//! across multiple nodes using MPI. The strategy is:
//!
//! - **Block-row distribution**: Each MPI rank owns a contiguous block of
//!   dipoles and assembles the corresponding rows of the interaction matrix.
//! - **Distributed GMRES**: Each rank performs its portion of the matrix-vector
//!   product; MPI_Allreduce handles the global dot products.
//! - **Hybrid MPI+Rayon**: Within each node, Rayon parallelises across cores.
//!
//! Gated behind the `distributed` feature flag.

// This module is intentionally left as a stub.
// Implementation is planned for v0.3+.
