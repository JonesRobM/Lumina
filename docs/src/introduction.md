# Lumina

**Lumina** is a high-performance Rust implementation of the Coupled Dipole Approximation (CDA), engineered to compute the linear and nonlinear optical responses of nanostructured, optically active materials.

## What It Does

Lumina bridges the gap between ab initio material properties and macroscopic electromagnetic observables. Given a nanostructure geometry and its material properties, it computes:

- **Extinction, absorption, and scattering cross-sections** as a function of wavelength.
- **Near-field intensity maps** showing the electric field enhancement around the structure.
- **Dipole moment distributions** for detailed analysis of the optical response.

## Who It's For

- **Researchers** studying plasmonic nanoparticles, chiral nanostructures, and metamaterials.
- **Engineers** designing optical sensors, filters, and enhancers.
- **Students** learning computational electromagnetics and nanophotonics.

## Key Features

- **Modular solver architecture** — the CDA is the first solver; future methods (BEM, FDTD, T-matrix) plug in via the same trait.
- **GPU acceleration** — compute-heavy operations offload to the GPU via wgpu.
- **Interactive GUI** — real-time parameter sweeping and visualisation.
- **HPC-ready CLI** — headless operation via TOML configuration files.
- **Extensible materials** — tabulated data (Johnson & Christy), cubic spline interpolation, and ab initio import.
