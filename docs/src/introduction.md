# Lumina

**Lumina** is a high-performance Rust implementation of the Coupled Dipole Approximation (CDA), engineered to compute the linear and nonlinear optical responses of nanostructured, optically active materials.

## What It Does

Lumina bridges the gap between ab initio material properties and macroscopic electromagnetic observables. Given a nanostructure geometry and its material properties, it computes:

- **Extinction, absorption, and scattering cross-sections** as a function of wavelength.
- **Near-field intensity maps** showing the electric field enhancement around the structure.
- **Far-field radiation patterns** on the unit sphere.
- **Circular dichroism** (differential extinction for left/right circularly polarised light).
- **Dipole moment distributions** for detailed analysis of the optical response.

## Who It's For

- **Researchers** studying plasmonic nanoparticles, chiral nanostructures, and metamaterials.
- **Engineers** designing optical sensors, filters, and enhancers.
- **Students** learning computational electromagnetics and nanophotonics.

## Key Features

- **Modular solver architecture** — the CDA is the first solver; future methods (BEM, FDTD, T-matrix) plug in via the same trait.
- **GPU acceleration** — matrix-vector products offload to the GPU via wgpu compute shaders, with automatic CPU fallback.
- **GMRES iterative solver** — matrix-free GMRES(m) with Rayon-parallel matrix assembly for systems with thousands of dipoles.
- **Interactive GUI** — real-time parameter sweeping and visualisation with spectra, near-field heatmaps, and far-field polar plots.
- **HPC-ready CLI** — headless operation via TOML configuration files with JSON and CSV export.
- **Extensible materials** — tabulated data (Johnson & Christy Au/Ag/Cu, Palik TiO\\(_2\\)/SiO\\(_2\\)), cubic spline interpolation, and custom dielectrics.
