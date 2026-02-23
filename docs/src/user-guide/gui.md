# GUI Dashboard

The Lumina GUI provides an interactive environment for configuring, running, and visualising simulations.

## Launching

```bash
# CPU only
cargo run --release -p lumina-gui

# With GPU acceleration toggle
cargo run --release -p lumina-gui --features gpu
```

## Panels

### Geometry

Configure shapes (sphere, cylinder, cuboid, ellipsoid, helix), import files (.xyz, .obj), adjust dipole spacing, and preview the dipole lattice as a 2D scatter plot (XY/XZ/YZ projections).

### Materials

Select from built-in materials (Au, Ag, Cu via Johnson & Christy; TiO\(_2\), SiO\(_2\) via Palik) or define custom dielectric functions. View the interpolated \(\epsilon(\lambda)\) plot showing \(\epsilon_1\) and \(\epsilon_2\) curves.

### Simulation

Set the wavelength range, medium refractive index, and solver parameters. Additional controls:

- **Incident polarisation** — select x-pol, y-pol, or circular
- **Circular dichroism** — checkbox to compute \(\Delta C_\text{ext}\) (requires two solver calls per wavelength)
- **GPU acceleration** — checkbox to enable GPU matvec (only visible when built with `--features gpu`)

Launch the simulation and monitor progress via the progress bar.

### Results

View interactive plots powered by `egui_plot`:

- **Spectra** — extinction, absorption, and scattering cross-sections vs wavelength (with optional CD curve)
- **Near-field heatmap** — \(|E|^2 / |E_0|^2\) on the XY plane at peak extinction, log-scale colourmap
- **Far-field polar plot** — E-plane and H-plane radiation pattern cuts

Export results to CSV via the file dialog (native file picker via `rfd`).
