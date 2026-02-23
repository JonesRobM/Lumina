# Quick Start

This guide walks through computing the extinction spectrum of a gold nanosphere.

## 1. Create a Configuration File

Save the following as `my_sphere.toml`:

```toml
[simulation]
wavelengths = { range = [400.0, 900.0], points = 100 }
environment_n = 1.0

[[geometry.object]]
name = "Gold_Sphere"
type = "sphere"
centre = [0.0, 0.0, 0.0]
radius = 20.0
material = "Au_JC"
dipole_spacing = 2.0

[output]
directory = "./output"
save_spectra = true
save_json = true
```

## 2. Run via CLI

```bash
cargo run --release -p lumina-cli -- run my_sphere.toml
```

To use GPU acceleration (requires `--features gpu` build):

```bash
cargo run --release -p lumina-cli --features gpu -- run my_sphere.toml
```

## 3. Or Use the GUI

Launch the GUI and configure the same parameters interactively:

```bash
cargo run --release -p lumina-gui
```

Select a shape, material, wavelength range, and click **Run Simulation**. Results appear in the Results panel with interactive spectra, near-field heatmaps, and far-field polar plots.

## 4. Examine the Results

The output directory will contain:

- `spectra.csv` — wavelength, extinction, absorption, and scattering cross-sections
- `spectra.json` — same data in JSON format (when `save_json = true`)
- `near_field/` — field intensity maps (when `save_near_field = true`)
