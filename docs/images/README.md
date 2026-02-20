# Lumina Documentation Images

This directory contains images and graphics for the README and documentation.

## Current Files

- **gui_overview.svg** - SVG mockup of the Lumina GUI interface
- **spectra_plot.svg** - Example extinction/absorption/scattering spectra plot
- **nearfield_heatmap.svg** - Example near-field intensity map

## Generating Actual Screenshots

To replace the SVG placeholders with real screenshots from the GUI:

### 1. Run the GUI

```bash
cargo run --release --bin lumina-gui
```

### 2. Set up a gold sphere simulation

- Geometry: Sphere, R=10nm, spacing=2nm
- Material: Au (J&C)
- Simulation: 400-800nm, 50 points

### 3. Capture screenshots

**GUI Overview** (`gui_overview.png`):
- Show the Geometry panel with sphere parameters and dipole preview
- Make sure the dipole count (515) is visible

**Spectra Plot** (`spectra_plot.png`):
- Navigate to Results panel after running simulation
- Ensure all three curves (extinction, absorption, scattering) are visible
- Show the legend

**Near-field Heatmap** (`nearfield_heatmap.png`):
- Navigate to Results → Near-field tab
- Capture the |E|² heatmap at peak wavelength
- Include the color scale

### 4. Export from simulation runs

You can also generate plots programmatically using the CLI:

```bash
cargo run --release --bin lumina-cli -- examples/gold_sphere.toml
```

Then use a plotting tool (Python/matplotlib, gnuplot, etc.) to visualize the CSV output.

## Image Guidelines

- **Format**: PNG for screenshots (better compression for photos)
- **Size**: Max 1200px width for README images
- **Resolution**: 144 DPI for retina displays
- **File size**: Keep under 500KB per image

## SVG vs PNG

The current SVG files are:
- ✅ Small file size (~10KB each)
- ✅ Scalable for any display
- ✅ Version control friendly (text-based)
- ✅ Professional appearance

Real PNG screenshots provide:
- ✅ Authentic look and feel
- ✅ Actual rendered results
- ⚠️ Larger file sizes (100-500KB)
- ⚠️ Fixed resolution

**Recommendation**: Keep both! SVGs for clean diagrams, PNGs for actual GUI screenshots when available.
