# Tutorial 3: Near-Field Intensity Map

**Goal:** Compute the electric field enhancement |E|²/|E₀|² around a gold nanosphere at its plasmon resonance and visualise the result as a 2D heatmap.

**Time:** ~10 minutes

---

## Background

When a nanoparticle is at resonance, the electromagnetic field is strongly confined near its surface — the so-called *plasmonic hotspot*. The near-field intensity enhancement |E|²/|E₀|² can reach factors of 10²–10⁴ at sharp tips and gaps, enabling applications such as SERS (surface-enhanced Raman scattering) and nanoscale optical trapping.

In the CDA, the near field at an observation point **r** is computed from the solved dipole moments {**p**_i}:

$$\mathbf{E}(\mathbf{r}) = \mathbf{E}_\text{inc}(\mathbf{r}) + \sum_i \mathbf{G}(\mathbf{r}, \mathbf{r}_i)\,\mathbf{p}_i$$

Lumina evaluates this on a user-defined observation plane.

---

## Step 1: Configure Near-Field Output

Create `examples/nearfield_tutorial.toml`:

```toml
# Tutorial 3 — near-field map around a 15 nm Au sphere at resonance
[simulation]
wavelengths   = { range = [510.0, 530.0], points = 5 }  # zoom around LSPR
environment_n = 1.0

[[geometry.object]]
name          = "Au_Sphere_NF"
type          = "sphere"
centre        = [0.0, 0.0, 0.0]
radius        = 15.0          # nm
material      = "Au_JC"
dipole_spacing = 2.0          # nm

[output]
directory       = "./output/tutorial_03"
save_spectra    = true
save_near_field = true        # ← enables near-field computation
save_dipoles    = false
```

> **Performance:** Near-field evaluation on a 40×40 grid adds roughly 10–20% overhead per wavelength. Lumina evaluates the near field only at the wavelength of peak extinction.

Run:

```bash
cargo run --release --bin lumina-cli -- examples/nearfield_tutorial.toml
```

---

## Step 2: Near-Field Output Files

After the run, `output/tutorial_03/` contains:

```
spectra.csv
near_field/
  near_field_xy_peak.csv    ← XY plane at z=0, |E|²/|E₀|²
```

The CSV format:

```
x_nm,y_nm,intensity
-30.0,-30.0,1.023
-30.0,-28.5,1.089
…
0.0,0.0,4.21e2
…
30.0,30.0,1.018
```

The observation grid spans ±2R (here ±30 nm) in both x and y with 40 points per axis (1600 total).

---

## Step 3: Plot the Heatmap

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle

df = pd.read_csv("output/tutorial_03/near_field/near_field_xy_peak.csv")

# Pivot to 2D grid
x = np.unique(df["x_nm"])
y = np.unique(df["y_nm"])
Z = df["intensity"].values.reshape(len(y), len(x))

fig, ax = plt.subplots(figsize=(6, 5))
im = ax.pcolormesh(x, y, Z,
                   cmap="inferno",
                   norm=LogNorm(vmin=1, vmax=Z.max()))
cbar = fig.colorbar(im, ax=ax, label="|E|² / |E₀|²")

# Outline the sphere
sphere = Circle((0, 0), radius=15, fill=False,
                edgecolor="white", linewidth=1.5, linestyle="--")
ax.add_patch(sphere)

ax.set_xlabel("x (nm)")
ax.set_ylabel("y (nm)")
ax.set_title("Near-field intensity — Au sphere R=15 nm at LSPR")
ax.set_aspect("equal")
plt.tight_layout()
plt.savefig("nearfield_map.png", dpi=150)
```

---

## Step 4: Interpret the Map

For an x-polarised plane wave and a sphere:

| Feature | Location | Explanation |
|---------|----------|-------------|
| Enhancement lobes | ±x poles | Field aligned with incident polarisation |
| Suppression | ±y | Dipole radiation pattern — nodes at 90° |
| Enhancement factor | ~100–1000× | Scales with dipole oscillator strength |

The enhancement pattern should be **dipolar** (two lobes) and symmetric about the x-axis. Any asymmetry indicates a discretisation artefact.

> **Field inside the sphere:** Grid points inside the nanoparticle are skipped during near-field evaluation (they would require the singular self-interaction). The CSV marks these as `NaN`.

---

## Step 5: Cross-Section Plane Sweep

To map the xz plane (along the propagation direction), modify the observation plane in the configuration — this is currently set via code when using the Rust API:

```rust
use lumina_core::fields::{NearFieldPlane, compute_near_field};

let plane = NearFieldPlane {
    centre:    [0.0, 0.0, 0.0],
    normal:    [0.0, 1.0, 0.0],   // xz plane (y=0)
    half_width: 40.0,
    resolution: 50,
};
let nf = compute_near_field(&plane, &dipoles, wavelength, &params)?;
```

CLI-based plane selection is planned for a future release.

---

## Step 6: GUI Near-Field View

In the GUI, the near-field heatmap renders automatically after a simulation:

1. Run the simulation in the GUI
2. Navigate to **Results → Near-Field** tab
3. The heatmap updates at the peak extinction wavelength
4. Hover over any pixel to read the local |E|²/|E₀|² value

The GUI uses a 16-bin hot colourmap (black → red → yellow → white) on a log scale.

---

## Understanding Enhancement Magnitudes

| System | Typical max |E|²/|E₀|² | Application |
|--------|------------------------|-------------|
| Au sphere (R=15 nm) | ~100–500 | Fluorescence enhancement |
| Au nanorod tip | ~1 000–5 000 | SERS |
| Dimer gap (~2 nm) | ~10⁴–10⁶ | Single-molecule SERS |

CDA systematically **underestimates** the peak enhancement at sharp features because a coarse dipole lattice cannot resolve the field singularity. For quantitative gap plasmonics, a BEM or FDTD calculation is more appropriate. The CDA remains reliable for the qualitative spatial distribution and for features larger than 2–3 dipole spacings.

---

## Next Steps

- [Tutorial 4: Chiral Helix](./04-chiral-helix.md) — near-field coupling between helix turns
- [Tutorial 5: Custom Material](./05-custom-material.md) — use a dielectric with your own ε values
