# Tutorial 2: Nanorod & Plasmon Red-Shift

**Goal:** Compute the extinction spectrum of a gold nanorod (cylinder geometry) and observe how the longitudinal plasmon resonance red-shifts with increasing aspect ratio.

**Time:** ~10 minutes

---

## Background

A nanorod supports two distinct plasmon modes:

- **Transverse mode** — electrons oscillate perpendicular to the rod axis (~520 nm for Au, similar to a sphere)
- **Longitudinal mode** — electrons oscillate along the rod axis; resonance wavelength depends strongly on aspect ratio AR = length/diameter

Empirically, the longitudinal LSPR red-shifts by roughly 90–100 nm per unit of aspect ratio for gold nanorods in water. This tutorial reproduces that trend using the CDA.

---

## Step 1: Single Nanorod

Create `examples/nanorod_ar3.toml`:

```toml
# Tutorial 2 — Au nanorod, AR=3, in water
[simulation]
wavelengths   = { range = [400.0, 900.0], points = 100 }
environment_n = 1.33        # water

[[geometry.object]]
name          = "Au_Rod"
type          = "cylinder"
base_centre   = [0.0, 0.0, -30.0]   # rod centred at origin
axis          = [0.0, 0.0, 1.0]      # along z
length        = 60.0                  # nm
radius        = 10.0                  # nm  → AR = length/(2r) = 3
material      = "Au_JC"
dipole_spacing = 2.0                  # nm

[output]
directory     = "./output/tutorial_02/ar3"
save_spectra  = true
```

Run it:

```bash
cargo run --release --bin lumina-cli -- examples/nanorod_ar3.toml
```

---

## Step 2: Aspect Ratio Sweep

Create separate TOML files for each aspect ratio, or use a shell loop:

```bash
for AR in 2 3 4 5; do
  LENGTH=$((AR * 20))     # keep diameter = 20 nm
  cat > /tmp/rod_ar${AR}.toml << EOF
[simulation]
wavelengths   = { range = [400.0, 1000.0], points = 120 }
environment_n = 1.33

[[geometry.object]]
name          = "Au_Rod_AR${AR}"
type          = "cylinder"
base_centre   = [0.0, 0.0, -$((LENGTH/2)).0]
axis          = [0.0, 0.0, 1.0]
length        = ${LENGTH}.0
radius        = 10.0
material      = "Au_JC"
dipole_spacing = 2.0

[output]
directory     = "./output/tutorial_02/ar${AR}"
save_spectra  = true
EOF

  cargo run --release --bin lumina-cli -- /tmp/rod_ar${AR}.toml
done
```

---

## Step 3: Plot the Aspect Ratio Dependence

```python
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(figsize=(9, 5))
colours = ["#2196F3", "#4CAF50", "#FF9800", "#F44336"]

peak_wavelengths = []
for ar, colour in zip([2, 3, 4, 5], colours):
    df = pd.read_csv(f"output/tutorial_02/ar{ar}/spectra.csv")
    ax.plot(df["wavelength_nm"], df["extinction_nm2"],
            label=f"AR = {ar}", color=colour)
    # find longitudinal peak (ignore transverse ~520 nm)
    mask = df["wavelength_nm"] > 600
    idx  = df.loc[mask, "extinction_nm2"].idxmax()
    peak_wavelengths.append(df.loc[idx, "wavelength_nm"])

ax.set_xlabel("Wavelength (nm)")
ax.set_ylabel("Extinction Cross-Section (nm²)")
ax.set_title("Au nanorod: aspect ratio dependence of longitudinal LSPR")
ax.legend()
plt.tight_layout()
plt.savefig("nanorod_spectra.png", dpi=150)

# Aspect ratio vs peak wavelength
fig2, ax2 = plt.subplots(figsize=(5, 4))
ax2.plot([2, 3, 4, 5], peak_wavelengths, "ko-")
ax2.set_xlabel("Aspect Ratio (length / diameter)")
ax2.set_ylabel("Longitudinal LSPR peak (nm)")
ax2.set_title("Plasmon red-shift with aspect ratio")
# Linear fit
m, b = np.polyfit([2, 3, 4, 5], peak_wavelengths, 1)
ax2_x = np.array([1.5, 5.5])
ax2.plot(ax2_x, m * ax2_x + b, "r--", label=f"slope ≈ {m:.0f} nm/AR")
ax2.legend()
plt.tight_layout()
plt.savefig("nanorod_redshift.png", dpi=150)
print(f"Longitudinal LSPR red-shifts by ~{m:.0f} nm per unit aspect ratio")
```

Expected output: roughly 80–100 nm red-shift per unit aspect ratio in water.

---

## Step 4: Polarisation Dependence

The longitudinal mode is excited when the incident electric field is **parallel to the rod axis**. By default, Lumina uses an x-polarised plane wave propagating along z. To excite the longitudinal mode of a z-oriented rod you need x-polarisation — this works by default. If your rod is oriented along x, rotate it using the `axis` parameter:

```toml
axis = [1.0, 0.0, 0.0]   # rod along x → x-polarised wave excites longitudinal mode
```

Since v0.1.3, the GUI simulation panel includes a polarisation selector (x-pol / y-pol / circular) that controls the incident field direction.

---

## Step 5: Verify with GUI

The GUI geometry panel shows a 2D scatter preview of the dipole lattice. For a cylinder you should see a filled circle in the XY projection and a filled rectangle in the XZ/YZ projections.

```bash
cargo run --release --bin lumina-gui
```

Switch the projection plane dropdown between **XY**, **XZ**, **YZ** to inspect the lattice from different angles.

---

## Physical Interpretation

| Aspect Ratio | Approximate longitudinal LSPR (water, 20 nm diameter) |
|-------------|-------------------------------------------------------|
| 2 | ~640 nm |
| 3 | ~730 nm |
| 4 | ~820 nm |
| 5 | ~910 nm |

These values are consistent with published experimental data for chemically synthesised nanorods (Nikoobakht & El-Sayed 2003).

> **Accuracy note:** CDA with d=2 nm gives ~15–25% error on the peak cross-section magnitude due to the staircase boundary. The peak *wavelength* is typically reproduced to within 5–10 nm because it is governed by the aspect ratio geometry rather than the fine surface discretisation.

---

## Next Steps

- [Tutorial 3: Near-Field Map](./03-near-field.md) — map the hotspots at the nanorod tips
- [Tutorial 4: Chiral Helix](./04-chiral-helix.md) — a multi-turn helix with chiroptical response
