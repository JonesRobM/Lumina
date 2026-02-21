# Tutorial 4: Chiral Helix

**Goal:** Simulate a gold helical nanoparticle and observe chiroptical extinction — the differential absorption of left- and right-circularly polarised light.

**Time:** ~15 minutes

---

## Background

A helix has no mirror symmetry and therefore interacts differently with left-circularly polarised (LCP) and right-circularly polarised (RCP) light. The differential extinction:

$$\Delta C_\text{ext} = C_\text{ext}^\text{LCP} - C_\text{ext}^\text{RCP}$$

is non-zero and reverses sign for the mirror-image helix (enantiomer). This is the hallmark of optical chirality and is relevant to chiral sensing, metamaterial design, and optically active photonics.

In the CDA, chiroptical response emerges naturally from the geometry — no extra physics is required beyond the standard dyadic Green's function.

---

## The Helix Geometry

Lumina's helix primitive is defined by:

| Parameter | Description |
|-----------|-------------|
| `centre` | Centre of the helix coil |
| `helix_radius` | Radius of the helical coil (nm) |
| `wire_radius` | Radius of the wire cross-section (nm) |
| `pitch` | Rise per full turn (nm) |
| `turns` | Number of turns |
| `handedness` | `"left"` or `"right"` |

---

## Step 1: Single Right-Handed Helix

Create `examples/helix_rh.toml`:

```toml
# Tutorial 4 — right-handed gold helix
[simulation]
wavelengths   = { range = [400.0, 900.0], points = 100 }
environment_n = 1.0

[[geometry.object]]
name          = "Helix_RH"
type          = "helix"
centre        = [0.0, 0.0, 0.0]
helix_radius  = 15.0          # nm
wire_radius   = 4.0           # nm
pitch         = 20.0          # nm per turn
turns         = 3
handedness    = "right"
material      = "Au_JC"
dipole_spacing = 2.0          # nm

[output]
directory     = "./output/tutorial_04/rh"
save_spectra  = true
```

Run:

```bash
cargo run --release --bin lumina-cli -- examples/helix_rh.toml
```

---

## Step 2: Left-Handed Enantiomer

Create `examples/helix_lh.toml` — identical except `handedness = "left"`:

```toml
# Tutorial 4 — left-handed gold helix (enantiomer)
[simulation]
wavelengths   = { range = [400.0, 900.0], points = 100 }
environment_n = 1.0

[[geometry.object]]
name          = "Helix_LH"
type          = "helix"
centre        = [0.0, 0.0, 0.0]
helix_radius  = 15.0
wire_radius   = 4.0
pitch         = 20.0
turns         = 3
handedness    = "left"
material      = "Au_JC"
dipole_spacing = 2.0

[output]
directory     = "./output/tutorial_04/lh"
save_spectra  = true
```

```bash
cargo run --release --bin lumina-cli -- examples/helix_lh.toml
```

---

## Step 3: Plot the Chiroptical Response

Lumina currently outputs spectra for a linearly polarised incident wave. Circularly polarised responses are obtained by computing LCP and RCP as linear combinations:

$$\mathbf{E}_\text{LCP} = \frac{1}{\sqrt{2}}(\hat{x} + i\hat{y})$$
$$\mathbf{E}_\text{RCP} = \frac{1}{\sqrt{2}}(\hat{x} - i\hat{y})$$

For a single linearly polarised result, a useful proxy for the circular dichroism is to compare the **right-handed** and **left-handed** enantiomers under the same (x-polarised) excitation. The spectra will differ in peak position and cross-section magnitude if the geometry has helical anisotropy.

```python
import pandas as pd
import matplotlib.pyplot as plt

rh = pd.read_csv("output/tutorial_04/rh/spectra.csv")
lh = pd.read_csv("output/tutorial_04/lh/spectra.csv")

fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

# Top: absolute extinction spectra
axes[0].plot(rh["wavelength_nm"], rh["extinction_nm2"],
             label="Right-handed", color="#E74C3C")
axes[0].plot(lh["wavelength_nm"], lh["extinction_nm2"],
             label="Left-handed",  color="#3498DB", linestyle="--")
axes[0].set_ylabel("C_ext (nm²)")
axes[0].set_title("Gold helix: enantiomer comparison")
axes[0].legend()

# Bottom: differential extinction ΔC = C(RH) − C(LH)
delta = rh["extinction_nm2"].values - lh["extinction_nm2"].values
axes[1].plot(rh["wavelength_nm"], delta, color="#27AE60")
axes[1].axhline(0, color="gray", linewidth=0.8, linestyle=":")
axes[1].set_xlabel("Wavelength (nm)")
axes[1].set_ylabel("ΔC_ext (nm²)")
axes[1].set_title("Differential extinction (chiroptical signal)")

plt.tight_layout()
plt.savefig("helix_chiroptical.png", dpi=150)
```

A non-zero differential signal confirms that the two enantiomers respond differently to x-polarised light — this is a geometric consequence of the broken mirror symmetry.

---

## Step 4: Multi-Object Chiral Dimer

The `examples/chiral_dimer.toml` in the repository shows a pair of nanorods twisted by 45°, the simplest chiral dimer:

```toml
[[geometry.object]]
name       = "Rod_Bottom"
type       = "cylinder"
base_centre = [0.0, 0.0, -15.0]
axis       = [1.0, 0.0, 0.0]
length     = 60.0
radius     = 10.0
material   = "Au_JC"
dipole_spacing = 2.0

[[geometry.object]]
name       = "Rod_Top"
type       = "cylinder"
base_centre = [0.0, 0.0, 15.0]
axis       = [0.707, 0.707, 0.0]   # 45° rotation
length     = 60.0
radius     = 10.0
material   = "Au_JC"
dipole_spacing = 2.0
```

The inter-rod gap is 30 nm (gap = 30 nm − 10 nm − 10 nm = 10 nm centre-surface). CDA-level coupling is valid for gaps > ~5 nm.

---

## Step 5: Inspect the Dipole Lattice

The GUI geometry panel lets you verify the helix discretisation visually:

```bash
cargo run --release --bin lumina-gui
```

Select **Helix** from the shape dropdown and set:
- Helix radius: 15
- Wire radius: 4
- Pitch: 20
- Turns: 3

Switch between XY/XZ/YZ projections to confirm the coil structure.

> The XY projection of a helix looks like a **filled annulus** (ring) since all turns project onto the same circle. The XZ projection shows the full helical ribbon.

---

## Physical Notes

- A 3-turn Au helix with ~15 nm coil radius is typical of DNA-origami assembled plasmonic helices (Kuzyk et al., *Nature* 2012).
- The CDA works well for helices because the wire cross-section (radius ~4 nm, d=2 nm → ~12 dipoles per cross-section) is small enough that near-field coupling between turns dominates over surface staircase artefacts.
- Chiroptical signals from single helices are weak; experimental spectra typically require large arrays.

---

## Next Steps

- [Tutorial 5: Custom Material](./05-custom-material.md) — supply your own ε(λ) data
- [Validation](../theory/validation.md) — accuracy benchmarks and known limitations
