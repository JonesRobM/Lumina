# Complex Assembly Builder Design

**Date:** 2026-03-02
**Status:** Approved
**Scope:** `lumina-geometry`, `lumina-core`, `lumina-cli`, `lumina-gui`

---

## Problem

The GUI supports only a single object with a single global material. Building
dimers, trimers, nanorod forests, or multi-species metallic assemblies requires
hand-editing TOML for the CLI — and even there, no transform (rotation/scale)
support exists. XYZ imports ignore element labels, assigning one material to all
atoms regardless of species.

---

## Goals

1. GUI scene panel with an object list: add/remove/duplicate objects, each with
   its own shape, material, dipole spacing, and full transform.
2. Full per-object transform: position (nm), ZYX Euler rotation (degrees),
   uniform scale.
3. XYZ multi-species mapping: manual element → material table in the GUI and
   TOML.
4. GUI ↔ TOML parity: export scene to TOML, import TOML into GUI.
5. CLI uses the same `SceneSpec` type — removes ~60 lines of duplicated
   assembly logic from `runner.rs`.

---

## Non-Goals

- Interactive drag-and-drop positioning (numeric inputs only).
- Anisotropic scale or shear transforms.
- Arrangement presets (dimer wizard, forest generator) — deferred to v0.3+.
- Undo/redo history in the GUI.

---

## Design

### 1. New Shared Types (`lumina-geometry/src/scene.rs`)

```rust
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ObjectTransform {
    pub position:     [f64; 3],   // nm offset (x, y, z)
    pub rotation_deg: [f64; 3],   // ZYX Euler angles in degrees
    pub scale:        f64,        // uniform scale (1.0 = identity)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ObjectSpec {
    pub name:           String,
    pub shape:          ShapeSpec,              // enum: primitive or file path
    pub material:       String,                 // "Au_JC", "Ag_JC", etc.
    pub dipole_spacing: f64,                    // nm
    pub transform:      ObjectTransform,
    pub species_map:    HashMap<String, String>,// element label → material (XYZ only)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SceneSpec {
    pub objects: Vec<ObjectSpec>,
}
```

`SceneSpec::build_dipoles(wavelength_nm, environment_n)` is the single entry
point that flattens the scene into `Vec<Dipole>` for the CDA solver.

### 2. Transform Pipeline (per object, per wavelength)

1. **Discretise** shape at `dipole_spacing` → local-space positions (centred at origin).
2. **Scale** — multiply by `transform.scale`.
3. **Rotate** — apply 3×3 matrix from ZYX Euler: `R = Rz(γ) · Ry(β) · Rx(α)`.
4. **Translate** — add `transform.position`.
5. **Resolve material** — for each point, use `species_map[label]` if present,
   else `object.material` → `ε(ω)` → Clausius-Mossotti → `Dipole`.

Steps 1–4 are geometry-only (wavelength-independent). Step 5 is per-wavelength.
The CDA solver receives a flat `Vec<Dipole>` — no solver changes required.

### 3. GUI Scene Panel (`lumina-gui`)

Replaces `GeometryPanel` + `MaterialsPanel` with a unified `ScenePanel`.

**Layout:**
- Left column: scrollable object list (name + shape icon + material tag).
  Buttons: `+ Add`, `Duplicate`, `Delete`, `↑ ↓`.
- Right column: property editor for the selected object.
  - Name, Shape selector, Material dropdown, Dipole spacing.
  - Transform: Position (x/y/z nm), Rotation (Rx/Ry/Rz °), Scale.
  - XYZ species map table (shown only for file imports).
- Top bar: `Import scene (TOML)` / `Export scene (TOML)`.

`app.rs` replaces `geometry_state: GeometryPanel` and `materials_state:
MaterialsPanel` with `scene_panel: ScenePanel`. The simulation launch reads
`scene_panel.scene_spec()` → `SceneSpec` → `build_dipoles(...)`.

### 4. XYZ Multi-Species Mapping

When an XYZ file is loaded, detected element labels are shown in a mapping
table in the object properties. Each row: label → material dropdown. In TOML:

```toml
[geometry.object.species_map]
Au = "Au_JC"
Ag = "Ag_JC"
C  = "SiO2_Palik"
```

`species_map` overrides `material` per atom. Empty map → uniform material.

### 5. TOML Parity

- **Export:** `toml::to_string(&scene_spec)` — one call, no conversion layer.
- **Import:** `toml::from_str::<SceneSpec>(...)` — deserialises directly into
  `ScenePanel` object list.
- **CLI migration:** `config.rs` `[[geometry.object]]` array deserialises into
  `SceneSpec`; `runner.rs` calls `scene.build_dipoles(...)` replacing the
  manual per-object loop.

TOML example (Au-Ag dimer):
```toml
[[geometry.object]]
name = "Left_sphere"
shape = { type = "sphere", radius = 20.0 }
material = "Au_JC"
dipole_spacing = 2.0
position = [-25.0, 0.0, 0.0]
rotation_deg = [0.0, 0.0, 0.0]
scale = 1.0

[[geometry.object]]
name = "Right_sphere"
shape = { type = "sphere", radius = 20.0 }
material = "Ag_JC"
dipole_spacing = 2.0
position = [25.0, 0.0, 0.0]
rotation_deg = [0.0, 0.0, 0.0]
scale = 1.0
```

---

## Testing

| Test | Description | Pass criterion |
|------|-------------|----------------|
| `test_object_transform_identity` | Zero rotation + zero translation | Positions unchanged to 1e-12 |
| `test_euler_rotation_x90` | 90° rotation about X: Y→Z | Correct coordinate permutation |
| `test_dimer_no_overlap` | Two 20nm spheres, 10nm gap | No dipole pair closer than 8nm |
| `test_scene_toml_round_trip` | Serialise → TOML → deserialise | All fields identical |
| `test_xyz_species_map_applies` | XYZ Au+Ag, mapped to Au_JC/Ag_JC | Dipole polarisabilities differ |
| `test_scene_build_dipoles_count` | Dimer: total = sum of individuals | Count correct |

---

## Files Affected

| File | Change |
|------|--------|
| `crates/lumina-geometry/src/scene.rs` | New: `ObjectTransform`, `ObjectSpec`, `SceneSpec`, `build_dipoles()` |
| `crates/lumina-geometry/src/lib.rs` | Export `scene` module |
| `crates/lumina-geometry/Cargo.toml` | Add `serde` + `toml` features |
| `crates/lumina-cli/src/config.rs` | Migrate `[[geometry.object]]` to use `SceneSpec` |
| `crates/lumina-cli/src/runner.rs` | Replace manual assembly loop with `scene.build_dipoles()` |
| `crates/lumina-gui/src/panels/scene.rs` | New: `ScenePanel` replacing geometry + materials panels |
| `crates/lumina-gui/src/app.rs` | Replace `geometry_state`/`materials_state` with `scene_panel` |
| `crates/lumina-geometry/src/tests/` | 6 new tests |
