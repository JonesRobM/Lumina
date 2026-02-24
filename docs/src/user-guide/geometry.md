# Geometry & Shapes

## Parametric Primitives

Define shapes directly in the TOML configuration. All dimensions are in nanometres.

### Sphere
```toml
type = "sphere"
centre = [0.0, 0.0, 0.0]
radius = 20.0
```

### Cylinder
```toml
type = "cylinder"
base_centre = [0.0, 0.0, 0.0]
axis = [0.0, 0.0, 1.0]
length = 60.0
radius = 10.0
```

### Cuboid
```toml
type = "cuboid"
centre = [0.0, 0.0, 0.0]
half_extents = [20.0, 10.0, 5.0]
```

### Ellipsoid
```toml
type = "ellipsoid"
centre = [0.0, 0.0, 0.0]
semi_axes = [20.0, 15.0, 10.0]
```

### Helix
```toml
type = "helix"
centre = [0.0, 0.0, 0.0]
helix_radius = 15.0
wire_radius = 4.0
pitch = 20.0
turns = 3
handedness = "right"
```

## File Import

Import atomic coordinates or mesh geometry into the GUI via the Geometry panel. Each atom or mesh vertex acts as a discrete dipole in the CDA.

Supported formats: `.xyz`, `.obj`.

### XYZ File Import

The `.xyz` format is a standard for atomic coordinates, widely used in computational chemistry (DFT outputs, MD snapshots, crystallographic generators). Lumina reads `.xyz` files and treats each atom as an interacting point dipole.

**Format:**
```
N
Comment line
Element x y z
Element x y z
...
```

Coordinates are in Angstroms. Lumina automatically converts to nanometres (1 A = 0.1 nm).

**Auto-detected parameters (v0.2.1):**

When importing an `.xyz` file, Lumina automatically detects the nearest-neighbour distance from the atomic positions. This distance is used for:

- **FCD cell size** — the volume over which the Green's tensor is integrated for near-field interactions
- **Clausius-Mossotti volume** — \\(V = d_\text{nn}^3\\), the effective dipole volume for computing the polarisability

This avoids the need to manually set a dipole spacing for atomistic structures, and ensures the FCD threshold (\\(2 d_\text{nn}\\)) correctly identifies near-field pairs.

### Supplied Example: Au Mackay Icosahedron

The file `examples/au_icosahedron_k6.xyz` contains a 923-atom Au Mackay icosahedron with k=6 shells. This is a close-packed structure with icosahedral (\\(I_h\\)) symmetry, commonly used as a model for metallic nanoparticle clusters.

**Properties:**
- 923 atoms (follows the Mackay formula: \\(N(k) = (10k^3 + 15k^2 + 11k + 3) / 3\\))
- Nearest-neighbour distance: 2.884 A (0.2884 nm), derived from Au lattice constant \\(a = 4.078\\) A
- Particle diameter: ~1.6 nm (bounding sphere)
- Full \\(I_h\\) symmetry: 12 vertices, 30 edges, 20 faces of the outer icosahedral shell

**GUI workflow:**
1. Geometry panel → **Import File** → **Set folder...** → select the `examples/` directory
2. Select `au_icosahedron_k6.xyz` from the dropdown
3. The dipole scatter preview shows the icosahedral symmetry
4. Materials → **Gold (J&C)**, Simulation → set wavelength range → **Run Simulation**

### OBJ Mesh Import

Wavefront `.obj` files define a triangulated surface mesh. Lumina fills the enclosed volume with a cubic lattice of dipoles using ray-casting (counting z-axis ray intersections to determine inside/outside).

```
# Simple cube
v 0.0 0.0 0.0
v 10.0 0.0 0.0
...
f 1 2 3
f 1 3 4
...
```

### Directory-Scoped File Browser (v0.2.1)

The GUI geometry panel offers a folder-based workflow for managing multiple structure files:

1. Click **Set folder...** to select a directory
2. All `.xyz` and `.obj` files appear in a dropdown
3. Select a file to load it immediately
4. **Refresh** re-scans the directory

The traditional file dialog buttons remain available for one-off picks.

## Discretisation

Each parametric shape is filled with a centred cubic lattice of dipole positions at the specified `dipole_spacing`. The spacing must satisfy the CDA validity criterion: \\(|m| k d < 1\\).

For `.xyz` imports, no discretisation is needed — the atom positions are used directly as dipole sites.

### Centred Lattice

The lattice is centred on the midpoint of the shape's bounding box, using integer offsets in each direction:

\\[
\mathbf{r}_{ijk} = \mathbf{r}_{\text{centre}} + d \cdot (i, j, k), \quad i, j, k \in \mathbb{Z}
\\]

Only points satisfying `shape.contains(r)` are retained.

Centring the lattice on the shape centre ensures that the dipole distribution is symmetric for symmetric shapes, regardless of the spacing. This is critical for metallic particles: an asymmetric lattice (e.g. anchored to a bounding box corner) can introduce 100–500% errors due to artificial polarisation charge at the uneven surface.

### Spacing Guidelines

| Material type | Recommended \\(d\\) | Expected accuracy |
|--------------|---------------------|-------------------|
| Dielectric (\\(\epsilon_2 > 0.5\\)) | \\(R/5\\) to \\(R/3\\) | \\(< 15\%\\) |
| Moderate metal (\\(|\epsilon_1|/\epsilon_2 < 3\\)) | \\(R/5\\) | \\(< 30\%\\) |
| Strong metal (\\(|\epsilon_1|/\epsilon_2 > 5\\)) | \\(\leq R/10\\) | Limited by staircase |

See [v0.1 Validation](../theory/validation.md) for convergence benchmarks.
