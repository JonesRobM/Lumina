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

Import atomic coordinates or mesh geometry:

```toml
geometry_file = "structure.xyz"
```

Supported formats: `.xyz`, `.cif`, `.obj`.

## Discretisation

Each shape is filled with a centred cubic lattice of dipole positions at the specified `dipole_spacing`. The spacing must satisfy the CDA validity criterion: \\(|m| k d < 1\\).

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
