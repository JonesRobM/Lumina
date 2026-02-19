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

## File Import

Import atomic coordinates or mesh geometry:

```toml
geometry_file = "structure.xyz"
```

Supported formats: `.xyz`, `.cif`, `.obj`.

## Discretisation

Each shape is filled with a cubic lattice of dipole positions at the specified `dipole_spacing`. The spacing must satisfy the CDA validity criterion: \\(|m| k d < 1\\).
