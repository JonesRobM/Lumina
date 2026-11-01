# Changelog

All notable changes to LuminaCDA are recorded here. Versions follow
[Semantic Versioning](https://semver.org/spec/v2.0.0.html); while the project is
pre-1.0, minor versions may contain breaking changes.

## v0.7.0

A correctness and performance release. Three changes alter numerical results, and
several public types changed shape. Read "Results that change" before comparing
output against an earlier version.

### Results that change

**FCD cell size now follows the dipole spacing (command-line interface only).**

The Filtered Coupled Dipole correction averages the Green's tensor over a cubic
cell and uses twice that cell size as its near-field cutoff. The command-line
interface previously left this at a fixed 3.0 nm regardless of the `dipole_spacing`
set in the job file, so every job with a different spacing integrated over a cell
that did not exist. The cell size is now taken from the geometry.

Measured effect on extinction for a gold sphere of radius 12 nm at 600 nm:

| dipole_spacing | before  | after   | change |
|----------------|---------|---------|--------|
| 1.90 nm        | 4.08e1  | 3.72e1  | 9.7 %  |
| 2.50 nm        | 4.77e1  | 4.60e1  | 3.7 %  |
| 4.00 nm        | 3.63e1  | 3.41e1  | 6.5 %  |

Jobs using a 3.0 nm spacing are unaffected, as are all runs made through the
graphical interface, which already passed the correct value.

**One-dimensional periodic (chain) calculations.**

Three separate defects in the chain Green's function are fixed:

1. The spectral sum's cross coefficient carried a spurious factor of 1/rho. It
   derives from `d_x d_alpha K0 = (i g_m)(-kappa K1(kappa rho) rho_hat_alpha)` and
   should carry a bare `kappa*K1`. Off-axis pairs were wrong by a factor of rho
   expressed in nanometres, so by 10x to 50x for typical geometries.
2. The number of spectral modes did not depend on rho, although the terms decay as
   `exp(-2*pi*m*rho/a1)`. Small transverse separations were under-resolved.
3. Pairs lying on the chain axis returned an exactly zero tensor. The Bessel
   representation is not valid at rho = 0, and no other routine covered that case,
   so a chain of on-axis dipoles experienced no lattice coupling at all and a
   discretised chain lost its self-image sum.

Any previous result from a one-dimensional periodic calculation should be recomputed.

**A Sommerfeld substrate combined with a periodic lattice is now refused.**

This combination previously ran and returned a number. The Sommerfeld reflected
tensor is evaluated pair by pair with no lattice sum, while the direct term is
Ewald-summed, so the result mixed a periodic and a non-periodic contribution. It now
returns `SolverError::UnsupportedConfiguration`. Use the Fresnel substrate
(`use_sommerfeld = false`) for periodic structures; its image-dipole term is
Ewald-summed correctly.

### Breaking API changes

- `Dipole` gains a `cell_size` field, and `Dipole::isotropic` takes it as a third
  argument. The cell size belongs to the dipole because FCD averages over the source
  cell, so a scene mixing dipole spacings needs one value per dipole.
- `CdaSolver::cell_size` is removed, as are the `cell_size` parameters of
  `CdaSolver::with_fcd`, `with_incident` and `with_backend`.
- `SolverError` is now `#[non_exhaustive]` and gains an `UnsupportedConfiguration`
  variant. Matches on this enum need a catch-all arm.
- `solver::cda::assembly` is now crate-private. `assemble_interaction_matrix` and
  `matvec_on_the_fly` were reachable from outside the crate, which allowed the
  substrate and lattice compatibility check to be bypassed.

### Performance

- The direct solver is now used up to 3000 dipoles rather than 1000, and is gated on
  the matrix memory budget. Between N = 1021 and N = 1791 it measured 8.2x to 13.5x
  faster than the iterative path, which had been slower than a direct solve
  throughout that range.
- Block-Jacobi preconditioning is applied on every GMRES path, including the
  second- and third-harmonic solves. Because the interaction matrix has
  `A_ii = alpha^-1`, the preconditioner is simply `diag(alpha)`. Note that it only
  improves convergence when polarisability varies between dipoles; for a single
  material it is a scalar multiple of the identity and has no effect.
- The O(N log N) FFT matrix-vector product can now be selected from a job file with
  `use_fcd = false`. On a 1331-dipole cube this measured 9.3x faster than the FCD
  path.

### Added

- `use_fcd` and `iterative_threshold` in the `[simulation]` section of a job file.
- `lumina_core::build_dipoles`, a single dipole-construction path shared by the
  command-line and graphical interfaces.
- `SceneDipoles::uniform_spacing()`.

### Fixed

- The workspace version had read 0.3.0 since before the 0.5 series and was written
  into the header of every exported data file, mislabelling results.
- None of the three bundled examples parsed, nor did the job files in the quick-start
  and geometry guides; all used a scene schema that predates the `shape` field.
- Integration tests built dipoles without the surface-averaging fill fraction, so the
  published accuracy figures described a configuration that production did not run.

### Changed

- Source, documentation and diagrams contain only ASCII characters.

### Known limitations

- Only right-handed helices can be constructed. There is no `handedness` field, and a
  negative pitch collapses the bounding box, so the circular-dichroism comparison in
  `tutorials/04-chiral-helix.md` cannot currently be followed.
- Periodic (Ewald) calculations have no entry point in either the command-line or
  graphical interface and are reachable only through the library API.
- The one-dimensional spectral sum asserts rather than degrading quietly when its
  mode cap binds, below a transverse separation of about 0.016 lattice periods.

## v0.6.0

Surface averaging. Recorded retrospectively; this work shipped without a version
entry.

- Per-dipole fill fractions from 5x5x5 sub-sampling of each cubic cell, carried
  through `LatticePoint` and `SceneDipoles`.
- Surface-averaged polarisability `alpha_SA = f * alpha_RRCM`, applied per dipole.
- Polarisability scaling is used rather than permittivity mixing, which produces a
  spurious Clausius-Mossotti resonance for Drude metals at fill fractions near 0.2
  to 0.3.

## v0.5.1 and earlier

See `docs/src/roadmap/versions.md` for the full history.
