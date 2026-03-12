# Nonlinear Optics (SHG/THG)

> **Status (v0.4.0):** Second-harmonic generation (SHG) and third-harmonic generation (THG) fully implemented for both CLI and GUI.

## Second-Harmonic Generation

For materials with broken inversion symmetry, the second-order nonlinear polarisation generates light at twice the fundamental frequency (SHG). Within the CDA framework, each dipole at position **r**_i generates a nonlinear source moment at \\(2\omega\\):

\\[
p^{\text{NL}}_{i,a}(2\omega) = \sum_{b,c} \chi^{(2)}_{abc}\,
  E_{\text{loc},b,i}(\omega)\, E_{\text{loc},c,i}(\omega)
\\]

where \\(\chi^{(2)}\\) is the second-order susceptibility tensor and \\(\mathbf{E}_{\text{loc}}(\omega)\\) is the local field from the linear CDA solve.

## Workflow

1. **Linear solve at ω** — obtain the dipole moments and local fields.
2. **Compute source terms** — contract \\(\chi^{(2)} : \mathbf{E}\mathbf{E}\\) per dipole.
3. **Assemble A(2ω)** — interaction matrix at the harmonic wavelength λ/2.
4. **Solve A(2ω) · p(2ω) = p_NL** — no incident field at the harmonic; the nonlinear sources drive the response.
5. **Compute SHG intensity** — \\(I_{\text{SHG}} \propto \sum_i |\mathbf{p}_i(2\omega)|^2\\).
6. **Optional** — compute far-field pattern at \\(2\omega\\).

## χ^(2) Symmetry Classes

### Isotropic surface (C_∞v)

For nanoparticles with surface-normal symmetry along **z**, the non-zero components are:

\\[
\chi_{zzz},\quad \chi_{zxx} = \chi_{zyy},\quad \chi_{xxz} = \chi_{xzx} = \chi_{yyz} = \chi_{yzy}
\\]

This gives exactly 7 non-zero elements. In code:

```rust
let chi2 = Chi2Tensor::isotropic_surface(
    Complex64::new(chi_zzz_re, chi_zzz_im),
    Complex64::new(chi_zxx_re, chi_zxx_im),
);
```

### Zero (centrosymmetric)

For structures with inversion symmetry, all \\(\chi^{(2)}\\) components are zero:

```rust
let chi2 = Chi2Tensor::zero();
```

## Scaling Law

The SHG intensity scales as the fourth power of the local field amplitude:

\\[
I_{\text{SHG}} \propto |\mathbf{E}_{\text{loc}}|^4
\\]

This is a key diagnostic: doubling the incident field should give 16× the SHG signal. This is verified by the `test_shg_intensity_scales_with_field` unit test.

## TOML Configuration (CLI)

```toml
[nonlinear]
enable_shg = true
symmetry = "isotropic_surface"   # or "zero"
chi_zzz = [1.0, 0.0]            # [real, imag] in nm³
chi_zxx = [0.3, 0.0]
far_field = false
```

Output: `output/shg_spectrum.csv` with columns `fundamental_nm`, `harmonic_nm`, `shg_intensity_nm6`.

## GUI

Enable SHG in the **Simulation** panel. When checked:

1. Select symmetry class from the dropdown.
2. Adjust \\(\chi_{zzz}\\) and \\(\chi_{zxx}\\) components with the drag inputs.
3. Click **Run Simulation** — the SHG spectrum appears as an additional tab in the results panel.

## Computational Cost

Each wavelength requires:
1. One linear solve at \\(\omega\\) (already needed for the linear spectrum).
2. One additional solve at \\(2\omega\\).

Total runtime is approximately \\(2\times\\) the linear-only cost.

## Third-Harmonic Generation

For centrosymmetric materials (bulk metals, glass), the first non-vanishing nonlinear term is χ^(3). The THG source moment at \\(3\omega\\) is:

\\[
p^{\text{NL}}_{i,a}(3\omega) = \sum_{b,c,d} \chi^{(3)}_{abcd}\,
  E_{\text{loc},b,i}(\omega)\, E_{\text{loc},c,i}(\omega)\, E_{\text{loc},d,i}(\omega)
\\]

The workflow is identical to SHG but solves at λ/3. The THG intensity scales as |E|⁶.

### χ^(3) Symmetry: Isotropic Bulk

For an isotropic centrosymmetric medium (Kleinman symmetry), the independent components are χ_xxxx and χ_xxyy. The full tensor is:

\\[
\chi^{(3)}_{abcd} = \chi_{\text{xxyy}} (\delta_{ab}\delta_{cd} + \delta_{ac}\delta_{bd} + \delta_{ad}\delta_{bc})
+ (\chi_{\text{xxxx}} - 3\chi_{\text{xxyy}}) \delta_{ab}\delta_{bc}\delta_{cd}
\\]

This gives χ_xxxx correctly (diagonal) and χ_xxyy = χ_xyxy = χ_xyyx (off-diagonal permutations). In code:

```rust
let chi3 = Chi3Tensor::isotropic_bulk(
    Complex64::new(chi3_xxxx_re, chi3_xxxx_im),
    Complex64::new(chi3_xxyy_re, chi3_xxyy_im),
);
```

### TOML Configuration (CLI)

```toml
[nonlinear]
enable_thg = true
chi3_symmetry = "isotropic_bulk"   # or "zero"
chi3_xxxx = [0.5, 0.0]            # [real, imag] in nm⁶
chi3_xxyy = [0.1, 0.0]
```

Output: `output/thg_spectrum.csv` with columns `fundamental_nm`, `harmonic_nm`, `thg_intensity_nm9`.

## Limitations

- χ^(2) and χ^(3) are uniform across all dipoles in a given object. Per-dipole or DFT-imported tensors are planned for v0.5.
- For SHG: only isotropic surface (C_∞v) and zero symmetry classes. For THG: only isotropic bulk and zero.
- The harmonic wavelength (λ/2 for SHG, λ/3 for THG) must lie within the material's tabulated range — points outside this range are silently skipped.
- Note: λ/3 is often in the UV (< 300 nm) and outside most material tables. Use a long-wavelength fundamental or provide a UV-range material.
