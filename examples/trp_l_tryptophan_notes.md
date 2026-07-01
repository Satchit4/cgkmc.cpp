# L-tryptophan example notes

`trp_l_tryptophan.json` is a literature-anchored coarse-grained starting
input for L-tryptophan (TRP), adapted to the current `cg++` core.

The current C++ implementation accepts orthorhombic periodic cells with
fractional basis coordinates, while the reported pure L-tryptophan structure is
triclinic P1. The JSON therefore uses a density-preserving orthorhombic
surrogate, not the full triclinic structure.

## Structural source

Gorbitz et al. report pure L-tryptophan as triclinic P1 at 123 K with:

- `a = 11.430 A`
- `b = 11.464 A`
- `c = 35.606 A`
- `alpha = 84.421 deg`
- `beta = 87.694 deg`
- `gamma = 60.102 deg`
- `V = 4025.6 A^3`
- `Z = 16`
- calculated density `1.348 g cm^-3`

Source: C. H. Gorbitz et al., "Single-crystal investigation of
L-tryptophan with Z' = 16", Acta Crystallographica Section B 68, 549-557
(2012), DOI: `10.1107/S0108768112033484`.

## Orthorhombic surrogate used here

The JSON keeps `a` and `b` from the literature and replaces `c` with:

```text
c_eff = V / (a * b) = 30.7219 A
```

This preserves the reported molecular volume:

```text
V / Z = 4025.6 / 16 = 251.6 A^3 per molecule
```

The 16 molecular sites in the cell are placed on a regular `2 x 2 x 4`
fractional grid. This preserves the number of molecules per unit cell and the
solid density but does not reproduce the true molecular centroids or triclinic
angles.

## Interaction energy source and coarse graining

The literature does not provide cgkmc shell energies for L-tryptophan. This
input uses a one-shell isotropic interaction chosen from the recommended
sublimation enthalpy implied by thermochemical data:

```text
Delta_sub H(298.15 K) ~= 167.5 kJ/mol = 1.736 eV/molecule
```

For a six-neighbor surrogate lattice:

```text
epsilon = 2 * E_coh / z = 2 * (-1.736 eV) / 6 = -0.578672 eV
```

This is a cohesive-energy normalization, not a fitted attachment-energy or
facet-specific model.

Thermochemistry source: V. P. Korolev et al., "Thermodynamic properties of
L-tryptophan", Journal of Chemical Thermodynamics 105, 44-49 (2017),
DOI: `10.1016/j.jct.2016.09.041`.

## Solvent and growth constants

- `beta = 38.9217 eV^-1` for 298.15 K.
- `diffusivity = 6.6e10 A^2/s`, converted from an estimated tryptophan
  diffusion coefficient of `6.6e-6 cm^2/s`.
- `solubility_limit = 3.3616e-5 A^-3`, converted from `11.4 g/L` at 25 C
  using molecular weight `204.2252 g/mol`.

Diffusion source: Lapidus et al., "Measuring the rate of intramolecular
contact formation in polypeptides", PNAS 97, 7220-7225 (2000),
DOI: `10.1073/pnas.97.13.7220`.

Solubility source: PubChem L-tryptophan compound record and the reviewed
proteinogenic amino-acid solubility literature.

## Intended use

Use this input to exercise the `cg++` workflow on a TRP-like molecular crystal.
Do not treat the morphology as validated until the triclinic geometry and
facet/shell interaction energies are parameterized explicitly.
