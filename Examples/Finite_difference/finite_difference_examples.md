# Finite-difference JSON examples

These examples read a complete VASP finite-difference calculation of wurtzite
ZnO through the program name `"finite_field"`. This is the reader's existing
name for the consolidated `finite_difference.json` format; the Raman tensors
and phonons in this dataset come from Cartesian atomic displacements.

One shared, self-contained data file is stored at
[`../Powder_Raman/Finite_difference/finite_difference.json`](../Powder_Raman/Finite_difference/finite_difference.json).
No access to the original calculation directories is needed to run the examples.

## Tests

| Directory relative to `Examples/` | Coverage | Reference |
| --- | --- | --- |
| `Finite_difference/preader` | Infrared output with average masses and Eckart projection | `command.ref.csv` |
| `Finite_difference/preader_raman` | Total, parallel and perpendicular Raman activities with average masses and Eckart projection | `command.ref.csv` |
| `Powder_Raman/Finite_difference` | ZnO spheres in PTFE, particle frequency shifts and local fields, EO off/on | `results.ref.xlsx` |

From any of these three directories, run:

```sh
../../../pdmake command.pdmake
```

For the powder example on a machine without a display, use:

```sh
QT_QPA_PLATFORM=offscreen ../../../pdmake -cpus 2 command.pdmake
```

From the repository root, `./pdmake test-preader` includes both reader tests and
`./pdmake test-powder_raman` includes the powder test. The powder test is also
registered in `test-pdgui`, so the standard `./pdmake test` includes all three.
The quoted `"finite_field"` argument is supported in `.pdmake` recipes.

The powder spectrum covers 1–700 cm⁻¹ at 1 cm⁻¹ intervals, using VV polarisation,
a 785 nm laser, 298 K, 5 cm⁻¹ broadening and a crystal volume fraction of 0.1.
Both scenarios retain particle frequency and local-field corrections; only the
electro-optic tensor contribution is toggled. Spherical particles use analytical
orientation averaging.

Reference files are regression baselines, not independent proofs of the physics.
When the examples were added, both reader comparisons and the powder comparison
passed, and the two powder curves were finite, nonnegative, nonzero and distinct.
The complete pytest suite passed (1,104 tests). The wider `test-powder_raman`
run passed the new example, then stopped because the existing
`Powder_Raman/AbInit/results.ref.xlsx` was absent in this checkout.

## Calculation provenance

The JSON is a byte-for-byte copy, obtained on 2026-09-10, of:

```text
~/Leeds/ZnO/Validation_Protocols/VASP/LDA_package_equilibrium/M3_cartesian_raman/finite_difference.json
SHA256: 33445519cf12e1b9ca0037add16e8934d22870e2396716eabde1a20480acefcf
```

It contains a four-atom LDA ZnO cell, electronic permittivity, Born charges,
Cartesian dielectric derivatives, force constants, normal modes, Raman tensors
and a finite-electric-field chi2 tensor, with the original audit metadata.
The VASP Raman displacement calculations use a Gamma-centred **21×21×14** mesh,
600 eV cutoff, `EDIFF=1e-6` eV, `ISYM=0` and `LEPSILON=.TRUE.`. This was the
finest complete VASP Raman JSON dataset found in the source tree. The finer
24×24×16 calculations in `M1_response_convergence` establish dielectric/Born
convergence; they do not supply a complete Cartesian Raman dataset.

The selected Raman displacement is 0.010 Å; the selected force displacement is
0.015 Å. Forces come from `M5_raman_displacement_plateau`, and the reference
dielectric/Born response comes from `M2_phonon_convergence/potim_0p015`.
The dielectric derivatives have a translation-sum projection, and the symmetric
force constants have an acoustic-sum projection. Relevant source audits are
`M3_cartesian_raman/status/M3_cartesian_raman.md`,
`M5_raman_displacement_plateau/status/M5_raman_displacement_plateau.md`, and
`M1_response_convergence/status/M1_response_convergence.md`, relative to the
VASP `LDA_package_equilibrium` directory.

The chi2 tensor comes from `M4_finite_field_chi2`. Its retained validation status
is `complete_plateau_supported_with_minor_soft_gate_cautions`; the source export
accepted the complete field grid with remaining dielectric consistency cautions.
The EO example exercises that supplied tensor without upgrading its validation
status. All diagnostic fields are retained in the JSON.

## Normalisation and output conventions

The stored Raman tensor is `R_epsilon = sqrt(Vcell) d epsilon / dQ`, with units
`sqrt(Å/amu)`. The JSON also contains activity tensors
`A = sqrt(Vcell)/(4*pi) R_epsilon`, whose squared activities have units `Å⁴/amu`.
`preader` reports activities in `R_epsilon²` units (`Å/amu`), so converting them
to the JSON's `Å⁴/amu` convention requires the factor `Vcell/(16*pi²)`.

The reader recalculates modes and tensor projections using the selected masses.
These examples use PDielec average masses; the source uses Zn 65.38 and O 15.999
amu. Small frequency/activity differences from the source's tabulated modes are
therefore expected. With source masses restored, the reader's optical frequencies
and Raman tensors were checked against the JSON, allowing eigenvector sign
changes; chi2 was preserved within floating-point roundoff.

The top-level `chi2` is already in PDielec's internal convention. The separate
`chi2_physical_pm_per_v_c6v` field uses physical susceptibility in
`P = epsilon0 chi2 EE`, with `d = chi2/2`; it must not be substituted directly
for the internal tensor.

Generic `preader` CSV columns for k points, cutoff and energies contain reader
defaults because the finite-field reader does not populate those fields. Use
the provenance above for the actual VASP calculation settings.
