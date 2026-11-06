# Examples/Powder_Raman/Particle_Shape

## ZnO particle-shape comparison (CASTEP 24.1)

The supplied `raman.castep` dataset uses PBE with the G06 dispersion correction.
`script.py` compares five particle scenarios using the same ZnO response:

| Scenario | Shape and orientation |
| --- | --- |
| Sphere | Spherical depolarisation tensor |
| Needle (001) | Needle axis along [001] |
| Plate(001) | Plate normal to (001) |
| Ellipse(001) | Ellipsoid with unique axis [001] and aspect ratio 1 |
| Plate(010) | Plate normal to (010) |

The ellipsoid's aspect ratio is **1**, so this scenario checks the sphere limit,
not an elongated ellipsoid. All scenarios use PTFE, a crystal volume fraction
of 0.1, a 785 nm laser, VV polarisation and 298 K. The EO term is disabled.
Particle frequency shifts and mode mixing remain enabled when the required
response data are present. All five shapes use analytical orientation averaging.

Run the saved comparison, including all five scenarios, from this directory:

```sh
../../../pdmake command.pdmake
../../../pdmake --view command.pdmake
```

To load only the underlying dataset, use
`../../../pdgui --spectroscopy powder_raman raman.castep`.

| File | Description |
| --- | --- |
| `raman.param` | CASTEP calculation parameters |
| `raman.cell` | CASTEP geometry |
| `raman.castep` | CASTEP output |
| `raman.efield` | Electric-field response |
| `raman.phonon` | Phonon response |
| `command.pdmake` | Example registered in `pdmake test-powder_raman` |
| `script.py` | Five saved powder Raman scenarios |
| `results.ref.xlsx` | Reference workbook for regression comparison |

The test writes `results.xlsx` and compares it with `results.ref.xlsx`.
The supplied reference currently contains only one `Scenario 1` spectrum over
1–200 cm⁻¹, whereas `script.py` produces the five scenarios above over
1–800 cm⁻¹ (0.2 cm⁻¹ spacing). It therefore needs a reviewed baseline update
before this new suite entry can pass.

Review the generated five-scenario workbook before using `../../../pdmake --regenerate command.pdmake`
to replace the reference.

[Back to Powder Raman](../README.md)
