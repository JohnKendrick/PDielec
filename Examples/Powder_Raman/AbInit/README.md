# Examples/Powder_Raman/AbInit
## ZnO Raman Calculation using AbInit
This calculation uses the PEAD method.  The raman.abi input file produces several databases.
The databases for steps 4 and 5 were merged before processing with anaddb.

```
mrgddb analysis.ddb ramano_DS4_DDB ramano_DS5_DDB
anaddb < analysis.files
```

The script reads `raman.abo`. PDielec projects its atomic
susceptibility derivatives onto the phonon modes and converts them to the
shared dielectric-derivative convention `R_epsilon = sqrt(Vcell) dε/dQ`.
No manual conversion of the printed `analysis.abo` Raman susceptibilities is
needed; anaddb results are supplementary analysis.

The saved scenarios compare VV, VH and unpolarised scattering for ZnO spheres
in PTFE, plus an unpolarised scenario with the EO term enabled. All use
analytical orientation averaging. Legacy sample-count entries are ignored.

| **File**                   | **Description**                                       |
| ------------------------------- | ----------------------------------------------------- |
| opt.abi                         | AbInit optimisation input                             |
| opt.abo                         | AbInit optimisation output                            |
| raman.abi                       | AbInit Raman input                                   |
| raman.abo                       | AbInit Raman output                                   |
| analysis.abi                    | anaddb input for the analysis step                    |
| analysis.abo                    | anaddb output                                         |
| analysis.ddb                    | the merged database                                   |
| analysis.files                  | anaddb input                                          |
| command.pdmake                  | part of the pdmake tests                              |
| results.ref.xlsx               | reference workbook for regression comparison          |
| script.py                       | part of the pdmake test system                        |
| O.psp8,Zn.psp8                  | the pseudo-potentials used                            |

Run the example from this directory:

```sh
../../../pdmake command.pdmake
../../../pdmake --view command.pdmake
```

The first command compares `results.xlsx` with `results.ref.xlsx`; the second
opens the saved scenarios interactively. Loading only the raw input with
`pdgui` does not restore all settings in `script.py`.

[Back to Powder Raman](../README.md)

