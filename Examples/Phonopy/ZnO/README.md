# Examples/Phonopy/ZnO
## Phonopy-Vasp calculation on ZnO

This calculation was performed using the VASP interface.  A supercell of 1x1x1 was used for a Gamma point calculation
The input for the VASP calculations can be found in the VASP.Born directory
The calculation of the displaced energies and forces was prepared with;

```
phonopy -d --dim="1 1 1";
```

After the Vasp calculations the FORCE_SETS file was calculated with;

```
phonopy -f */vasprun.xml
```

The Born charges were calculated using Vasp on a non-supercell and the output is in the VASP.born/OUTCAR file. 
Generally speaking the KPOINTS file for the Born charge calculation will be different to that used
for the forces calculation.  There should be a KPOINTS file in the directory which reflects the KPOINTS 
used in the the supercell forces calculations. 

phonopy was used to calculate and store the dynamical matrix at the Gamma point with;

```
phonopy --dim="1 1 1" --qpoints-"0 0 0" --writedm
```

The output from this command is stored in qpoints.yaml and phonopy.yaml

| **File**          | **Description**                              |
| ----------------- | -------------------------------------------- |
| command.pdmake    | pdmake commands for test suite               |
| disp.yaml         | phonopy displacement info                    |
| phonopy_disp.yaml | phonopy displacement info                    |
| phonopy.yaml      | phonopy output file                          |
| qpoints.yaml      | phonopy qpoints file                         |
| FORCE_SETS        | phonopy Force sets file                      |
| KPOINTS           | VASP KPOINTS file                            |
| OUTCAR            | File linked to OUTCAR.born                   |
| OUTCAR.born       | VASP output file for Born charges            |
| results.ref.xlsx  | Reference results for test suite             |
| script.py         | pdgui test script                            |
| VASP.Born/        | Directory containing born charge calculation |

[Back](..)
