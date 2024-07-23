# Phonopy-VASP calculation of Na2(SO4)2

This calculation was performed using the VASP interface.  A supercell of 1x1x1 was used for a Gamma point calculation.
The input for the VASP calculations can be found in INCAR.supercell KPOINTS.supercell POSCAR and POTCAR-index in the VASP/ directory
The calculation of the displaced energies and forces was prepared with;

```
phonopy -d --dim="1 1 1";
```

After the Vasp calculations the FORCE_SETS file was calculated with;

```
phonopy -f */vasprun.xml
```

The Born charges were calculated using Vasp on a non-supercell and the output is in the VASP/OUTCAR file.  The appropriate KPOINTS and INCAR files are KPOINTS.born and INCAR.born also in the VASP/ directory
Generally speaking the KPOINTS file for the Born charge calculation will be different to that used for the forces calculation. 

phonopy was used to calculate and store the dynamical matrix at the Gamma point with;

```
phonopy --dim="1 1 1" --qpoints-"0 0 0" --writedm
```

The output from this command is stored in qpoints.yaml and phonopy.yaml

phonopy band.conf was used to generate the band.yaml file.  This contains the displacements of the phonon modes and can be visualised using;
phononwebsite at http://henriquemiranda.github.io/phononwebsite/

| **File**           | **Description**                                           |
| ------------------ | --------------------------------------------------------- |
| band.conf          | phonopy file generated running phonopy                    |
| band.yaml          | phonopy file generated running phonopy                    |
| phonopy_disp.yaml  | phonopy file generated running phonopy                    |
| phonopy.yaml       | phonopy file generated running phonopy                    |
| qpoints.yaml       | phonopy file generated running phonopy                    |
| disp.yaml          | phonopy file generated running phonopy                    |
| FORCE_SETS         | phonopy file generated running phonopy                    |
| KPOINTS            | The VASP KPOINTS file                                     |
| OUTCAR             | a link to OUTCAR.born                                     |
| OUTCAR.born        | VASP output containing born charges                       |
| SPOSCAR            | Supercell atom positions for VASP                         |
| VASP/              | Directory used to calculat Born charges (see OUTCAR.born) |
| command.pdmake     | pdmake commands for test suite                            |
| results.ref.xlsx   | Reference results for test suite                          |
| script.py          | pdmake commands for test suite                            |

