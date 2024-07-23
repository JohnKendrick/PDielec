# ZnO VASP calculation of powder IR

1% mass fraction of ZnO is supported in a ptfe matrix.  Maxwell-Garnett and Bruggeman effective medium theories are used to examine the IR absorption of spheres, plates and needles.


| **File**                   | **Description**                                           |
| -------------------------- | --------------------------------------------------------- |
| zno_opt.cif                | Crystal structure from the VASP optimised calculation     |
| INCAR                      | VASP input settings                                       |
| KPOINTS                    | VASP k-point settings                                     |
| POSCAR                     | VASP cell and atom geometry                               |
| POTCAR.index               | A pointer to the pseudo potentials used in the calculation|
| OUTCAR                     | VASP output, analysed by PDielec                          |
| results.ref.xlsx           | Reference results for the test suite                      |
| script.py                  | pdgui script to create the reference results              |

