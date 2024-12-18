# Examples/Vasp/Na2SO42
## Na2(SO4)2 VASP calculation of powder IR

10% volume fraction of Na2(SO4)2 is supported in a ptfe matrix.  Maxwell-Garnett and Bruggeman effective medium theories are used to examine the IR absorption of spheres, plates and needles.


| **File**                   | **Description**                                           |
| -------------------------- | --------------------------------------------------------- |
| Na2SO42_vasp_optimised.cif | Crystal structure from the VASP optimised calculation     |
| INCAR                      | VASP input settings                                       |
| KPOINTS                    | VASP k-point settings                                     |
| POSCAR                     | VASP cell and atom geometry                               |
| POTCAR.index               | A pointer to the pseudo potentials used in the calculation|
| OUTCAR                     | VASP output, analysed by PDielec                          |
| results.ref.xlsx           | Reference results for the test suite                      |
| script.py                  | pdgui script to create the reference results              |


[Back](..)
