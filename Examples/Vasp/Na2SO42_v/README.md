# Na2(SO4)2 VASP calculation of powder IR

10% volume fraction of Na2(SO4)2 is supported in a ptfe matrix.  Maxwell-Garnett effective medium theory is used to examine the IR absorption of spheres in different support matrices.  In particular frequency dependent matrices are considered, as well as absorbing media.

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
| TestMaterialsDataBase.xlsx | Database containing permittivity of 'test' material       |

