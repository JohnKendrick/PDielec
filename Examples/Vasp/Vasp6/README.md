# Examples/Vasp/Vasp6
## Na2(SO4)2 VASP6 calculation of powder IR

The VASP 6 calculation for this test was provided by David Santos-Carballal.
10% volume fraction of Na2(SO4)2 is supported in a ptfe matrix.  Maxwell-Garnett effective medium theory is used to examine the IR absorption.

| **File**                   | **Description**                                           |
| -------------------------- | --------------------------------------------------------- |
| INCAR                      | VASP input settings                                       |
| KPOINTS                    | VASP k-point settings                                     |
| POSCAR                     | VASP cell and atom geometry                               |
| OUTCAR                     | VASP output, analysed by PDielec                          |
| results.ref.xlsx           | Reference results for the test suite                      |
| script.py                  | pdgui script to create the reference results              |


[Back](..)
