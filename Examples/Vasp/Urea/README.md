# Examples/Vasp/Urea
## Urea VASP calculation using xml input

An "elastic" calculation of the permittivity, Born charges and frequencies.
The pdgui calculation in script.py include an powder and a crystal calculation
The input to pdgui comes entirely from vasprun.xml.


| **File**                   | **Description**                                           |
| -------------------------- | --------------------------------------------------------- |
| INCAR                      | VASP input settings                                       |
| KPOINTS                    | VASP k-point settings                                     |
| POSCAR                     | VASP cell and atom geometry                               |
| POTCAR.index               | A pointer to the pseudo potentials used in the calculation|
| OUTCAR.for_reference       | VASP output, kept only for reference purposes             |
| results.ref.xlsx           | Reference results for the test suite                      |
| script.py                  | pdgui script to create the reference results              |


[Back](..)
