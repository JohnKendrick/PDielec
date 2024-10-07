# Examples/QE/Urea
## Quantum Espresso 7.3.1 Powder & Single Crystal IR Calculation on Urea
The is example uses the pwscf.xml file to provide the information needed for pdgui.
The program ignores the log files completely.
The program uses pwscf.xml and urea.ir.dynG to provide the information it needs.
urea.ir.dynG provides the Born charges, the dynamical matrix and the optical permittivity,
pwscf.xml provides the rest.


| **File**          | **Description**                                                           |
| ----------------- | --------------------------------------------------------------------------|
| urea.scf.in       | Input to a single point SCF at the optimised geometry                     |
| urea.scf.log      | This is the log of the output files from the SCF                          |
| pwscf.xml         | The xml file created by the scf run                                       |
| urea.ir.in        | This is a phonon calculation input                                        |
| urea.ir.log       | This is a phonon calculation output                                       |
| urea.ir.dynG      | This is the dynamical matrix from the phonon calculation                  |
| command.pdmake    | pdmake command for test suite                                             |
| results.ref.xlsx  | Reference results for test suite                                          |
| script.py         | pdgui script to generate test results                                     |

urea.ph.dynG is used by PDielec to generate the absorption spectrum.

[Back](..)
