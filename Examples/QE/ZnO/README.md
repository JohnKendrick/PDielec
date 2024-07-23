# Quantum Espresso 5.1 Calculation on ZnO

| **File**          | **Description**                                                           |
| ----------------- | --------------------------------------------------------------------------|
| zno.scf.in        | Input to a single point SCF at the optimised geometry                     |
| zno.scf.out       | This is the log of the output files from the SCF                          |
| zno.ph.in         | This is a phonon calculation input                                        |
| zno.ph.out        | This is a phonon calculation output                                       |
| zno.dynG          | This is the dynamical matrix from the phonon calculation                  |
| command.pdmake    | pdmake command for test suite                                             |
| results.ref.xlsx  | Reference results for test suite                                          |
| script.py         | pdgui script to generate test results                                     |

zno.ph.dynG is used by PDielec to generate the absorption spectrum.
