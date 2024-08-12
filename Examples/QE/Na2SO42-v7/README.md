# Examples/QE/Na2SO42-v7
## Quantum Espresso 7.3.1 Powder IR Calculation on Na2(SO4)2

| **File**          | **Description**                                                           |
| ----------------- | --------------------------------------------------------------------------|
| Na2SO42.scf.in    | Input to a single point SCF at the optimised geometry                     |
| Na2SO42.ph.in     | This is a phonon calculation using the wavefunction from the run above    |
| Na2SO42.log       | This is the log of the output files from the SCF and Phonon calculations  |
| Na2SO42.dynG      | This is the dynamical matrix from the phonon calculation                  |
| command.pdmake    | pdmake command for test suite                                             |
| results.ref.xlsx  | Reference results for test suite                                          |
| script.py         | pdgui script to generate test results                                     |

Na2SO42.dynG is used by PDielec to generate the absorption spectrum.

[Back](..)
