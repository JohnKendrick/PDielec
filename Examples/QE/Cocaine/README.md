# Cocaine using Quantum Espresso

| **File**        | **Description**                                                           |
| --------------- | ------------------------------------------------------------------------- |
| phonon.scf.in   |Input to a single point SCF at the optimised geometry                      |
| phonon.ph.in    |This is a phonon calculation using the wavefunction from the run above     |
| phonon.log      |This is the log of the output files from the SCF and Phonon calculations   |
| phonon.dynG     |This is the dynamical matrix from the phonon calculation                   |

phonon.dynG is used by PDielec to generate the absorption spectrum.

Warning: the PWSCF output indicates that the mass of Cl is 1.0.  
